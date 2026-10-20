"""Continuation search over 13, 26, 52 and 104 finite-S3 pulses.

First train *pointer-plus-system-sector capture*, which is insensitive to
within-syndrome rotations and avoids a misleading target-overlap plateau.
Then refine the actual coherent transfer on all 14 pointer-ready physical
directions. An independent PR68 certificate checks the resulting actual
42D pulse product. Neither stage inserts a free code-projection gate.

A larger depth always includes the preceding schedule padded by zero-angle
pulses. Selection maximizes the *worst* coherent branch transfer, then
dimension-weighted average, with identity as an explicit negative control.
This is numerical research, not an SU(3)/QCD or hardware decoder.
"""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize

from colour_s3_branch_aware_optimizer import (
    BRANCHES, BranchAwareCompiler,
)
from colour_s3_recovery_pulse_compiler import (
    TOTAL, native_operation_expansion,
)
from colour_s3_direct_logical_compiler import (
    no_jump_lifetime_comparison,
)


class SectorCaptureCompiler:
    """Subspace-capture warmup, with exact analytic reverse-mode gradient.

    q_b = ||T_b^dagger U(theta)X_b||_F^2 / rank_b counts the correct
    SYSTEM subspace and pointer simultaneously, WITHOUT postselection.

    This is only an optimization surrogate; high q_b may scramble the
    encoded logical state. Final selection uses *coherent* transfer and
    independent certification, not q_b.
    """

    def __init__(self, label="AS", cycles=1, min_capture=.80,
                 balance_weight=2., effort_lambda=.0002):
        for name, value in (("min_capture", min_capture),
                            ("balance_weight", balance_weight),
                            ("effort_lambda", effort_lambda)):
            if not np.isfinite(value):
                raise ValueError(name + " must be finite")
        if not 0 <= min_capture <= 1 or balance_weight < 0 or effort_lambda < 0:
            raise ValueError("invalid sector capture configuration")
        self.coherent = BranchAwareCompiler(
            label, cycles, effort_lambda=effort_lambda)
        self.base = self.coherent.base
        self.initial = self.coherent.initial
        self.target = self.coherent.target
        self.spans = self.coherent.spans
        self.n_angles = self.coherent.n_angles
        self.min_capture = float(min_capture)
        self.balance_weight = float(balance_weight)
        self.effort_lambda = float(effort_lambda)
        self.calls = 0

    def _capture_and_adjoint(self, output):
        if output.shape != (TOTAL, 14):
            raise ValueError("must propagate all 14 pointer-ready inputs")
        gradient_output = np.zeros_like(output)
        branch = {}
        mean_capture = 0.
        deficit = 0.
        for name, cols, rank in self.spans:
            target = self.target[:, cols]
            physical = output[:, cols]
            overlap = target.conj().T @ physical
            projected = target @ overlap
            capture = float(np.linalg.norm(overlap, "fro")**2 / rank)
            if not -1e-9 <= capture <= 1 + 1e-9:
                raise AssertionError("nonphysical syndrome-sector capture")
            weight = rank / 14.
            shortfall = max(0., self.min_capture - capture)
            mean_capture += weight * capture
            deficit += shortfall ** 2
            derivative_weight = weight + (2 * self.balance_weight / 3
                                           * shortfall)
            gradient_output[:, cols] = -derivative_weight * projected / rank
            branch[name] = capture
        loss = 1 - mean_capture + (self.balance_weight / 3) * deficit
        return float(loss), gradient_output, {
            "per_branch_correct_sector_and_pointer_probability": branch,
            "dimension_weighted_capture_probability": float(mean_capture),
            "worst_branch_capture_probability": float(min(branch.values())),
        }

    def objective_and_gradient(self, angles):
        self.calls += 1
        gates = self.base.matrices(angles)
        forward = [self.initial]
        for gate in gates:
            forward.append(gate @ forward[-1])
        loss, adjoint, _ = self._capture_and_adjoint(forward[-1])
        grad = np.zeros(self.n_angles)
        for k in range(self.n_angles - 1, -1, -1):
            derivative = -1j * (self.base.H[k] @ forward[k + 1])
            grad[k] = 2 * np.real(np.vdot(adjoint, derivative))
            adjoint = gates[k].conj().T @ adjoint
        effort, deffort = self.coherent._effort(angles)
        return (float(loss + self.effort_lambda * effort),
                grad + self.effort_lambda * deffort)

    def diagnostics(self, angles):
        output = self.base.propagate(angles, input_states=self.initial)
        loss, _, diagnostics = self._capture_and_adjoint(output)
        return dict(diagnostics, unregularized_capture_loss=loss)


def verify_sector_gradient(label="AS", cycles=1, seed=702,
                           step=1e-6):
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be finite positive")
    comp = SectorCaptureCompiler(label, cycles)
    angles = np.random.default_rng(seed).normal(0, .2, comp.n_angles)
    _, analytic = comp.objective_and_gradient(angles)
    rows = []
    for idx in sorted(set((0, 2, 7, comp.n_angles - 1))):
        d = np.zeros(comp.n_angles)
        d[idx] = step
        difference = (comp.objective_and_gradient(angles+d)[0] -
                      comp.objective_and_gradient(angles-d)[0])/(2*step)
        rows.append({
            "pulse_index":idx,
            "analytic_gradient":float(analytic[idx]),
            "finite_difference":float(difference),
            "absolute_error":float(abs(difference-analytic[idx])),
        })
    return {"label":label, "cycles":cycles, "sampled_derivatives":rows,
            "max_absolute_error":max(r["absolute_error"] for r in rows)}


def extend_angles(previous, next_size):
    """Appending identity-angle gates retains the EXACT earlier unitary."""
    previous = np.asarray(previous, dtype=float)
    if (previous.ndim != 1 or not 0 < len(previous) < next_size or
            not np.all(np.isfinite(previous)) or
            np.max(np.abs(previous)) > np.pi + 1e-12):
        raise ValueError("prior pulse vector must be shorter, finite, bounded")
    return np.pad(previous, (0, next_size-len(previous)))


def _branch_metrics(compiler, angles):
    """Selection strictly prefers worst-branch coherence, then total."""
    record = compiler.transfer_diagnostics(angles)
    scores = record["branch_coherent_transfer"]
    return (float(min(scores.values())),
            float(record["dimension_weighted_branch_transfer"]))


def _run_minimize(fun, initial, iterations):
    result = minimize(fun, initial, jac=True, method="L-BFGS-B",
                      bounds=[(-np.pi, np.pi)]*len(initial),
                      options={"maxiter":iterations, "maxls":20,
                               "gtol":1e-7, "ftol":1e-10})
    return np.asarray(result.x, float), {
        "solver_success":bool(result.success),
        "solver_message":str(result.message),
        "iterations":int(result.nit),
        "reported_objective":float(result.fun),
    }


def search_depth(label, cycles, previous=None, capture_iter=35,
                 coherent_iter=50, seed=20261019):
    """Search depth with two stages; return a fully audited real pulse product.

    Stochastic candidates are deterministic given seed. We retain all
    input candidates, including the previous smaller-depth pulse padded
    with identities. No solver can erase the previous worst-branch score.
    """
    if (type(capture_iter) is not int or capture_iter < 1 or
            type(coherent_iter) is not int or coherent_iter < 1):
        raise ValueError("positive integer iteration budgets required")
    capture = SectorCaptureCompiler(label, cycles)
    coherent = capture.coherent
    zero = np.zeros(coherent.n_angles)
    rng = np.random.default_rng(seed + 101*cycles + (label=="BS")*11)
    starts = [("zero_gates", zero)]
    if previous is not None:
        starts.append(("padded_previous_depth",
                       extend_angles(previous, coherent.n_angles)))
    starts.append(("random_exploration",
                   rng.normal(0, .30, coherent.n_angles)))
    if previous is not None:
        jitter = starts[1][1] + rng.normal(
            0, .10, coherent.n_angles)
        starts.append(("jittered_previous_depth",
                       np.clip(jitter, -np.pi, np.pi)))

    # Max-min score is a valid diagnostic selection criterion, NOT
    # a QEC theorem or a proof of convergence to the global optimum.
    candidates = []
    trajectory = []
    for name, init in starts:
        candidates.append((name + ":initial", init))
        sector_angles, stage1 = _run_minimize(
            capture.objective_and_gradient, init, capture_iter)
        candidates.append((name + ":sector", sector_angles))
        # Stage two protects the quantum state map, not just the pointer.
        coherent_angles, stage2 = _run_minimize(
            coherent.objective_and_gradient, sector_angles, coherent_iter)
        candidates.append((name + ":refined", coherent_angles))
        trajectory.append({"start":name, "sector_stage":stage1,
                           "coherent_stage":stage2})

    choice = None
    winner = None
    best_metrics = (-1., -1.)
    candidate_metrics = []
    for name, angles in candidates:
        m = _branch_metrics(coherent, angles)
        candidate_metrics.append({
            "candidate":name,
            "worst_branch_coherent_transfer":m[0],
            "dimension_weighted_coherent_transfer":m[1],
        })
        if m > best_metrics:
            choice = angles.copy()
            winner = name
            best_metrics = m
    assert choice is not None
    transfer, cert = coherent.independently_audit(choice)
    sector = capture.diagnostics(choice)
    for name in BRANCHES:
        from_certificate = cert["branches"][name][
            "intended_pointer_and_system_sector_probability"]
        if abs(sector["per_branch_correct_sector_and_pointer_probability"][name]
               - from_certificate) > 2e-9:
            raise AssertionError("independent complete-sector audit mismatch")
        if transfer["branch_coherent_transfer"][name] > from_certificate+2e-9:
            raise AssertionError("coherent transfer exceeds captured sector")

    native = native_operation_expansion({
        "pulse_names":coherent.names,
        "angles_rad":choice.tolist(),
        "pulse_count":coherent.n_angles,
    }, readout_ns=300. if np.any(np.abs(choice)>1e-9) else 0.)
    logical = cert["independent_seven_Kraus_single_fault_channel"][
        "unconditional_seven_fault_mixture_logical_F_e"]
    lifetime = no_jump_lifetime_comparison(
        logical, .63+.37/9,
        native["hypothetical_serial_duration_ns_including_readout"])
    return {
        "group":label, "cycles":cycles,
        "pulse_count":coherent.n_angles,
        "pulse_names":coherent.names,
        "angles_rad":choice.tolist(),
        "selected_candidate":winner,
        "candidate_metrics":candidate_metrics,
        "optimization_trajectory":trajectory,
        "capture_only_diagnostics":sector,
        "coherent_branch_diagnostics":transfer,
        "worst_coherent_branch_transfer":best_metrics[0],
        "mean_coherent_branch_transfer":best_metrics[1],
        "independent_instrument_certificate":cert,
        "three_outcome_instrument_certified":cert[
            "mathematical_three_outcome_instrument_pass"],
        "seven_fault_projected_logical_F_e_diagnostic_only":float(logical),
        "no_gate_projected_F_e":.63+.37/9,
        "assumed_native_costs":native,
        "optimistic_no_jump_lifetime_sensitivity":lifetime,
        "all_three_branches_audited_without_postselection":True,
        "full_physical_feedback_and_reset_implemented":False,
        "no_ideal_code_projection_inserted_as_pulse":True,
    }


def run_continuation(depths=(1, 2, 4, 8), labels=("AS", "BS"),
                     capture_iter=35, coherent_iter=50, seed=20261019):
    if not depths or any(type(n) is not int or n < 1 for n in depths):
        raise ValueError("depths must be a nonempty sequence of positive ints")
    if any(a >= b for a, b in zip(depths, depths[1:])):
        raise ValueError("depths must be strictly increasing")
    if any(x not in ("AS", "BS") for x in labels) or not labels:
        raise ValueError("labels must be a nonempty subset of AS, BS")
    rows = []
    for label in labels:
        previous = None
        previous_worst = -1.
        for depth in depths:
            r = search_depth(label, depth, previous,
                             capture_iter=capture_iter,
                             coherent_iter=coherent_iter, seed=seed)
            if previous is not None:
                if r["worst_coherent_branch_transfer"] < previous_worst-2e-9:
                    raise AssertionError("deeper warmstarted search regressed")
            previous_worst = r["worst_coherent_branch_transfer"]
            previous = r["angles_rad"]
            rows.append(r)
    return {
        "experiment":"two-stage sector capture then coherent instrument search",
        "depth_cycles":list(depths),
        "physical_model":"14D finite-S3 Gauss sector x neutral qutrit pointer",
        "all_fourteen_directions":True,
        "independent_certificate_threshold":.99,
        "gradient_validation":verify_sector_gradient(),
        "results":rows,
        "claim_limits":[
            "Only short fixed-order numerical pulse products were searched.",
            "Subspace capture does not imply preserved quantum coherence.",
            "The projected seven-fault F_e diagnostic retains an ideal final code projection.",
            "Device-specific CP22, T1/T2, fourth-level leakage, pointer readout/reset and feedback are not simulated.",
            "No claim of continuous SU(3), nonperturbative QCD or a functioning hardware code.",
        ],
    }


def main():
    data = run_continuation()
    out = Path("colour_s3_continuation_results.json")
    out.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n")
    for row in data["results"]:
        print("CONTINUATION", json.dumps({
            "group":row["group"],"cycles":row["cycles"],
            "pulses":row["pulse_count"],
            "worst_branch":row["worst_coherent_branch_transfer"],
            "mean_branch":row["mean_coherent_branch_transfer"],
            "sector_capture":row["capture_only_diagnostics"][
                "per_branch_correct_sector_and_pointer_probability"],
            "passed":row["three_outcome_instrument_certified"],
            "projected_single_fault_Fe_diagnostic_only":
                row["seven_fault_projected_logical_F_e_diagnostic_only"],
            "hypothetical_native_count":row["assumed_native_costs"][
                "two_register_native_or_abstract_primitive_count"],
            "hypothetical_serial_ns":row["assumed_native_costs"][
                "hypothetical_serial_duration_ns_including_readout"],
        },sort_keys=True),flush=True)


if __name__ == "__main__":
    main()
