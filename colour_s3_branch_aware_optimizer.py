"""Branch-aware local-control search for the complete finite-S3 instrument.

Unlike PR67's seven-fault logical fidelity objective, this optimizer uses
ALL FOURTEEN Gauss-physical pointer-ready input directions, and scores
code -> code/pointer0, electric -> code/pointer1, complement -> flag/pointer2.
It maximizes branch-coherent target transfer, allowing a separate phase
per measured pointer branch, and penalizes the weakest branches.

No ideal final code projection is applied to the evolving states in the
objective. PR67's projected seven-Kraus logical fidelity is reported
separately and clearly remains optimistic until a native decoder exists.
Every optimized schedule is checked independently with PR68's certificate.

This is a finite-S3 model, NOT continuous SU(3) or a device-ready decoder.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from colour_s3_recovery_pulse_compiler import (
    POINTER, TOTAL, SpectralCompiler, native_operation_expansion,
)
from colour_s3_instrument_certificate import certify_three_outcome
from colour_s3_direct_logical_compiler import (
    effort_weight_for_name, no_jump_lifetime_comparison,
)

BRANCHES = ("code", "electric", "flag")


class BranchAwareCompiler:
    """Optimize three target transfers with an analytic adjoint gradient.

    Score = sum_b (dim_b/14) * (1-F_b)
          + branch_penalty/3 * sum_b max(0, min_transfer-F_b)^2
          + effort_lambda * sum_k cost_k*(sqrt(theta_k^2+eps^2)-eps).

    F_b = |Tr(T_b^dagger U(theta) X_b)/dim_b|^2.
    F_b <= intended pointer probability: no off-pointer or off-code
    amplitude can masquerade as successful syndrome extraction.

    Branch penalties are search heuristics, not strict constraints.
    Only PR68's independent 0.99 certificate determines acceptance.
    """

    def __init__(self, label="AS", cycles=1, min_transfer=.90,
                 branch_penalty=2., effort_lambda=.0005,
                 effort_epsilon=.02):
        for name, value in (
            ("min_transfer", min_transfer),
            ("branch_penalty", branch_penalty),
            ("effort_lambda", effort_lambda),
            ("effort_epsilon", effort_epsilon),
        ):
            if not np.isfinite(value):
                raise ValueError(name + " must be finite")
        if not 0 <= min_transfer <= 1:
            raise ValueError("min_transfer must lie in [0,1]")
        if branch_penalty < 0 or effort_lambda < 0 or effort_epsilon <= 0:
            raise ValueError("penalty must be nonnegative and epsilon positive")
        self.base = SpectralCompiler(label, cycles)
        self.label = label
        self.cycles = cycles
        self.min_transfer = float(min_transfer)
        self.branch_penalty = float(branch_penalty)
        self.effort_lambda = float(effort_lambda)
        self.effort_epsilon = float(effort_epsilon)
        self.pointer0 = np.eye(POINTER)[:, [0]]
        self.spans = []
        sources = []
        for name in BRANCHES:
            basis = self.base.branches[name]
            start = sum(x.shape[1] for x in sources)
            source = np.kron(basis, self.pointer0)
            self.spans.append((name, slice(start, start + basis.shape[1]),
                               basis.shape[1]))
            sources.append(source)
        self.initial = np.column_stack(sources)
        self.target = self.base.target_unitary @ self.initial
        assert self.initial.shape == (TOTAL, 14)
        self.order = self.base.order
        self.n_angles = len(self.order)
        self.names = [self.base.named[i][0] for i in self.order]
        self.effort_weights = np.array(
            [effort_weight_for_name(n) for n in self.names], dtype=float)
        self.calls = 0

    def _terminal(self, achieved):
        """Independent per-pointer-branch overlaps and score adjoint."""
        adjoint = np.zeros_like(achieved)
        scores = {}
        deficit = 0.
        weighted_transfer = 0.
        for name, cols, rank in self.spans:
            y = achieved[:, cols]
            t = self.target[:, cols]
            amplitude = np.vdot(t, y)
            f = float(abs(amplitude / rank) ** 2)
            # Numerical tolerance only; no silent normalization of failures.
            if not -1e-9 <= f <= 1 + 1e-9:
                raise AssertionError("invalid branch coherent transfer")
            scores[name] = f
            weight = rank / 14.
            weighted_transfer += weight * f
            shortfall = max(0., self.min_transfer - f)
            deficit += shortfall ** 2
            coefficient = weight + 2 * self.branch_penalty / 3 * shortfall
            adjoint[:, cols] = coefficient * amplitude / rank ** 2 * t
        loss = (1 - weighted_transfer
                + self.branch_penalty / 3 * deficit)
        return float(loss), adjoint, {
            "branch_coherent_transfer": scores,
            "dimension_weighted_branch_transfer": float(weighted_transfer),
            "squared_shortfall_penalty": float(
                self.branch_penalty / 3 * deficit),
            "worst_branch_transfer": float(min(scores.values())),
        }

    def _effort(self, angles):
        a = np.asarray(angles, dtype=float)
        root = np.sqrt(a*a + self.effort_epsilon**2)
        amount = float(np.dot(
            self.effort_weights, root - self.effort_epsilon))
        deriv = self.effort_weights * a / root
        return amount, deriv

    def objective_and_gradient(self, angles):
        self.calls += 1
        gates = self.base.matrices(angles)
        forward = [self.initial]
        for gate in gates:
            forward.append(gate @ forward[-1])
        loss, back, _ = self._terminal(forward[-1])
        gradient = np.zeros(self.n_angles)
        for k in range(self.n_angles - 1, -1, -1):
            dy = -1j * (self.base.H[k] @ forward[k + 1])
            gradient[k] = -2 * np.real(np.vdot(back, dy))
            back = gates[k].conj().T @ back
        effort, effort_grad = self._effort(angles)
        return (float(loss + self.effort_lambda * effort),
                gradient + self.effort_lambda * effort_grad)

    def transfer_diagnostics(self, angles):
        final = self.base.propagate(angles, input_states=self.initial)
        loss, _, records = self._terminal(final)
        effort, _ = self._effort(angles)
        return dict(records,
                    instrument_only_loss=loss,
                    smooth_effort_proxy=effort,
                    regularized_objective=loss + self.effort_lambda * effort)

    def independently_audit(self, angles, threshold=.99):
        """Certify actual pulse product; do not swap in ideal U42."""
        full_U = self.base.propagate(
            angles, input_states=np.eye(TOTAL, dtype=complex))
        certificate = certify_three_outcome(
            full_U, self.label, pointer_threshold=threshold,
            transfer_threshold=threshold)
        direct = self.transfer_diagnostics(angles)
        for name in BRANCHES:
            certified = certificate["branches"][name][
                "branch_coherent_transfer_fidelity"]
            if abs(direct["branch_coherent_transfer"][name] - certified)>2e-9:
                raise AssertionError(
                    "optimizer differs from independent branch certificate")
        return direct, certificate


def finite_gradient_audit(label="AS", cycles=1, seed=812,
                          epsilon=1e-6):
    compiler = BranchAwareCompiler(label, cycles)
    angles = np.random.default_rng(seed).normal(
        0, .2, compiler.n_angles)
    _, grad = compiler.objective_and_gradient(angles)
    rows = []
    for k in sorted(set((0, 1, 4, compiler.n_angles-1))):
        delta = np.zeros(compiler.n_angles)
        delta[k] = epsilon
        fd = (compiler.objective_and_gradient(angles+delta)[0] -
              compiler.objective_and_gradient(angles-delta)[0]
              )/(2*epsilon)
        rows.append({"index": k, "adjoint": float(grad[k]),
                     "finite_difference": float(fd),
                     "absolute_error": float(abs(fd-grad[k]))})
    return {"group":label,"cycles":cycles,"checks":rows,
            "max_absolute_error":max(x["absolute_error"] for x in rows)}


def compile_branches(label="AS", cycles=1, maxiter=35, seed=20261019,
                     min_transfer=.90, branch_penalty=2.,
                     effort_lambda=.0005, warmstart=None):
    """Search bounded native-control angles, retain zero-angle policy.

    A lower loss does NOT imply a passing certificate. In particular the
    no-gate control is never mislabeled as a three-outcome measurement.
    """
    if type(maxiter) is not int or maxiter < 1:
        raise ValueError("maxiter must be a positive integer")
    comp = BranchAwareCompiler(label, cycles, min_transfer,
                               branch_penalty, effort_lambda)
    zero = np.zeros(comp.n_angles)
    random = np.random.default_rng(seed+cycles+97*(label=="BS"))
    starts = [("identity", zero),
              ("seeded_small", random.normal(0,.25, comp.n_angles)),
              ("seeded_medium", random.normal(0,.60, comp.n_angles))]
    if warmstart is not None:
        prior = np.asarray(warmstart, dtype=float)
        if (prior.ndim != 1 or len(prior)>=comp.n_angles
                or not np.all(np.isfinite(prior))):
            raise ValueError("warmstart must be a shorter finite schedule")
        starts.append(("previous_depth", np.pad(
            prior, (0, comp.n_angles-len(prior)))))
    baseline = comp.objective_and_gradient(zero)[0]
    best = (baseline, zero.copy(), "no_gates")
    attempts = []
    for name, init in starts:
        result = minimize(comp.objective_and_gradient, init,
                          method="L-BFGS-B", jac=True,
                          bounds=[(-np.pi, np.pi)]*comp.n_angles,
                          options={"maxiter":maxiter, "maxls":20,
                                   "gtol":1e-7, "ftol":1e-10})
        value = comp.objective_and_gradient(result.x)[0]
        attempts.append({
            "initialization":name, "objective":float(value),
            "iterations":int(result.nit),
            "solver_success":bool(result.success),
            "solver_message":str(result.message),
        })
        if np.isfinite(value) and value < best[0] - 1e-11:
            best = (value, result.x.copy(), name)
    loss, angles, winner = best
    transfer, certificate = comp.independently_audit(angles)
    zero_transfer = comp.transfer_diagnostics(zero)
    native = native_operation_expansion({
        "pulse_names":comp.names, "angles_rad":angles.tolist(),
        "pulse_count":comp.n_angles,
    }, readout_ns=(300. if np.any(np.abs(angles)>1e-9) else 0.))
    logical = certificate["independent_seven_Kraus_single_fault_channel"][
        "unconditional_seven_fault_mixture_logical_F_e"]
    identity_logical = .63 + .37 / 9
    time = native["hypothetical_serial_duration_ns_including_readout"]
    lifetime = no_jump_lifetime_comparison(
        logical, identity_logical, time)
    return {
        "group":label, "cycles":cycles,
        "pulse_count":comp.n_angles, "angles_rad":angles.tolist(),
        "pulse_names":comp.names,
        "winning_start":winner,
        "initial_no_gate_objective":float(baseline),
        "selected_objective":float(loss),
        "improved_branch_search_objective_over_no_gate":bool(
            loss < baseline-1e-10),
        "zero_gate_branch_transfer":zero_transfer,
        "actual_branch_transfer":transfer,
        "complete_instrument_certificate":certificate,
        "all_three_branch_acceptance_passed":certificate[
            "mathematical_three_outcome_instrument_pass"],
        "attempts":attempts,
        "objective_evaluations":comp.calls,
        "native_costs_hypothetical":native,
        "optimistic_seven_fault_logical_Fe_with_ideal_code_projection":
            float(logical),
        "lifetime_sensitivity_NOT_physical_fidelity":lifetime,
        "no_ideal_projector_applied_as_control_operation":True,
        "no_physical_device_model_or_hardware_success_claim":True,
        "no_uncompiled_feedback_gate_counted":True,
    }


def main():
    results = {
        "experiment":"full finite-S3 14D three-branch-aware pulse search",
        "method":"dimension-weighted branch transfer + shortfall penalty",
        "minimum_certificate_threshold":.99,
        "gradient_audit":finite_gradient_audit(),
        "runs":[compile_branches(label, cycles, maxiter=35)
                for label in ("AS", "BS") for cycles in (1, 2)],
        "limitations":[
            "A numerical pulse-product 42D unitary, not a device compilation.",
            "The 0.99 certificate is a target threshold, not an error threshold.",
            "Projected seven-fault F_e is diagnostic and assumes an ideal final projection.",
            "Timing, charged transport, readout and T1 inputs are hypothetical.",
            "No pointer reset, feedback latency, fourth-level leakage or device CP22 included.",
            "Finite S3 is not continuous SU3 or QCD.",
        ],
    }
    Path("colour_s3_branch_aware_results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n")
    for r in results["runs"]:
        print(json.dumps({
            "group":r["group"],"cycles":r["cycles"],
            "best_loss":r["selected_objective"],
            "branch_transfer":r["actual_branch_transfer"][
                "branch_coherent_transfer"],
            "passed_certificate":r["all_three_branch_acceptance_passed"],
            "optimistic_projected_logical_Fe":
                r["optimistic_seven_fault_logical_Fe_with_ideal_code_projection"],
            "hypothetical_serial_duration_ns":
                r["native_costs_hypothetical"][
                    "hypothetical_serial_duration_ns_including_readout"],
        }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
