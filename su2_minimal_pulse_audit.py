"""Decision test: one smooth SU(2) gate versus long composite controls.

A single SU(2)-invariant smooth two-body exchange pulse uses no physical
local X/Y echo or SK1 correction. A composite suppresses cross-register
exchange drift but needs many additional, precisely calibrated local
SU(2)-invariant rotations. Compare *all physically required controls*.

Methods:
- Exact time-ordered evolution of the eight-spin, SU(2)-singlet sector
  (14D block of the full 256D microscopic Hamiltonian).
- A reproducible nominal single pulse calibrated to a CZ-class gate
  modulo fixed, nominally tracked local/global Z phases.
- Single and composite faulty 4x4 projected operations compared with
  *their own nominal operation*, no fault-dependent logical rephasing.
- Separate +/-1%, +/-5%, +/-10% coherent error on either physical
  cross coupling, both cross couplings, and (composite only) all required
  intra-block rotations. Max final product-singlet-code leakage measured.
- A small, bounded search across alternate positive smooth profiles
  (not a general waveform optimization or a global no-go theorem).
- Exact symbolic weak-coupling obstruction: nonzero second-order ZZ
  phase proportional to J1 J2 has a nonzero first-order derivative
  with respect to either calibration error, for any common positive
  fixed-shape envelope.

Outputs deterministic, dimensionless operator-distance diagnostics,
NOT average gate fidelities, stochastic noise thresholds, hardware
time-to-error claims or full Standard Model quantum error correction.
"""
import numpy as np
import sympy as sp
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

from su2_shaped_pulse import metric
from su2_sk1_composite_pulse import PhysicalComposite, aligned_error

NOMINAL_SINGLE_DURATION = 38.234506218
TEST_BINS = 256
VALIDATE_BINS = 512
LEVELS = (0.01, 0.05, 0.10)


def weak_coupling_obstruction():
    """Scope: second-order Schrieffer-Wolff, fixed smooth envelope.

    If J1=a1*f1(t), J2=a2*f2(t) with a fixed positive/negative
    waveform shape, leading ZZ phase is -const*a1*a2.
    Any independent fractional rescaling causes first-order phase
    error whenever the intended entangling phase is nonzero.
    It does NOT rule out nonperturbative sweet spots, signed
    composite controls, additional axes or measurement/recovery.
    """
    a, b, scale = sp.symbols("a b C", positive=True, real=True)
    eps_a, eps_b = sp.symbols("eps_a eps_b", real=True)
    theta = scale*a*b
    disturbed = scale*a*b*(1+eps_a)*(1+eps_b)
    assert sp.diff(disturbed, eps_a).subs({eps_a:0,eps_b:0}) == theta
    assert sp.diff(disturbed, eps_b).subs({eps_a:0,eps_b:0}) == theta
    assert sp.simplify(disturbed-theta-theta*(
        eps_a + eps_b + eps_a*eps_b)) == 0
    print("Second-order fixed-envelope ZZ phase: each coupler has "
          "nonzero first derivative when entangling phase is nonzero.")


def one_pulse(model, duration, a=0.0, fault1=0.0,
              fault2=0.0, bins=TEST_BINS):
    """One real time-dependent physical exchange waveform, no local gates.

    Modulation parameter a is bounded to [-0.25,0.25]. The couplings
    are nonnegative and peak-normalized, so this is a *small sampled
    alternative* to sin^2, not a change to the physical allowed
    interaction operators. Coherent errors rescale the two couplers
    individually. Each fixed control shape is used both nominally
    and when faulty.
    """
    assert abs(a)<=0.25
    if a==0:
        return model.smooth(duration,fault1,fault2,bins)
    U=np.eye(14,dtype=complex)
    for i in range(bins):
        u=(i+0.5)/bins
        f=(np.sin(np.pi*u)**2)*(1+a*np.cos(2*np.pi*u))/(1+abs(a))
        energies, vec=eigh(model.h0+f*((1+fault1)*model.h1+
                                        (1+fault2)*model.h2),
                           check_finite=False)
        U=(vec*np.exp(-1j*energies*duration/bins)) @ vec.T @ U
    return U


def fault_cases(pct):
    return {
        "first+": (pct,0), "first-": (-pct,0),
        "second+": (0,pct), "second-": (0,-pct),
        "both+": (pct,pct), "both-": (-pct,-pct),
    }


def compare_controls(model):
    t,c=model.calibrate()
    composite_duration=2*(t+2*c)
    short=one_pulse(model,NOMINAL_SINGLE_DURATION,bins=VALIDATE_BINS)
    comp=model.composite(t,c,bins=VALIDATE_BINS)
    short_nom=metric(short,model.code)
    comp_nom=metric(comp,model.code)
    assert short_nom["deviation"]<1e-4
    assert short_nom["leakage"]<1e-6
    assert comp_nom["deviation"]<1e-4
    assert comp_nom["leakage"]<5e-6
    assert composite_duration > 8*NOMINAL_SINGLE_DURATION

    observed={}
    for e in LEVELS:
        short_cross={}
        comp_cross={}
        for name,(e1,e2) in fault_cases(e).items():
            short_cross[name]=aligned_error(
                model.code,
                one_pulse(model,NOMINAL_SINGLE_DURATION,
                          fault1=e1,fault2=e2,bins=TEST_BINS),
                one_pulse(model,NOMINAL_SINGLE_DURATION,bins=TEST_BINS))
            comp_cross[name]=aligned_error(
                model.code,
                model.composite(t,c,e1,e2,bins=TEST_BINS),
                model.composite(t,c,bins=TEST_BINS))
        local_errors={}
        for sign in (-1,1):
            local_errors["local%+d"%(sign,)]=aligned_error(
                model.code,
                model.composite(t,c,bins=TEST_BINS,
                                local_error=sign*e),
                model.composite(t,c,bins=TEST_BINS))
        short_worst=max(v[0] for v in short_cross.values())
        comp_worst_cross=max(v[0] for v in comp_cross.values())
        comp_worst_local=max(v[0] for v in local_errors.values())
        comp_worst_all=max(comp_worst_local,comp_worst_cross)
        observed[e]=(short_worst,comp_worst_cross,comp_worst_local,
                     comp_worst_all)
        print(f"Calibration {e:.0%}: short worst six cross scenarios "
              f"{short_worst:.8f}; composite worst cross "
              f"{comp_worst_cross:.8f}; composite worst extra-local "
              f"{comp_worst_local:.8f}; composite max sampled "
              f"{comp_worst_all:.8f}")
        print("  Short cross-coupler cases:",short_cross)
        print("  Composite local-rotation cases:",local_errors)
        assert 0<short_worst<1
        assert comp_worst_cross<short_worst
        assert comp_worst_local>short_worst
        assert comp_worst_all>short_worst
    print("Nominal single CZ-class pulse:",short_nom,
          "T Delta=",NOMINAL_SINGLE_DURATION,
          "zero EXTRA local exchange rotations")
    print("Nominal composite:",comp_nom,
          "cross-exchange T Delta=",composite_duration,
          "PLUS uncounted local rotation durations")
    print("For declared independent cross+local calibration scenarios, "
          "single physical pulse has smaller worst sampled deviation.")
    return observed


def limited_shape_search(model):
    """Compare three positive single-waveform candidates at bounded peak.

    Calibrate each duration to its own nominal CZ-class target, then
    test fixed-calibration +/-5% errors in each cross coupler and
    both together. Searches are local and cannot prove a global
    optimum nor absence of nonperturbative sweet spots.
    """
    outcome={}
    for a in (-0.2, 0.0, 0.2):
        if a==0:
            duration=NOMINAL_SINGLE_DURATION
        else:
            objective=lambda t: metric(
                one_pulse(model,t,a,bins=TEST_BINS),model.code)["deviation"]
            opt=minimize_scalar(objective,bounds=(30.0,75.0),
                                method="bounded",
                                options={"xatol":1e-5,"maxiter":100})
            assert opt.success
            duration=float(opt.x)
        nominal=one_pulse(model,duration,a,bins=VALIDATE_BINS)
        values=metric(nominal,model.code)
        cases={}
        for name,(e1,e2) in fault_cases(.05).items():
            dev,leak=aligned_error(
                model.code,
                one_pulse(model,duration,a,e1,e2,bins=VALIDATE_BINS),
                nominal)
            cases[name]=(dev,leak)
        worst=max(x[0] for x in cases.values())
        print(f"Bounded profile a={a:+.1f}, duration={duration:.6f}, "
              f"nominal deviation={values['deviation']:.8g}, "
              f"leakage={values['leakage']:.8g}, "
              f"worst +/-5% coupler deviation={worst:.8f}")
        outcome[a]=(duration,values,worst)
        assert np.isfinite(worst)
        assert values["leakage"]<.05
    print("Three-point waveform scan only; no claim of a universal "
          "no-go or globally best simple waveform.")
    return outcome


def test():
    weak_coupling_obstruction()
    model=PhysicalComposite()
    compare_controls(model)
    limited_shape_search(model)
    print("ALL MINIMAL-PULSE SU2 ROBUSTNESS AUDITS PASSED")


if __name__=="__main__":
    test()
