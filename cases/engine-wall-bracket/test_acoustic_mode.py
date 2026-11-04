"""tests for the shared acoustic mode's fit and removal.

the fit has to recover a mode it is given, and the removal has to take
the mode out completely, including when the mode carries a phase.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import acoustic_mode  # noqa: E402


def synth(period, amplitude, phase, constant, drift, n=80, start=0):
    """a leg's series carrying a known mode."""
    points = []
    for index in range(n):
        step = start + index * 50
        omega = 2.0 * math.pi / period
        value = (
            constant
            + drift * step
            + amplitude * math.cos(omega * step + phase)
        )
        points.append((step, value))
    return points


def test_a_pure_drift_is_not_called_a_shared_mode():
    """a monotone climb is not a common oscillation.

    the legs march from the same start and the chambers develop the same
    way, so a drift LOOKS common across the legs: the same timeline
    produces the same climb in each. the sharing test must reject it.
    a true mode arrives in step and stays in step, while a drift
    leaves the phase scattered and leaves one leg unable to predict
    another.
    """
    steps = [250 * index for index in range(12)]
    series = {}
    for offset, leg in enumerate(("a", "b", "c")):
        # each leg climbs at its OWN rate and carries its OWN wobble, so
        # the legs are not scalar multiples of one another. legs that
        # were proportional would be trivially predictable from one
        # another, which is not what a drift looks like in practice.
        wobble = [0.001, -0.002, 0.0005][offset % 3]
        series[leg] = [
            (step,
             1.0 + 0.002 * (1.0 + 0.3 * offset) * step
             + wobble * ((index * (offset + 1)) % 4 - 1.5))
            for index, step in enumerate(steps)
        ]
    fit = acoustic_mode.fit_shared_mode(series)
    sharing = acoustic_mode.mode_sharing(fit, series)
    assert not sharing["shares"], sharing


def test_a_period_shorter_than_the_sampling_limit_is_refused():
    """the fit must not select a period the samples cannot resolve.

    a sinusoid sampled every DELTA steps is indistinguishable from one
    sampled at any period below two DELTA steps. if the search were
    allowed there it could select a harmonic of the true mode purely
    because of the sampling grid, and the fit would look confident while
    being an artefact.
    """
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, 0.0, 12.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    spacing = 25.0
    assert fit["period_steps"] >= 2.0 * spacing, fit


def test_the_fit_recovers_an_injected_period():
    """the joint fit finds the period it was given."""
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, 0.0, 12.0, 0.001),
        "c": synth(3400.0, 2.0, 0.0, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    assert abs(fit["period_steps"] - 3400.0) < 20.0, fit["period_steps"]
    assert abs(fit["amplitude"] - 2.0) < 0.05, fit["amplitude"]
    assert fit["explained_fraction"] > 0.999, fit["explained_fraction"]


def test_the_removal_takes_a_phased_mode_out_completely():
    """the removal leaves only the leg's own level when the mode has a phase.

    this is the case a cosine-only removal gets wrong: with a phase away
    from zero the residual retains the sine component, and the gate would
    read that residue as physics.
    """
    phase = 1.1
    series = {
        "a": synth(3400.0, 2.0, phase, 10.0, 0.001),
        "b": synth(3400.0, 2.0, phase, 12.0, 0.001),
        "c": synth(3400.0, 2.0, phase, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    residual = acoustic_mode.remove_mode(fit, "a", series["a"])
    residual = acoustic_mode.detrend_line(residual)

    values = [value for _, value in residual]
    rms = math.sqrt(sum(v * v for v in values) / len(values))
    # the constant and drift survive the removal, so the varying part must
    # be what vanishes
    assert rms < 0.05, "residual rms %s after removing a phased mode" % rms


def test_the_removal_leaves_a_different_phase_quiet():
    """a mode fitted at one phase still removes cleanly from another leg.

    the legs march together, so their phases agree; this guards the
    weaker property that a phase offset between legs does not leave a
    visible residue.
    """
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, 0.4, 12.0, 0.001),
        "c": synth(3400.0, 2.0, 0.8, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    residual = acoustic_mode.remove_mode(fit, "b", series["b"])
    residual = acoustic_mode.detrend_line(residual)
    values = [value for _, value in residual]
    rms = math.sqrt(sum(v * v for v in values) / len(values))
    assert rms < 0.2, "residual rms %s" % rms


def test_the_legs_genuinely_share_the_mode():
    """cross prediction must succeed: each leg predicts the others."""
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, 0.2, 12.0, 0.001),
        "c": synth(3400.0, 2.0, 0.4, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    sharing = acoustic_mode.mode_sharing(fit, series)
    assert sharing["shares"], sharing
    assert sharing["weakest_prediction"] > 0.9, sharing
    assert sharing["phase_coherence"] > 0.9, sharing


def test_two_unrelated_periods_are_not_called_shared():
    """legs with DIFFERENT periods must not pass the sharing test.

    this is the discriminator that separates a common chamber mode from
    legs that merely both oscillate. every leg here has its own
    frequency, so no pair shares one.
    """
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(5100.0, 2.0, 0.0, 12.0, 0.001),
        "c": synth(2300.0, 2.0, 0.0, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    sharing = acoustic_mode.mode_sharing(fit, series)
    assert not sharing["shares"], sharing


def test_a_common_frequency_without_common_phase_is_rejected():
    """the same frequency but a different phase is not a shared mode.

    the legs march from the same initial condition, so a genuine chamber
    mode arrives in step. legs that oscillate at the right frequency out
    of step are two separate events that happen to coincide in
    frequency, and subtracting one from the other would corrupt both.
    """
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, math.pi, 12.0, 0.001),
        "c": synth(3400.0, 2.0, 0.0, 14.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    sharing = acoustic_mode.mode_sharing(fit, series)
    assert sharing["phase_coherence"] < 0.9, sharing


def test_the_adjusted_metric_penalises_extra_parameters():
    """the dof-adjusted explained fraction must not flatter the model.

    a model with free per-leg phases fits better than one without, so the
    naive explained fraction rises for reasons that have nothing to do
    with the mode being real. the adjusted figure is the honest one.
    """
    series = {
        "a": synth(3400.0, 2.0, 0.0, 10.0, 0.001),
        "b": synth(3400.0, 2.0, 0.0, 12.0, 0.001),
    }
    fit = acoustic_mode.fit_shared_mode(series)
    assert fit["adjusted_explained_fraction"] > 0.9, fit["adjusted_explained_fraction"]
    assert fit["adjusted_explained_fraction"] <= fit["explained_fraction"] + 1e-9


def statistics_mean(values):
    return sum(values) / len(values)


def main():
    failures = 0
    for name, function in sorted(globals().items()):
        if not name.startswith("test_") or not callable(function):
            continue
        try:
            function()
        except AssertionError as exc:
            failures += 1
            print("FAIL %s: %s" % (name, exc))
        else:
            print("pass %s" % name)
    if failures:
        print("%d test(s) failed" % failures)
        return 1
    print("all mode tests pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
