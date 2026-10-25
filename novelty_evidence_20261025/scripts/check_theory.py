"""Small numerical checks of restricted shrinkage claims, not full-method guarantees.

Only NumPy and the Python standard library are required. No MATLAB invocation.
The script is intentionally independent of the seismic processing pipeline.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


SEED = 20261025
N = 50_000
BATCH = 2_000


def json_default(value):
    """Convert NumPy scalar results without altering the numerical checks."""
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported JSON value: {type(value).__name__}")


class Moments:
    def __init__(self):
        self.n = 0
        self.total = 0.0
        self.square_total = 0.0

    def add(self, values):
        values = np.asarray(values, dtype=float)
        self.n += values.size
        self.total += float(values.sum())
        self.square_total += float(np.square(values).sum())

    def summary(self):
        mean = self.total / self.n
        variance = max(0.0, (self.square_total - self.n * mean * mean) / (self.n - 1))
        return {"n": self.n, "mean": mean, "standard_error": (variance / self.n) ** 0.5}


def gain(energy, variance, lam):
    """Nonnegative inputs; zero denominator convention is gain=1."""
    energy = np.asarray(energy, dtype=float)
    denominator = energy + np.asarray(lam) * np.asarray(variance)
    return np.divide(energy, denominator, out=np.ones_like(energy), where=denominator != 0)


def fixed_gain_check(rng, signal_energy, lam, d=4, sigma=1.0, fixed_g=None):
    variance = d * sigma * sigma
    g = float(gain(signal_energy, variance, lam)) if fixed_g is None else float(fixed_g)
    mu = np.zeros(d)
    mu[0] = signal_energy ** 0.5
    losses, paired = Moments(), Moments()
    for start in range(0, N, BATCH):
        eta = rng.normal(0.0, sigma, size=(min(BATCH, N - start), d))
        z = mu + eta
        loss = np.square(g * z - mu).sum(axis=1)
        losses.add(loss)
        paired.add(loss - np.square(eta).sum(axis=1))
    empirical, difference = losses.summary(), paired.summary()
    expected = (1.0 - g) ** 2 * signal_energy + g * g * variance
    expected_difference = expected - variance
    tolerance = 6.0 * empirical["standard_error"] + 1e-10
    return {
        "signal_energy_A": signal_energy, "noise_energy_V": variance,
        "lambda": lam if fixed_g is None else None,
        "gain_kind": "true_signal_energy_oracle" if fixed_g is None else "fixed_external_gain",
        "gain": g, "dimension": d, "expected_risk": expected,
        "expected_difference_from_unshrunk": expected_difference,
        "empirical_risk": empirical, "paired_risk_difference": difference,
        "agreement_pass": abs(empirical["mean"] - expected) <= tolerance,
        "analytical_nonincrease": expected_difference <= 1e-12,
    }


def stein_check(rng, d, mu_norm, lam, sigma=1.0):
    """Isotropic Gaussian, externally fixed support, constant lambda only."""
    mu = np.zeros(d)
    mu[0] = mu_norm
    variance = d * sigma * sigma
    a = lam * variance
    paired, stein, discrepancy = Moments(), Moments(), Moments()
    for start in range(0, N, BATCH):
        eta = rng.normal(0.0, sigma, size=(min(BATCH, N - start), d))
        z = mu + eta
        t = np.square(z).sum(axis=1)
        g = gain(t, variance, lam)
        delta = np.square(g[:, None] * z - mu).sum(axis=1) - np.square(eta).sum(axis=1)
        estimate = a * ((a + (4 - 2 * d) * sigma * sigma) * t - 2 * d * sigma * sigma * a) / np.square(t + a)
        paired.add(delta)
        stein.add(estimate)
        discrepancy.add(delta - estimate)
    p, s, e = paired.summary(), stein.summary(), discrepancy.summary()
    sufficient = d >= 3 and 0.0 < a <= 2 * (d - 2) * sigma * sigma
    counterexample = d == 1 and mu_norm == 5.0 and lam == 1.0
    return {
        "dimension": d, "mu_norm": mu_norm, "sigma": sigma, "lambda": lam,
        "a": a, "energy_kind": "observed_squared_norm",
        "paired_risk_difference": p, "stein_risk_difference": s,
        "paired_minus_stein": e,
        "agreement_pass": abs(e["mean"]) <= 6.0 * e["standard_error"] + 1e-10,
        "sufficient_uniform_improvement_condition": sufficient,
        "observed_improvement_pass": p["mean"] < 0 if sufficient else None,
        "scalar_high_signal_counterexample": counterexample,
        "counterexample_pass": p["mean"] > 3 * p["standard_error"] if counterexample else None,
    }


def main():
    rng = np.random.default_rng(SEED)
    energies=np.exp(rng.uniform(-6,3,100));variances=np.exp(rng.uniform(-6,3,100));coherence=rng.uniform(0,1,100)
    lam=.5+15.5*(1-coherence)**2
    optimal=energies/(energies+lam*variances)
    proxy=(1-optimal)**2*energies+optimal**2*lam*variances
    closed=energies*lam*variances/(energies+lam*variances)
    stationarity=2*((energies+lam*variances)*optimal-energies)
    proxy_check=dict(max_stationarity=float(np.max(np.abs(stationarity))),max_objective_difference=float(np.max(np.abs(proxy-closed))))
    proxy_check['pass']=proxy_check['max_stationarity']<1e-12 and proxy_check['max_objective_difference']<1e-12
    a = 1.0
    z0 = np.sqrt(3.0 * a)
    eps = 1e-6

    def transform(z):
        return z ** 3 / (z * z + a)

    finite_difference = (transform(z0 + eps) - transform(z0 - eps)) / (2 * eps)
    derivative = {
        "a": a, "z": float(z0), "finite_difference": float(finite_difference),
        "expected": 9.0 / 8.0,
        "pass": abs(finite_difference - 9.0 / 8.0) < 1e-7,
        "interpretation": "Gain <= 1 does not imply a nonexpansive input-output map.",
    }

    synthesis = np.array([[1.0, 1.0]])
    coefficients = np.array([1.0, -1.0])
    gains = np.array([0.8, 0.2])
    original = synthesis @ coefficients
    shrunk = synthesis @ (gains * coefficients)
    interference = {
        "synthesis_matrix": synthesis.tolist(), "coefficients": coefficients.tolist(),
        "gains": gains.tolist(), "before": original.tolist(), "after": shrunk.tolist(),
        "coefficient_energy_before": float(np.square(coefficients).sum()),
        "coefficient_energy_after": float(np.square(gains * coefficients).sum()),
        "waveform_energy_before": float(np.square(original).sum()),
        "waveform_energy_after": float(np.square(shrunk).sum()),
        "pass": bool(np.square(shrunk).sum() > np.square(original).sum()),
    }
    fixed = [fixed_gain_check(rng, A, lam) for A in (0.25, 4.0, 25.0) for lam in (0.5, 1.0, 2.0, 5.0)]
    fixed.append(fixed_gain_check(rng, 25.0, None, fixed_g=0.05))
    stein_cases = [stein_check(rng, 1, 5.0, 1.0), stein_check(rng, 3, 3.0, 0.5), stein_check(rng, 4, 3.0, 1.0)]
    zero_case = {"E": 0.0, "V": 0.0, "lambda": 1.0, "gain": float(gain(0.0, 0.0, 1.0))}
    zero_case["pass"] = zero_case["gain"] == 1.0
    all_pass = proxy_check['pass'] and derivative["pass"] and interference["pass"] and zero_case["pass"]
    all_pass = all_pass and all(x["agreement_pass"] for x in fixed)
    all_pass = all_pass and all(x["agreement_pass"] for x in stein_cases)
    all_pass = all_pass and all(x["observed_improvement_pass"] is not False for x in stein_cases)
    all_pass = all_pass and all(x["counterexample_pass"] is not False for x in stein_cases)
    result = {
        "seed": SEED, "draws_per_monte_carlo_case": N, "batch_size": BATCH,
        "scope": "Restricted numerical checks only; no literature-first claim or full-pipeline guarantee.",
        "limitations": [
            "Gaussian cases use externally fixed support and constant lambda.",
            "Oracle gain uses true signal energy; Stein cases use observed energy.",
            "An arbitrary correlated source-domain covariance is not isotropic Gaussian noise.",
            "Same-data coherence, support selection, Q selection, and model bias are not tested here.",
        ],
        "zero_denominator_convention": zero_case,
        "observable_quadratic_proxy": proxy_check,
        "derivative_9_over_8": derivative, "synthesis_interference": interference,
        "fixed_gain_risk": fixed, "stein_observed_energy": stein_cases,
        "all_checks_pass": bool(all_pass),
    }
    output = Path(__file__).resolve().parents[1] / "analysis" / "theory_checks.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=json_default), encoding="utf-8")
    print(json.dumps({"output": str(output), "all_checks_pass": bool(all_pass)}, ensure_ascii=False))
    if not all_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
