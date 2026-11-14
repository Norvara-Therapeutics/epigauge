#!/usr/bin/env python3
"""Separate real annotation heterogeneity from reference-sampling noise in Ceiling-2.

WHAT THIS TESTS. The published Ceiling-2 of 0.9119 is a plug-in: it ranks by p-hat_u and
weights by the same p-hat_u. p-hat_u is a mean of k reference votes, and most audited
units have k = 2, so p-hat is a very noisy estimate of the true disagreement rate p_u.
Noise in p-hat pushes it toward 0 and 1 -- 61.3 per cent of k = 2 units sit exactly at an
endpoint -- and endpoints are maximally separable. A plug-in ceiling therefore measures
partly how much annotators really disagree and partly how few of them were asked.

The quantity the headroom argument needs is the ceiling under the TRUE p_u distribution:
the AUROC a predictor could reach if it learned each unit's real disagreement rate. This
script estimates that distribution by beta-binomial maximum likelihood, which treats each
unit's votes as Binomial(k_u, p_u) with p_u ~ Beta(a, b), and then evaluates Ceiling-1 and
Ceiling-2 on the fitted latent distribution instead of on the empirical p-hat.

Ties are excluded here rather than counted as half a vote. A half-vote is not a binomial
outcome, and the published p_u definition -- (pos + tie/2)/k -- has no beta-binomial
likelihood. The restriction is stated in the artifact as a scope condition: this estimates
the latent distribution of the DECISIVE vote rate.

TWO INDEPENDENT CHECKS must agree before the result is reportable:
  1. the fitted model must reproduce the observed variance of p-hat at each k, and
  2. the deconvolved ceiling must land near the half-vote cross-fit of
     supplementary/scripts/ceiling_crossfit.py, which reaches the same question by splitting
     votes rather than by modelling them.
Disagreement between them means the beta-binomial is the wrong latent family, and the
script says so rather than reporting a number.

    python3 supplementary/scripts/ceiling_deconvolved.py [--grid 4000]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VOTES = ROOT / "supplementary" / "artifacts" / "annotation" / "iedb_unit_votes.json"
CROSSFIT = ROOT / "supplementary" / "artifacts" / "annotation" / "ceiling_crossfit.json"
OUT = ROOT / "supplementary" / "artifacts" / "annotation" / "ceiling_deconvolved.json"

PLUGIN_C1 = 0.8237
PLUGIN_C2 = 0.9119


def lbeta(a: float, b: float) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def nll(a: float, b: float, data: list[tuple[int, int]]) -> float:
    """Negative log-likelihood of Beta(a,b)-Binomial over (successes, trials)."""
    if a <= 0 or b <= 0:
        return float("inf")
    base = lbeta(a, b)
    tot = 0.0
    for s, n in data:
        tot -= lbeta(a + s, b + n - s) - base
    return tot


def fit(data: list[tuple[int, int]]) -> tuple[float, float]:
    """Nelder-Mead on (log a, log b); no SciPy, so the simplex is written out."""
    import itertools

    def f(v):
        return nll(math.exp(v[0]), math.exp(v[1]), data)

    # Start from the method-of-moments fit, which is already close.
    ph = [s / n for s, n in data]
    m, var = st.fmean(ph), st.pvariance(ph)
    mean_inv_k = st.fmean([1.0 / n for _, n in data])
    v_true = max(var - m * (1 - m) * mean_inv_k, 1e-4)
    conc = max(m * (1 - m) / v_true - 1.0, 0.5)
    simplex = [[math.log(m * conc), math.log((1 - m) * conc)]]
    simplex += [[simplex[0][0] + 0.5, simplex[0][1]], [simplex[0][0], simplex[0][1] + 0.5]]
    vals = [f(p) for p in simplex]
    for _ in range(600):
        order = sorted(range(3), key=lambda i: vals[i])
        simplex = [simplex[i] for i in order]; vals = [vals[i] for i in order]
        if abs(vals[2] - vals[0]) < 1e-10:
            break
        cen = [(simplex[0][j] + simplex[1][j]) / 2 for j in (0, 1)]
        for coef in (2.0, 1.5, 0.5, -0.5):   # reflect, expand, contract, shrink
            trial = [cen[j] + coef * (cen[j] - simplex[2][j]) for j in (0, 1)]
            tv = f(trial)
            if tv < vals[2]:
                simplex[2], vals[2] = trial, tv
                break
    return math.exp(simplex[0][0]), math.exp(simplex[0][1])


def ceilings_from_density(xs: list[float], ws: list[float]) -> tuple[float, float]:
    """Ceiling-1 (majority-label ranker) and Ceiling-2 (Bayes-optimal) on a p density."""
    def auroc(score: list[float]) -> float:
        pts = sorted(zip(score, xs, ws))
        P = sum(w * x for _, x, w in pts)
        N = sum(w * (1 - x) for _, x, w in pts)
        g: dict[float, list[float]] = {}
        for s, x, w in pts:
            acc = g.setdefault(s, [0.0, 0.0])
            acc[0] += w * x
            acc[1] += w * (1 - x)
        tot = nb = 0.0
        for s in sorted(g):
            gp, gn = g[s]
            tot += gp * (nb + 0.5 * gn)
            nb += gn
        return tot / (P * N)

    c2 = auroc(list(xs))
    c1 = auroc([1.0 if x > 0.5 else 0.0 if x < 0.5 else 0.5 for x in xs])
    return c1, c2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", type=int, default=4000)
    args = ap.parse_args(argv)

    v = json.loads(VOTES.read_text())
    k, pos, neg = v["k"], v["n_pos_refs"], v["n_neg_refs"]
    data = [(int(pos[i]), int(pos[i] + neg[i]))
            for i in range(v["n_units"])
            if k[i] >= 2 and (pos[i] + neg[i]) >= 2]
    if len(data) < 500:
        print("could not run: too few units with two decisive votes", file=sys.stderr)
        return 2

    a, b = fit(data)
    mean = a / (a + b)
    var = a * b / ((a + b) ** 2 * (a + b + 1))

    # CHECK 1: the fitted model must reproduce the observed spread of p-hat at each k.
    rows = []
    worst = 0.0
    for kk in range(2, 7):
        grp = [(s, n) for s, n in data if n == kk]
        if len(grp) < 50:
            continue
        obs = st.pvariance([s / n for s, n in grp])
        pred = var + mean * (1 - mean) / kk - var / kk      # beta-binomial p-hat variance
        rows.append({"k": kk, "n": len(grp), "observed_var": round(obs, 5),
                     "predicted_var": round(pred, 5)})
        worst = max(worst, abs(obs - pred))

    xs = [(i + 0.5) / args.grid for i in range(args.grid)]
    logw = [(a - 1) * math.log(x) + (b - 1) * math.log(1 - x) for x in xs]
    mx = max(logw)
    ws = [math.exp(w - mx) for w in logw]
    c1, c2 = ceilings_from_density(xs, ws)

    # CHECK 2: agreement with the split-based cross-fit, which models nothing.
    xf = None
    if CROSSFIT.is_file():
        xf = json.loads(CROSSFIT.read_text()).get("ceiling_2_crossfit_mean")

    ok_var = worst < 0.01
    ok_xf = xf is None or abs(c2 - xf) < 0.08
    if not (ok_var and ok_xf):
        print("REFUSING TO REPORT a deconvolved ceiling.", file=sys.stderr)
        if not ok_var:
            print(f"  the fitted beta-binomial misses the observed variance of p-hat by "
                  f"{worst:.4f} at some k; the latent family is wrong.", file=sys.stderr)
        if not ok_xf:
            print(f"  deconvolved Ceiling-2 {c2:.4f} disagrees with the model-free "
                  f"cross-fit {xf:.4f} by more than 0.08.", file=sys.stderr)
        for r in rows:
            print(f"    k={r['k']:<2} n={r['n']:<5} observed {r['observed_var']:.5f} "
                  f"predicted {r['predicted_var']:.5f}", file=sys.stderr)
        return 1

    out = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": "supplementary/scripts/ceiling_deconvolved.py",
        "source_sha256": v.get("source_sha256"),
        "basis": f"beta-binomial MLE, decisive votes only, {args.grid}-point density",
        "scope": "latent distribution of the DECISIVE vote rate; tie votes excluded",
        "n_units_fitted": len(data),
        "beta_a": round(a, 4),
        "beta_b": round(b, 4),
        "latent_mean_p": round(mean, 4),
        "latent_var_p": round(var, 5),
        "latent_sd_p": round(math.sqrt(var), 4),
        "variance_check": rows,
        "variance_check_worst_abs_gap": round(worst, 5),
        "ceiling_1_plugin": PLUGIN_C1,
        "ceiling_2_plugin": PLUGIN_C2,
        "ceiling_1_deconvolved": round(c1, 4),
        "ceiling_2_deconvolved": round(c2, 4),
        "ceiling_2_crossfit_independent": xf,
        "plugin_inflation_c2": round(PLUGIN_C2 - c2, 4),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    print(f"units fitted (>=2 decisive votes)   {len(data)}")
    print(f"latent Beta(a, b)                   a={a:.4f}  b={b:.4f}")
    print(f"latent p_u: mean {mean:.4f}, sd {math.sqrt(var):.4f}\n")
    print("variance check (fitted vs observed spread of p-hat):")
    for r in rows:
        print(f"    k={r['k']:<2} n={r['n']:<5} observed {r['observed_var']:.5f} "
              f"predicted {r['predicted_var']:.5f}")
    print(f"  worst gap {worst:.5f}\n")
    print(f"Ceiling-1  plug-in {PLUGIN_C1}   deconvolved {c1:.4f}")
    print(f"Ceiling-2  plug-in {PLUGIN_C2}   deconvolved {c2:.4f}")
    if xf is not None:
        print(f"           independent half-vote cross-fit {xf:.4f}")
    print(f"\nplug-in inflation of Ceiling-2       {PLUGIN_C2 - c2:+.4f}")
    print(f"\nARTIFACT WRITTEN: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
