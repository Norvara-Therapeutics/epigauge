#!/usr/bin/env python3
"""Bracket the true Ceiling-2 by stratifying on reference count.

The plug-in ceiling ranks by p-hat_u and weights by the same p-hat_u, so it inherits the
estimation noise in p-hat. For Ceiling-2 the half-vote cross-fit scores from one half of a
unit's votes and targets the other, so it pays that noise twice. The first is therefore
biased UP and the second biased DOWN, and the true Ceiling-2 lies between them: the half-B
label masses are unbiased for p and independent of the half-A scores, so each unordered
pair's cross-fit contribution is at most its true contribution, term by term, with equality
when p_u = p_v. That inequality is a pairwise algebraic fact and it is specific to
Ceiling-2.

Both biases shrink as a unit's reference count grows. Stratifying on minimum k and
watching the two bounds move is model-free: no latent family is assumed, which matters
because supplementary/scripts/ceiling_deconvolved.py shows a beta-binomial does not fit these
counts. If the plug-in figure measured real annotation heterogeneity, the bracket would
already be narrow at k = 2 and would stay put. If it measures how few references were
asked, the bracket closes as k grows.

The population changes with the stratum -- heavily referenced units are not a random
sample -- so a stratum's numbers describe that stratum. For Ceiling-2 the narrowing of the
bracket is the evidence, and a population shift does not drive two bounds toward each
other from opposite sides. Ceiling-1 is reported here with a plug-in upper bound only; the
half-vote statistic computed alongside it is not a bound and no Ceiling-1 bracket is
claimed. See the comment below.

    python3 supplementary/scripts/ceiling_convergence.py [--splits 400]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VOTES = ROOT / "supplementary" / "artifacts" / "annotation" / "iedb_unit_votes.json"
OUT = ROOT / "supplementary" / "artifacts" / "annotation" / "ceiling_convergence.json"

PLUGIN_C1 = 0.8237   # the published figures, which are the min-k = 2 stratum's plug-ins
PLUGIN_C2 = 0.9119
REFERENCE = 0.6125   # the in-distribution cross-validated figure the headroom measures from
STRATA = (2, 3, 4, 5, 6, 8, 10)

# Ceiling-1 CARRIES NO LOWER BOUND, and the half-vote statistic computed for it below is
# NOT one. Two reasons, either sufficient. (1) No proof: for two units with p_u > p_v
# sharing a majority label the cross-fit contribution strictly exceeds the true one,
# because maj(p-hat) is stochastically increasing in p, so a noisier but finer-grained
# score can order better than a coarse majority label. (2) The construction degenerates:
# maj() on a single vote in {0, 0.5, 1} is the identity, so on a two-vote split the
# Ceiling-1 and the Ceiling-2 cross-fit scorers are the same function, and 62.9 per cent of
# audited units carry two references. That is why the Ceiling-1 half-vote figure sits so
# close to the Ceiling-2 one: it is largely the same number relabelled. The field is
# therefore named `c1_crossfit_half_vote`, it is a descriptive statistic of one particular
# split construction, and Ceiling-1 is reported as a plug-in UPPER bound alone, 0.8237, on
# the score of one particular predictor. All headroom is reported on Ceiling-2.


def auroc(points: list[tuple[float, float]]) -> float:
    P = sum(w for _, w in points)
    N = sum(1.0 - w for _, w in points)
    if P <= 0 or N <= 0:
        return float("nan")
    g: dict[float, list[float]] = {}
    for s, w in points:
        acc = g.setdefault(s, [0.0, 0.0])
        acc[0] += w
        acc[1] += 1.0 - w
    tot = nb = 0.0
    for s in sorted(g):
        gp, gn = g[s]
        tot += gp * (nb + 0.5 * gn)
        nb += gn
    return tot / (P * N)


def population_shares(units: list[list[float]]) -> dict:
    """Descriptive shares of the audited subpopulation by reference depth.

    These are counts, not a model. An earlier version of this script also regressed
    Var(p-hat | k) on 1/k and reported the intercept as true heterogeneity and the slope as
    per-draw variance. That decomposition is withdrawn: the fitted slope exceeded its own
    model's ceiling E[p](1-E[p]) at the measured mean, and leave-one-out over the strata
    moved the implied noise share outside [0, 1]. Nothing here depends on it, and the
    bracket on Ceiling-2 is the evidence.
    """
    at2 = [sum(u) / len(u) for u in units if len(u) == 2]
    return {
        "share_of_units_with_two_refs": round(100 * len(at2) / len(units), 1),
        "share_of_two_ref_units_at_an_endpoint": round(
            100 * sum(1 for x in at2 if x in (0.0, 1.0)) / len(at2), 1),
    }


def majority(x: float) -> float:
    return 1.0 if x > 0.5 else (0.0 if x < 0.5 else 0.5)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--splits", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20261114)
    args = ap.parse_args(argv)

    v = json.loads(VOTES.read_text())
    k, pos, neg, tie = v["k"], v["n_pos_refs"], v["n_neg_refs"], v["n_tie_refs"]
    units = [[1.0] * int(pos[i]) + [0.0] * int(neg[i]) + [0.5] * int(tie[i])
             for i in range(v["n_units"]) if k[i] >= 2]

    rows = []
    for mink in STRATA:
        sel = [u for u in units if len(u) >= mink]
        if len(sel) < 60:
            continue
        ph = [sum(u) / len(u) for u in sel]
        plug2 = auroc([(x, x) for x in ph])
        plug1 = auroc([(majority(x), x) for x in ph])
        rng = random.Random(args.seed)
        v1, v2 = [], []
        for _ in range(args.splits):
            a1, a2 = [], []
            for u in sel:
                vv = list(u)
                rng.shuffle(vv)
                h = len(vv) // 2
                sa = sum(vv[:h]) / h
                tb = sum(vv[h:]) / (len(vv) - h)
                a1.append((majority(sa), tb))
                a2.append((sa, tb))
            v1.append(auroc(a1))
            v2.append(auroc(a2))
        xf1, xf2 = st.fmean(v1), st.fmean(v2)
        rows.append({"min_k": mink, "n_units": len(sel),
                     "c1_plugin_upper": round(plug1, 4),
                     "c1_crossfit_half_vote": round(xf1, 4),
                     "c1_crossfit_half_vote_is_a_bound": False,
                     "c1_plugin_minus_crossfit": round(plug1 - xf1, 4),
                     "c2_plugin_upper": round(plug2, 4),
                     "c2_crossfit_lower": round(xf2, 4),
                     "c2_bracket": round(plug2 - xf2, 4)})

    # The plug-in bound must be the published figure in the min-k = 2 stratum, or this is
    # not the artifact behind the paper.
    base = next(r for r in rows if r["min_k"] == 2)
    for got, want, nm in ((base["c1_plugin_upper"], PLUGIN_C1, "Ceiling-1"),
                          (base["c2_plugin_upper"], PLUGIN_C2, "Ceiling-2")):
        if abs(got - want) > 0.0005:
            print(f"REFUSING TO REPORT: the min-k=2 {nm} plug-in is {got}, not the "
                  f"published {want}. This is not the artifact behind the paper.",
                  file=sys.stderr)
            return 1

    widths = [r["c2_bracket"] for r in rows]

    # Every monotonicity statement below is MEASURED off the table and none is asserted.
    # A hardcoded claim about the shape of these curves is exactly what went wrong before:
    # the plug-in bound does not in fact fall at every step, and the width does not in fact
    # narrow at every step, so a fixed sentence contradicted the rows printed beside it.
    def steps_where(vals, up: bool):
        return [(rows[i]["min_k"], rows[i + 1]["min_k"])
                for i in range(len(vals) - 1)
                if (vals[i + 1] > vals[i]) == up and vals[i + 1] != vals[i]]

    def fmt(pairs):
        return ", ".join(f"min k {a} to {b}" for a, b in pairs) if pairs else "no step"

    c2_plugin = [r["c2_plugin_upper"] for r in rows]
    c2_lower = [r["c2_crossfit_lower"] for r in rows]
    width_rises = steps_where(widths, up=True)
    plugin_rises = steps_where(c2_plugin, up=True)
    lower_falls = steps_where(c2_lower, up=False)
    shares = population_shares(units)

    out = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": "supplementary/scripts/ceiling_convergence.py",
        "source_sha256": v.get("source_sha256"),
        "basis": f"stratified on min reference count, {args.splits} splits, seed-{args.seed}",
        "strata": rows,
        "population_shares": shares,
        "c2_width_narrows_overall": widths[-1] < widths[0],
        "c2_width_rises_at": [list(p) for p in width_rises],
        "c2_plugin_rises_at": [list(p) for p in plugin_rises],
        "c2_crossfit_falls_at": [list(p) for p in lower_falls],
        "widest_c2_bracket": round(max(widths), 4),
        "narrowest_c2_bracket": round(min(widths), 4),
        "reference_in_distribution": REFERENCE,
        # Headroom is reported on Ceiling-2 only. Ceiling-1 has no lower bound, so it
        # supports no headroom figure at a lower bound and none is written here.
        "headroom_at_c2_plugin": round(base["c2_plugin_upper"] - REFERENCE, 4),
        "headroom_at_c2_crossfit": round(base["c2_crossfit_lower"] - REFERENCE, 4),
        "c1_plugin_upper_audited_subpopulation": base["c1_plugin_upper"],
        "c1_lower_bound": None,
        "c1_lower_bound_note": (
            "Withdrawn. The half-vote construction degenerates on Ceiling-1 -- maj() on a "
            "single vote in {0, 0.5, 1} is the identity, so on a two-vote split the two "
            f"cross-fit scorers are the same function and "
            f"{shares['share_of_units_with_two_refs']} per cent of audited units carry two "
            "references -- and no proof bounds Ceiling-1 below. `c1_crossfit_half_vote` is "
            "a descriptive statistic of that construction and is not a bound."
        ),
        "c2_bracket_audited_subpopulation": [base["c2_crossfit_lower"],
                                             base["c2_plugin_upper"]],
        # Derived, never hardcoded: a fixed string drifts from the numbers beside it
        # the first time the split count or the strata change.
        "interpretation": (
            f"Ceiling-2 is bracketed on the audited subpopulation by "
            f"{base['c2_crossfit_lower']:.4f} below and "
            f"{base['c2_plugin_upper']:.4f} above. The "
            f"plug-in bound is biased up because it ranks and weights by the same p-hat; "
            f"the half-vote cross-fit bound is biased down, and provably so, because the "
            f"half-B label masses are unbiased for p and independent of the half-A scores. "
            f"The width narrows overall, from {widths[0]:.4f} at min k "
            f"{rows[0]['min_k']} to {widths[-1]:.4f} at min k {rows[-1]['min_k']}, "
            f"and it does so non-monotonically: "
            f"the width rises at {fmt(width_rises)}, the plug-in bound rises at "
            f"{fmt(plugin_rises)}, and the cross-fit bound falls at {fmt(lower_falls)}. "
            f"The published {PLUGIN_C2} is the plug-in bound in the noisiest stratum, "
            f"where {shares['share_of_units_with_two_refs']} per cent of units carry two "
            f"references, so it is an upper bound on the true ceiling and not an estimate "
            f"of it. Ceiling-1 carries an upper bound only: {PLUGIN_C1} is a plug-in upper "
            f"bound on the score of one particular predictor, not a bound on model "
            f"performance, and `c1_crossfit_half_vote` is not a lower bound on it. "
            f"Headroom above the in-distribution cross-validated {REFERENCE} is therefore "
            f"reported on Ceiling-2 alone: at most "
            f"{base['c2_plugin_upper'] - REFERENCE:.4f}, with "
            f"{base['c2_crossfit_lower'] - REFERENCE:.4f} not excluded. The level at "
            f"which Ceiling-2's two bounds approach each other on deeply measured units is "
            f"a property of those more-studied epitopes and does not transfer to the "
            f"audited subpopulation."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    print("  min k  units   C1 upper  C1 half-vote (not a bound)   "
          "C2 upper  C2 lower  C2 width")
    for r in rows:
        print(f"  >={r['min_k']:<4} {r['n_units']:<6}  {r['c1_plugin_upper']:.4f}    "
              f"{r['c1_crossfit_half_vote']:.4f}                     "
              f"{r['c2_plugin_upper']:.4f}    {r['c2_crossfit_lower']:.4f}   "
              f"{r['c2_bracket']:+.4f}")
    print(f"\n  {shares['share_of_units_with_two_refs']}% of audited units carry two "
          f"references, and {shares['share_of_two_ref_units_at_an_endpoint']}% of those "
          f"sit exactly at an endpoint")
    print(f"  headroom from {REFERENCE} to Ceiling-2: "
          f"{base['c2_plugin_upper'] - REFERENCE:+.4f} at the plug-in upper bound, "
          f"{base['c2_crossfit_lower'] - REFERENCE:+.4f} at the cross-fit lower bound")
    print("  Ceiling-1 has no lower bound; its half-vote figure is not one")
    print(f"\nC2 width narrows overall: {widths[-1] < widths[0]}   "
          f"(rises at {fmt(width_rises)})")
    print(f"widest {max(widths):.4f} -> narrowest {min(widths):.4f}")
    print(f"\nARTIFACT WRITTEN: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
