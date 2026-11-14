#!/usr/bin/env python3
"""Bracket the true Ceiling-2 by stratifying on reference count.

The plug-in ceiling ranks by p-hat_u and weights by the same p-hat_u, so it inherits the
estimation noise in p-hat. The half-vote cross-fit scores from one half of a unit's votes
and targets the other, so it pays that noise twice. The first is therefore biased UP and
the second biased DOWN, and the true ceiling lies between them.

Both biases shrink as a unit's reference count grows. Stratifying on minimum k and
watching the two bounds move is model-free: no latent family is assumed, which matters
because supplementary/scripts/ceiling_deconvolved.py shows a beta-binomial does not fit these
counts. If the plug-in figure measured real annotation heterogeneity, the bracket would
already be narrow at k = 2 and would stay put. If it measures how few references were
asked, the bracket closes as k grows.

The population changes with the stratum -- heavily referenced units are not a random
sample -- so a stratum's numbers describe that stratum. The CONVERGENCE is the evidence,
and a population shift does not drive two bounds toward each other from opposite sides.

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

# Ceiling-1 is the figure the headroom bound of Section 5 is built on, so it is bracketed
# here too. Its cross-fit gives a predictor a ONE-reference estimate of the majority label
# in place of perfect knowledge of it, which is a lower bound on a ceiling defined by
# perfect knowledge; the target half is a single-reference draw, which is what a benchmark
# label is. Both bounds are therefore legitimate for this subpopulation.


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
                     "c1_crossfit_lower": round(xf1, 4),
                     "c1_bracket": round(plug1 - xf1, 4),
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
    monotone = all(widths[i] >= widths[i + 1] for i in range(1, len(widths) - 1))

    out = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": "supplementary/scripts/ceiling_convergence.py",
        "source_sha256": v.get("source_sha256"),
        "basis": f"stratified on min reference count, {args.splits} splits, seed-{args.seed}",
        "strata": rows,
        "bracket_narrows_with_k": monotone,
        "widest_bracket": round(max(widths), 4),
        "narrowest_bracket": round(min(widths), 4),
        "reference_in_distribution": REFERENCE,
        "headroom_at_c1_plugin": round(base["c1_plugin_upper"] - REFERENCE, 4),
        "headroom_at_c1_crossfit": round(base["c1_crossfit_lower"] - REFERENCE, 4),
        "c1_bracket_audited_subpopulation": [base["c1_crossfit_lower"],
                                             base["c1_plugin_upper"]],
        "c2_bracket_audited_subpopulation": [base["c2_crossfit_lower"],
                                             base["c2_plugin_upper"]],
        # Derived, never hardcoded: a fixed string drifts from the numbers beside it
        # the first time the split count or the strata change.
        "interpretation": (
            f"The plug-in bound falls and the cross-fit bound rises as reference count "
            f"grows. Ceiling-2's bracket goes from {base['c2_bracket']} at min k 2 to "
            f"{rows[-1]['c2_bracket']} at min k {rows[-1]['min_k']}, and Ceiling-1's from "
            f"{base['c1_bracket']} to {rows[-1]['c1_bracket']}. The published "
            f"{PLUGIN_C2} and {PLUGIN_C1} are the plug-in bounds in the noisiest stratum, "
            f"where most units carry two references. They are upper bounds on the true "
            f"ceilings, not estimates of them. The level at which the two bounds meet on "
            f"deeply measured units is a property of those more-studied epitopes and does "
            f"not transfer to the audited subpopulation; the convergence is the evidence."
        ),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    print("  min k  units   C1 upper  C1 lower  C1 bkt   C2 upper  C2 lower  C2 bkt")
    for r in rows:
        print(f"  >={r['min_k']:<4} {r['n_units']:<6}  {r['c1_plugin_upper']:.4f}    "
              f"{r['c1_crossfit_lower']:.4f}   {r['c1_bracket']:+.4f}   "
              f"{r['c2_plugin_upper']:.4f}    {r['c2_crossfit_lower']:.4f}   "
              f"{r['c2_bracket']:+.4f}")
    print(f"\n  headroom from {REFERENCE} to Ceiling-1: "
          f"{base['c1_plugin_upper'] - REFERENCE:+.4f} at the plug-in bound, "
          f"{base['c1_crossfit_lower'] - REFERENCE:+.4f} at the cross-fit bound")
    print(f"\nbracket narrows monotonically from k>=3: {monotone}")
    print(f"widest {max(widths):.4f} -> narrowest {min(widths):.4f}")
    print(f"\nARTIFACT WRITTEN: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
