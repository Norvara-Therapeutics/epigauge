#!/usr/bin/env python3
"""Cross-fit Ceiling-2, which the published figure estimates in-sample.

THE OBJECTION THIS ANSWERS, disclosed in the preprint rather than settled: Ceiling-2 is
the AUROC of a ranker that scores by p_u, evaluated against a label distribution that is
also p_u. The same estimate supplies the score and defines the target, so the figure is
optimistically biased. With two decisive votes p_u lands in {0, 0.5, 1} and the ranker is
correlated with its own target by construction.

The fix is a split. For each unit, partition its decisive reference votes into two halves,
estimate p_u from half A, and evaluate against the label mass of half B. No vote supplies
both a score and its own target. Averaged over random partitions this is an out-of-sample
Ceiling-2, and the difference from the in-sample figure is the optimism.

Ceiling-1 is NOT cross-fit here. It scores by the majority label, which on a two-vote
split is the single vote in half A -- the same quantity as the score, so the construction
degenerates rather than correcting anything. Section 6q already establishes that Ceiling-1
is the draw-sensitive term and that the headroom argument rests on Ceiling-2.

Standard library only, matching the shipped supplementary scripts, so this runs on a bare
interpreter with no environment to reconstruct.

    python3 supplementary/scripts/ceiling_crossfit.py [--splits 400] [--seed 20261114]

The default split count is 400, which is the count recorded in the shipped artifact's
basis string, so the command above with no arguments regenerates the shipped record.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VOTES = ROOT / "supplementary" / "artifacts" / "annotation" / "iedb_unit_votes.json"
OUT = ROOT / "supplementary" / "artifacts" / "annotation" / "ceiling_crossfit.json"

EXPECT_C2 = 0.9119   # the audited in-sample figure this must reproduce before reporting
TOL = 0.0001         # p_u is reconstructed exactly, so demand an exact match

# p_u IS NOT pos/k AND IS NOT pos/(pos+neg). A reference whose rows split evenly is a tie
# vote counting half: p_u = (pos + tie/2) / k, which the table also stores as `p_ref`. The
# two wrong definitions give 0.9080 and 0.9523 against a published 0.9119, and the gate
# below refused both before they could be reported. Ties are votes, not missing data, so
# they are carried into the split rather than dropped.


def auroc_weighted(points: list[tuple[float, float]]) -> float:
    """Exact weighted Mann-Whitney: each unit carries positive mass w and negative 1-w."""
    P = sum(w for _, w in points)
    N = sum(1.0 - w for _, w in points)
    if P <= 0 or N <= 0:
        return float("nan")
    by_score: dict[float, list[float]] = {}
    for s, w in points:
        acc = by_score.setdefault(s, [0.0, 0.0])
        acc[0] += w
        acc[1] += 1.0 - w
    total = 0.0
    neg_below = 0.0
    for s in sorted(by_score):
        gp, gn = by_score[s]
        total += gp * (neg_below + 0.5 * gn)   # ties take half credit
        neg_below += gn
    return total / (P * N)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--splits", type=int, default=400,
                    help="number of random half-vote partitions to average over "
                         "(default 400, the count behind the shipped artifact)")
    ap.add_argument("--seed", type=int, default=20261114)
    args = ap.parse_args(argv)

    if not VOTES.is_file():
        print(f"could not run: {VOTES} is missing", file=sys.stderr)
        return 2
    v = json.loads(VOTES.read_text())
    k, pos, neg = v["k"], v["n_pos_refs"], v["n_neg_refs"]

    tie = v["n_tie_refs"]
    # Every unit in the audited subpopulation has k >= 2 references, so every one of them
    # has two votes to split. Each unit becomes a multiset of vote values in {0, 0.5, 1}.
    units = [[1.0] * int(pos[i]) + [0.0] * int(neg[i]) + [0.5] * int(tie[i])
             for i in range(v["n_units"]) if k[i] >= 2]
    n_audited = len(units)

    # Gate: the in-sample statistic on this subset must sit near the audited figure. The
    # subset is not the whole audited set, so an exact match is not expected; a wide drift
    # would mean this vote table is not the one behind the published number.
    p_full = [sum(u) / len(u) for u in units]
    c2_in = auroc_weighted([(x, x) for x in p_full])
    if abs(c2_in - EXPECT_C2) > TOL:
        print(f"REFUSING TO REPORT: in-sample Ceiling-2 on the decisive-vote subset is "
              f"{c2_in:.4f}, more than {TOL} from the audited {EXPECT_C2}. This vote table "
              "is not the one behind the published figure.", file=sys.stderr)
        return 1

    rng = random.Random(args.seed)
    vals: list[float] = []
    for _ in range(args.splits):
        pts = []
        for u in units:
            votes = list(u)
            rng.shuffle(votes)
            half = len(votes) // 2
            a, b = votes[:half], votes[half:]
            if not a or not b:
                continue
            pts.append((sum(a) / len(a), sum(b) / len(b)))  # score from A, target from B
        vals.append(auroc_weighted(pts))

    mean = statistics.fmean(vals)
    sd = statistics.stdev(vals)
    out = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": "supplementary/scripts/ceiling_crossfit.py",
        "source_sha256": v.get("source_sha256"),
        "basis": f"split/{args.splits}/fresh/seed-{args.seed}/half-vote-partition",
        "n_units_audited": n_audited,
        "n_units_crossfittable": len(units),
        "ceiling_2_in_sample_published": EXPECT_C2,
        "ceiling_2_in_sample_on_subset": round(c2_in, 4),
        "ceiling_2_crossfit_mean": round(mean, 4),
        "ceiling_2_crossfit_sd": round(sd, 4),
        "ceiling_2_crossfit_range": [round(min(vals), 4), round(max(vals), 4)],
        "optimism_in_sample_minus_crossfit": round(c2_in - mean, 4),
    }
    OUT.write_text(json.dumps(out, indent=2) + "\n")

    print(f"audited subpopulation           {n_audited} units")
    print(f"cross-fittable (k>=2, all of them) {len(units)} units")
    print(f"basis                           {out['basis']}\n")
    print(f"Ceiling-2 in-sample, published       {EXPECT_C2}")
    print(f"Ceiling-2 in-sample, this subset     {c2_in:.4f}")
    print(f"Ceiling-2 CROSS-FIT, mean of splits  {mean:.4f}"
          f"   (sd {sd:.4f}, range {min(vals):.4f}-{max(vals):.4f})")
    print(f"\noptimism of the in-sample figure      {c2_in - mean:+.4f}")
    print(f"\nARTIFACT WRITTEN: {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
