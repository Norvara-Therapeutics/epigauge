#!/usr/bin/env python3
"""Why allele is the stratification variable, measured rather than assumed.

A referee's objection, and it is correct: the preprint argues carefully for the WEIGHTS
given allele strata and treats the choice of allele as self-evident. The evidence it
offers -- an allele lookup scoring 0.7122 pooled and 0.5000 within allele -- is equally
available for any imbalanced covariate, and the paper itself notes that peptide length
wins the neoantigen arm, so length carries rankable signal.

This answers it in three parts.

  1. WHICH VARIABLES CAN SERVE AT ALL. The shipped schema offers five candidates. Three
     are refused by the estimator's own degeneracy guards, which is a measurement, not an
     argument: an all-singleton grouping returns exactly 0.5000 regardless of the scores,
     and a constant grouping equals pooled AUROC to machine precision. The guards are
     demonstrated firing, not asserted.

  2. THE LOOKUP DIAGNOSTIC, RUN FOR EVERY SURVIVOR. For each candidate V, build a
     predictor that sees only V -- it emits the positive rate of the stratum a row falls
     in, so it is constant within stratum by construction -- and report its pooled and
     within-V AUROC. A large pooled figure with within-V at exactly 0.5000 says V carries
     confounding that a pooled estimator will reward. The variable with the largest gap is
     the one worth removing.

  3. RANK-ORDER SENSITIVITY. Re-rank every shipped score column under each candidate's
     strata and report how far entries move. If the ranking is insensitive to the choice,
     the choice does not matter and the paper should say so; if it is sensitive, the
     justification has to be the diagnostic above and not taste.

The lookup is fitted on the arm it scores, which inflates it. That is deliberate and
stated: the diagnostic asks how much signal V carries about the label at all, not how
well a V-only model would generalise. An out-of-sample version would understate the
confounding a pooled estimator rewards.

Needs nothing beyond the shipped benchmark rows and a standard-library interpreter.

    python3 -I supplementary/scripts/stratum_choice.py
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

# The released package root: `supplementary/scripts/x.py` -> two levels up. Inserted
# explicitly rather than relying on PYTHONPATH, so `python3 -I` works from any cwd.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SCRIPT_REL = "supplementary/scripts/stratum_choice.py"
OUT = ROOT / "supplementary" / "artifacts" / "stratum_choice.json"
ARMS = ("titanian_infectious", "titanian_neoantigen")
CANDIDATES = ("allele", "length", "leakage_group", "mhc_class", "label_origin")


def lookup_scores(rows, field):
    """Positive rate of the row's own stratum: constant within stratum by construction."""
    agg: dict[str, list[int]] = {}
    for r in rows:
        agg.setdefault(str(r.get(field)), []).append(int(r["label"]))
    rate = {k: sum(v) / len(v) for k, v in agg.items()}
    return [rate[str(r.get(field))] for r in rows]


def main() -> int:
    try:
        from epigauge import metrics
    except ImportError as exc:  # pragma: no cover - environment, not a data finding
        sys.stderr.write(
            f"cannot import epigauge from {ROOT}: {exc}\n"
            "Run this from an unpacked epigauge source distribution, or install the "
            "package first.\n"
        )
        return 2

    results = {}
    for arm in ARMS:
        from epigauge import data
        rows = data.load(arm, scores=True)
        labels = [int(r["label"]) for r in rows]
        arm_rec = {"n_rows": len(rows), "candidates": {}, "rankings": {}}

        for field in CANDIDATES:
            groups = [str(r.get(field)) for r in rows]
            n_strata = len(set(groups))
            sc = lookup_scores(rows, field)
            entry = {"n_strata": n_strata, "pooled_auroc_of_lookup": None,
                     "within_auroc_of_lookup": None, "guard": None}
            entry["pooled_auroc_of_lookup"] = round(metrics.auroc(labels, sc), 4)
            try:
                st = metrics.stratified_auroc(labels, sc, groups)
                entry["within_auroc_of_lookup"] = round(st.value, 4)
                entry["strata_admissible"] = (
                    f"{st.n_strata_contributing}/{st.n_strata_total}")
                entry["rows_evaluable"] = f"{st.n_rows_evaluable}/{st.n_rows_total}"
                entry["gap_pooled_minus_within"] = round(
                    entry["pooled_auroc_of_lookup"] - st.value, 4)
                entry["usable_as_stratum"] = True
            except Exception as e:
                # A refused grouping is the finding, not an error to route around.
                entry["guard"] = f"{type(e).__name__}: {e}"
                entry["usable_as_stratum"] = False
            arm_rec["candidates"][field] = entry

        # Rank-order sensitivity over the shipped score columns.
        cols = sorted({c for r in rows for c in r.get("scores", {})})
        usable = [f for f, e in arm_rec["candidates"].items()
                  if e["usable_as_stratum"]]
        per_field: dict[str, dict[str, float]] = {}
        for field in usable:
            groups = [str(r.get(field)) for r in rows]
            vals = {}
            for c in cols:
                sc = [r.get("scores", {}).get(c) for r in rows]
                if any(v is None for v in sc):
                    continue
                try:
                    vals[c] = metrics.stratified_auroc(
                        labels, [float(v) for v in sc], groups).value
                except Exception:
                    continue
            per_field[field] = vals

        base = "allele"
        if base in per_field:
            order_base = sorted(per_field[base], key=lambda c: -per_field[base][c])
            rank_base = {c: i + 1 for i, c in enumerate(order_base)}
            for field, vals in per_field.items():
                if field == base:
                    continue
                shared = [c for c in order_base if c in vals]
                order = sorted(shared, key=lambda c: -vals[c])
                rank = {c: i + 1 for i, c in enumerate(order)}
                shifts = {c: abs(rank[c] - rank_base[c]) for c in shared}
                arm_rec["rankings"][f"{base}_vs_{field}"] = {
                    "n_entries": len(shared),
                    "max_abs_shift": max(shifts.values()) if shifts else None,
                    "mean_abs_shift": (round(sum(shifts.values()) / len(shifts), 2)
                                       if shifts else None),
                    "n_moving_3_or_more": sum(1 for v in shifts.values() if v >= 3),
                    "largest": sorted(
                        ({"entry": c, "rank_allele": rank_base[c],
                          f"rank_{field}": rank[c], "shift": shifts[c]}
                         for c in shared), key=lambda d: -d["shift"])[:4],
                }
        results[arm] = arm_rec

        print(f"=== {arm} ({len(rows)} rows) ===")
        print(f"{'candidate':16s} {'strata':>7s} {'lookup pooled':>14s} "
              f"{'lookup within':>14s} {'gap':>8s}")
        for field, e in arm_rec["candidates"].items():
            if not e["usable_as_stratum"]:
                print(f"{field:16s} {e['n_strata']:7d}   REFUSED BY GUARD: "
                      f"{e['guard'][:58]}")
                continue
            print(f"{field:16s} {e['n_strata']:7d} {e['pooled_auroc_of_lookup']:14.4f} "
                  f"{e['within_auroc_of_lookup']:14.4f} "
                  f"{e['gap_pooled_minus_within']:+8.4f}")
        for k, v in arm_rec["rankings"].items():
            print(f"  ranking {k}: {v['n_entries']} entries, max shift "
                  f"{v['max_abs_shift']}, mean {v['mean_abs_shift']}, "
                  f">=3 places: {v['n_moving_3_or_more']}")
        print()

    rec = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": SCRIPT_REL,
        "note": ("the lookup is fitted on the arm it scores; the diagnostic asks how much "
                 "label signal the variable carries, not how a V-only model generalises"),
        "arms": results,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"ARTIFACT WRITTEN: {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
