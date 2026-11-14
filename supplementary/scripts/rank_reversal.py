#!/usr/bin/env python3
"""How far the two estimators disagree about the ranking, on the published leaderboards.

The preprint argues that pooled and within-allele AUROC are not two presentations of one
quantity but two different measurements. This quantifies that on the released rows: rank
every entry by each estimator and report how far entries move.

Standard library only. Reads the published leaderboard CSVs, so it measures exactly what
a reader sees.

Needs nothing beyond the two shipped leaderboard CSVs and a standard-library interpreter.

    python3 -I supplementary/scripts/rank_reversal.py [--emit-pgfplots DIR]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys
from pathlib import Path

# The released package root: `supplementary/scripts/x.py` -> two levels up.
ROOT = Path(__file__).resolve().parents[2]

SCRIPT_REL = "supplementary/scripts/rank_reversal.py"
LEADERBOARDS = ROOT / "leaderboards"
OUT = ROOT / "supplementary" / "artifacts" / "rank_reversal.json"
ARMS = ("infectious", "neoantigen")


def num(x: str | None) -> float | None:
    """A withheld figure is not a number and must not sort as one."""
    try:
        return float(x)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def analyse(arm: str) -> dict:
    path = LEADERBOARDS / f"{arm}.csv"
    if not path.is_file():
        print(f"could not run: {path} is missing", file=sys.stderr)
        raise SystemExit(2)
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    ranked = [r for r in rows
              if num(r.get("within_allele")) is not None
              and num(r.get("pooled")) is not None]
    withheld = len(rows) - len(ranked)

    by_wa = sorted(ranked, key=lambda r: -num(r["within_allele"]))
    by_pool = sorted(ranked, key=lambda r: -num(r["pooled"]))
    wa_rank = {r["entry"]: i + 1 for i, r in enumerate(by_wa)}
    pool_rank = {r["entry"]: i + 1 for i, r in enumerate(by_pool)}

    shifts = {e: pool_rank[e] - wa_rank[e] for e in wa_rank}
    absshift = {e: abs(v) for e, v in shifts.items()}
    top = sorted(absshift, key=lambda e: -absshift[e])[:5]

    return {
        "arm": arm,
        "n_entries_total": len(rows),
        "n_ranked": len(ranked),
        "n_withheld": withheld,
        "max_abs_shift": max(absshift.values()),
        "mean_abs_shift": round(sum(absshift.values()) / len(absshift), 2),
        "median_abs_shift": sorted(absshift.values())[len(absshift) // 2],
        "n_moving_5_or_more": sum(1 for v in absshift.values() if v >= 5),
        "n_moving_10_or_more": sum(1 for v in absshift.values() if v >= 10),
        "largest_moves": [
            {"entry": e, "within_allele_rank": wa_rank[e],
             "pooled_rank": pool_rank[e], "shift": shifts[e],
             }
            for e in top
        ],
        "points": [
            {"entry": r["entry"], "within_allele": num(r["within_allele"]),
             "pooled": num(r["pooled"]),
             "baseline": r["entry"].startswith("BASELINE")}
            for r in ranked
        ],
    }


def emit_pgfplots(res: list[dict], dest: Path) -> None:
    """One whitespace table per arm, for pgfplots to read with \\addplot table."""
    dest.mkdir(parents=True, exist_ok=True)
    for r in res:
        lines = ["within pooled baseline"]
        for p in r["points"]:
            lines.append(f"{p['within_allele']:.4f} {p['pooled']:.4f} "
                         f"{1 if p['baseline'] else 0}")
        f = dest / f"reversal-{r['arm']}.dat"
        f.write_text("\n".join(lines) + "\n")
        print(f"  wrote {f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--emit-pgfplots", type=Path, default=None,
                    help="directory to write plot data into")
    args = ap.parse_args(argv)

    res = [analyse(a) for a in ARMS]
    rec = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": SCRIPT_REL,
        "source": "leaderboards",
        "arms": res,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=2) + "\n")

    for r in res:
        print(f"{r['arm']}: {r['n_ranked']} ranked ({r['n_withheld']} withheld)")
        print(f"  max |shift| {r['max_abs_shift']}   mean {r['mean_abs_shift']}   "
              f"median {r['median_abs_shift']}")
        print(f"  moving >=5 places: {r['n_moving_5_or_more']}   "
              f">=10: {r['n_moving_10_or_more']}")
        for m in r["largest_moves"][:3]:
            print(f"    {m['entry'][:48]:50s} #{m['within_allele_rank']:2d} -> "
                  f"#{m['pooled_rank']:2d}")
        print()

    if args.emit_pgfplots:
        emit_pgfplots(res, args.emit_pgfplots)
    print(f"ARTIFACT WRITTEN: {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
