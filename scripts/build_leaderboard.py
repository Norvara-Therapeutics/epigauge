#!/usr/bin/env python3
"""Regenerate the pages and CSVs under leaderboards/ from the reference rows plus the ledger.

Closes the gap between merging a submission and publishing it. Merging appended a row to
`submissions/ledger.csv` and nothing rewrote the published table, so a submitter saw a
scorecard on their PR, a row in the ledger, and no change on the site.

`reference.csv` holds the shipped columns, which change only when the benchmark does.
`submissions/ledger.csv` holds every scored submission. This merges them, sorts on
within-allele, and writes both outputs. The markdown prose between the markers is
preserved verbatim, so editing the narrative and regenerating the table are independent.

    python3 -I scripts/build_leaderboard.py          # rewrite both files
    python3 -I scripts/build_leaderboard.py --check  # exit 1 if they are stale

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "submissions" / "ledger.csv"

#: One entry per published arm. Each has its own reference table, page and CSV; the
#: ledger is shared and filtered by the `benchmark` column a submission declares.
ARMS = {
    "infectious": {
        "reference": ROOT / "leaderboards/sources/infectious.csv",
        "csv": ROOT / "leaderboards/infectious.csv",
        "md": ROOT / "leaderboards/infectious.md",
        "set_id": "titanian_infectious",
    },
    "neoantigen": {
        "reference": ROOT / "leaderboards/sources/neoantigen.csv",
        "csv": ROOT / "leaderboards/neoantigen.csv",
        "md": ROOT / "leaderboards/neoantigen.md",
        "set_id": "titanian_neoantigen",
    },
}

BEGIN = "<!-- BEGIN GENERATED TABLE -->"
END = "<!-- END GENERATED TABLE -->"
COLUMNS = ["rank", "within_allele", "pooled", "entry", "flag", "submitted_by", "date"]


def _f(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def load_reference(arm: dict) -> list[dict]:
    ref = arm["reference"]
    if not ref.is_file():
        raise SystemExit(f"could not run: {ref} is missing; it is the shipped reference "
                         "table and nothing can be rebuilt without it")
    with ref.open(newline="") as fh:
        return list(csv.DictReader(fh))


def load_ledger(arm: dict) -> list[dict]:
    """Submissions scored by CI for THIS arm. Absent before the first merge."""
    if not LEDGER.is_file():
        return []
    rows = []
    with LEDGER.open(newline="") as fh:
        for r in csv.DictReader(fh):
            if (r.get("benchmark") or "").strip() != arm["set_id"]:
                continue  # a submission belongs to exactly one arm
            wa = _f(r.get("within_allele", ""))
            if wa is None:
                continue  # a submission that could not be scored is not a leaderboard row
            name = f"{r.get('method', '').strip()} {r.get('version', '').strip()}".strip()
            rows.append({
                "rank": "", "within_allele": f"{wa:.4f}",
                "pooled": f"{_f(r.get('pooled', '')) or 0:.4f}",
                "entry": name or r.get("team", "unnamed"),
                "flag": "", "submitted_by": r.get("team", ""),
                "date": (r.get("scored_utc", "") or "")[:10],
            })
    return rows


def merged(arm: dict) -> list[dict]:
    """Reference columns and submissions in one table, ranked on within-allele.

    Withheld rows carry no figure and cannot be sorted, so they keep the rank recorded in
    the reference table and are appended in that order.
    """
    rows = load_reference(arm) + load_ledger(arm)
    rankable = [r for r in rows if _f(r.get("within_allele", "")) is not None]
    withheld = [r for r in rows if _f(r.get("within_allele", "")) is None]
    rankable.sort(key=lambda r: float(r["within_allele"]), reverse=True)
    for i, r in enumerate(rankable, 1):
        r["rank"] = str(i)
    return rankable + withheld


def to_csv(rows: list[dict]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLUMNS, lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in COLUMNS})
    return buf.getvalue()


def to_markdown(rows: list[dict]) -> str:
    out = ["| # | within-allele | pooled | entry |", "|---|---|---|---|"]
    for r in rows:
        entry = r["entry"]
        if r.get("flag", "").startswith("not-a-capability"):
            entry += " ⚠"
        if entry.startswith("BASELINE"):
            entry = f"**{entry}**"
            cells = (f"**{r['rank']}**", f"**{r['within_allele']}**", r["pooled"])
        else:
            cells = (r["rank"], r["within_allele"], r["pooled"])
        out.append(f"| {cells[0]} | {cells[1]} | {cells[2]} | {entry} |")
    return "\n".join(out)


def splice(md: str, table: str, name: str) -> str:
    if BEGIN not in md or END not in md:
        raise SystemExit(f"could not run: {name} has no {BEGIN} / {END} markers, so the "
                         "generated table has nowhere to go and the prose would be overwritten")
    head, rest = md.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    # Blank lines are REQUIRED on both sides. kramdown treats an HTML comment as a
    # block-level element and absorbs what follows as raw HTML, so without them the
    # table is never parsed as markdown: GitHub renders it, GitHub Pages does not,
    # and `---` falls through to typographic substitution as an em-dash.
    return f"{head}{BEGIN}\n\n{table}\n\n{END}{tail}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", choices=sorted(ARMS), action="append",
                    help="arm to rebuild; repeatable, default all")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any published file is stale; write nothing")
    args = ap.parse_args(argv)

    stale, wrote = [], []
    for name in (args.arm or sorted(ARMS)):
        arm = ARMS[name]
        rows = merged(arm)
        new_csv = to_csv(rows)
        new_md = splice(arm["md"].read_text(), to_markdown(rows), arm["md"].name)
        if args.check:
            stale += [p.name for p, new in ((arm["csv"], new_csv), (arm["md"], new_md))
                      if p.read_text() != new]
        else:
            arm["csv"].write_text(new_csv)
            arm["md"].write_text(new_md)
            from_ledger = len(load_ledger(arm))
            wrote.append(f"{name}: {len(rows)} rows, {from_ledger} from the ledger")

    if args.check:
        if stale:
            print(f"STALE: {', '.join(stale)} do not match the reference tables and the ledger.")
            print("Run scripts/build_leaderboard.py and commit the result.")
            return 1
        print("current: every published file matches its sources.")
        return 0

    for line in wrote:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
