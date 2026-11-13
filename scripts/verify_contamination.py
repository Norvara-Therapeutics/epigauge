#!/usr/bin/env python3
"""Verify, from the two files themselves, that a public dataset's `im_test.csv` and a
published 833-row neoantigen benchmark are the same evaluation table.

Run it. It is a check, so it exits non-zero when the contamination is present.

    python3 scripts/verify_contamination.py
    python3 scripts/verify_contamination.py --a /your/im_test.csv --b /your/TITANiAN_IM_Neoantigen.csv

WHY THIS SCRIPT EXISTS AS A SCRIPT. A contamination finding stated in prose is an
accusation and the reader has to take it on trust. The same finding as a command the
reader runs on their own copies is an exercise, and the reader's own machine answers.
Everything printed below is computed from the two files at run time; nothing is
hard-coded except the names of the columns the files themselves declare.

WHAT IT CHECKS, IN THE ORDER A READER SHOULD SEE IT.

  Step 1  The NAIVE join, on the raw allele strings exactly as each file spells them.
          It returns 0 of 833. This step is not padding -- it is the finding. One file
          spells the allele `HLA-A*01:01` and the other spells it `A0101`, so the
          obvious check reports a clean-looking zero and the reader stops. The naive
          result is printed FIRST so the reader sees the trap before the answer.

  Step 2  The same join after `epigauge.canon.canonicalise_allele`, which is this
          project's single allele normaliser and the fixture for its own documented
          instrument-failure of exactly this kind.

  Step 3  Label agreement on the matched pairs (`tgt` in one file, `label` in the other).

  Step 4  Every predictor column shared by the two headers, compared BYTE FOR BYTE as
          the CSV text stands -- not as parsed floats, not within a tolerance. Two
          independently produced score tables do not agree to the last printed digit.

No third-party dependencies. Python 3.8+. Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

# scripts/ sits beside epigauge/, so the package resolves without installation.
_PKG_PARENT = Path(__file__).resolve().parents[1]
if str(_PKG_PARENT) not in sys.path:
    sys.path.insert(0, str(_PKG_PARENT))

try:
    from epigauge.canon import CanonError, canonicalise_allele, canonicalise_peptide
except ImportError as exc:  # pragma: no cover - environment problem, not a data finding
    sys.stderr.write(
        f"cannot import epigauge.canon from {_PKG_PARENT}: {exc}\n"
        "Run this script from a checkout that contains the epigauge/ package.\n"
    )
    raise SystemExit(2)


# Each file names its own columns differently. Candidates are tried in order; the
# script refuses rather than guessing if none is present.
ALLELE_COLS = ("mhc", "allele", "Allele", "HLA", "hla")
PEPTIDE_COLS = ("pep", "peptide", "Peptide", "epitope")
LABEL_COLS = ("tgt", "label", "Label", "target", "immunogenicity")

# The seven predictors named in the written finding. Each is a prefix; a predictor may
# contribute more than one column (`_Scores` and `_Ranks`).
NAMED_SEVEN = (
    "NetMHCpan-4.1",
    "MHCflurry-2.0",
    "PRIME-2.0",
    "BigMHC_IM",
    "BigMHC_EL",
    "TransPHLA",
    "HLAthena",
)

RULE = "=" * 78
THIN = "-" * 78


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pick(header, candidates, what, path):
    for c in candidates:
        if c in header:
            return c
    raise SystemExit(
        f"{path}: no {what} column found. Tried {list(candidates)}; header is {header}"
    )


def load(path: Path):
    """Return (header, rows) with every field left as the exact text in the file."""
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"{path}: file is empty")
        rows = [r for r in reader if r and any(f.strip() for f in r)]
    return header, rows


def as_dicts(header, rows):
    return [dict(zip(header, r)) for r in rows]


def canon_key(allele_raw, peptide_raw):
    """Canonical (allele, peptide). Returns None if either field is unparseable.

    `canon` raises rather than coercing, by design: a plausible fallback key is how a
    mismatch becomes a published number. Here a raise means the row is reported as
    uncanonicalisable, never silently folded into some default bucket.
    """
    try:
        return (canonicalise_allele(allele_raw), canonicalise_peptide(peptide_raw))
    except CanonError:
        return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Verify whether two immunogenicity evaluation tables are the same table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--a", type=Path, required=True,
                    help="the larger/community-shipped file, e.g. bigmhc im_test.csv")
    ap.add_argument("--b", type=Path, required=True,
                    help="the published benchmark file (default: TITANiAN_IM_Neoantigen.csv)")
    ap.add_argument("--show", type=int, default=3,
                    help="how many matched rows to print in full as a spot check (default: 3)")
    args = ap.parse_args(argv)

    for p in (args.a, args.b):
        if not p.is_file():
            sys.stderr.write(f"not a file: {p}\n")
            return 2

    print(RULE)
    print("CONTAMINATION CHECK: are these two files the same evaluation table?")
    print(RULE)

    # ---------------------------------------------------------------- inputs
    header_a, rows_a = load(args.a)
    header_b, rows_b = load(args.b)
    dicts_a = as_dicts(header_a, rows_a)
    dicts_b = as_dicts(header_b, rows_b)

    print()
    print("INPUTS")
    for tag, path, header, rows in (("A", args.a, header_a, rows_a),
                                    ("B", args.b, header_b, rows_b)):
        print(f"  {tag}  {path}")
        print(f"     sha256   {sha256_of(path)}")
        print(f"     bytes    {path.stat().st_size}")
        print(f"     rows     {len(rows)} data rows, {len(header)} columns")

    al_a = pick(header_a, ALLELE_COLS, "allele", args.a)
    pe_a = pick(header_a, PEPTIDE_COLS, "peptide", args.a)
    lb_a = pick(header_a, LABEL_COLS, "label", args.a)
    al_b = pick(header_b, ALLELE_COLS, "allele", args.b)
    pe_b = pick(header_b, PEPTIDE_COLS, "peptide", args.b)
    lb_b = pick(header_b, LABEL_COLS, "label", args.b)

    print()
    print(f"  key columns   A: ({al_a}, {pe_a})  label: {lb_a}")
    print(f"                B: ({al_b}, {pe_b})  label: {lb_b}")
    print(f"  A spells its first allele   {dicts_a[0][al_a]!r}")
    print(f"  B spells its first allele   {dicts_b[0][al_b]!r}")

    n_b = len(dicts_b)

    # ------------------------------------------------- step 1: the naive join
    print()
    print(THIN)
    print("STEP 1  NAIVE JOIN on the raw allele strings, as most reviewers would run it")
    print(THIN)
    raw_a = {(r[al_a], r[pe_a]) for r in dicts_a}
    naive_hits = sum(1 for r in dicts_b if (r[al_b], r[pe_b]) in raw_a)
    print(f"  (allele, peptide) pairs of B found in A:  {naive_hits} / {n_b}")
    if naive_hits == 0:
        print()
        print("  ^ ZERO. On this evidence the two files share nothing and the benchmark")
        print("    is independent. That conclusion is wrong, and the next step shows why.")
        print("    The allele field is spelled differently in the two files, so the join")
        print("    key never matches. The check did not fail loudly; it failed clean.")

    # ------------------------------------------- step 2: canonicalised join
    print()
    print(THIN)
    print("STEP 2  SAME JOIN after epigauge.canon.canonicalise_allele")
    print(THIN)
    print(f"  canonicalise_allele({dicts_a[0][al_a]!r}) -> "
          f"{canonicalise_allele(dicts_a[0][al_a])!r}")
    print(f"  canonicalise_allele({dicts_b[0][al_b]!r}) -> "
          f"{canonicalise_allele(dicts_b[0][al_b])!r}")

    index_a = {}
    bad_a = 0
    dup_a = 0
    for r in dicts_a:
        k = canon_key(r[al_a], r[pe_a])
        if k is None:
            bad_a += 1
            continue
        if k in index_a:
            dup_a += 1
        index_a[k] = r

    matched = []          # (b_row, a_row)
    unmatched_b = []
    bad_b = 0
    seen_b = set()
    dup_b = 0
    for r in dicts_b:
        k = canon_key(r[al_b], r[pe_b])
        if k is None:
            bad_b += 1
            continue
        if k in seen_b:
            dup_b += 1
        seen_b.add(k)
        a_row = index_a.get(k)
        if a_row is None:
            unmatched_b.append(r)
        else:
            matched.append((r, a_row))

    n_match = len(matched)
    print()
    print(f"  (allele, peptide) pairs of B found in A:  {n_match} / {n_b}")
    print(f"  rows of B with no counterpart in A:       {len(unmatched_b)}")
    print(f"  rows of A not in B (A is the superset):   {len(index_a) - n_match}")
    if bad_a or bad_b:
        print(f"  rows dropped as uncanonicalisable:        A {bad_a}, B {bad_b}")
    if dup_a or dup_b:
        print(f"  duplicate canonical keys within a file:   A {dup_a}, B {dup_b}")

    # ----------------------------------------------- step 3: label agreement
    print()
    print(THIN)
    print("STEP 3  LABEL AGREEMENT on the matched pairs")
    print(THIN)
    label_agree = sum(1 for rb, ra in matched if rb[lb_b].strip() == ra[lb_a].strip())
    label_disagree = [(rb, ra) for rb, ra in matched
                      if rb[lb_b].strip() != ra[lb_a].strip()]
    pos_b = sum(1 for rb, _ in matched if rb[lb_b].strip() == "1")
    print(f"  B.{lb_b} == A.{lb_a}:  {label_agree} / {n_match}"
          + (f"  ({100.0 * label_agree / n_match:.3f}%)" if n_match else ""))
    print(f"  positives among the matched rows:  {pos_b}")
    for rb, ra in label_disagree[:5]:
        print(f"    disagree: {rb[al_b]} {rb[pe_b]}  B={rb[lb_b]!r} A={ra[lb_a]!r}")

    # ------------------------------- step 4: predictor columns, byte for byte
    print()
    print(THIN)
    print("STEP 4  SHARED PREDICTOR COLUMNS, compared BYTE FOR BYTE (not as floats)")
    print(THIN)
    key_cols = {al_a, pe_a, lb_a, al_b, pe_b, lb_b}
    shared_cols = [c for c in header_a if c in header_b and c not in key_cols]
    order_a = [c for c in header_a if c in shared_cols]
    order_b = [c for c in header_b if c in shared_cols]

    print(f"  columns present in BOTH headers (excluding key columns):  {len(shared_cols)}")
    print(f"  and they appear in the SAME ORDER in both headers:        "
          f"{'yes' if order_a == order_b else 'NO'}")
    print()

    if n_match == 0:
        print("  No matched rows, so there is nothing to compare column-wise. A column")
        print("  reported as 'identical' over zero rows would be vacuous, so this step")
        print("  reports only the header overlap above.")

    identical_cols = []
    differing_cols = []
    for col in shared_cols:
        diffs = sum(1 for rb, ra in matched if rb[col] != ra[col])
        if diffs == 0:
            identical_cols.append(col)
        else:
            differing_cols.append((col, diffs))

    if n_match:
        width = max((len(c) for c in shared_cols), default=10)
        for col in shared_cols:
            d = dict(differing_cols).get(col, 0)
            mark = "identical" if d == 0 else f"{d} differing rows"
            print(f"    {col:<{width}}  {mark}")

        print()
        print(f"  byte-identical across all {n_match} matched rows:  "
              f"{len(identical_cols)} / {len(shared_cols)} shared predictor columns")

        # The seven predictors named in the written finding, checked by name.
        print()
        print("  the seven predictors named in CONTAMINATION.md, by name:")
        for name in NAMED_SEVEN:
            cols = [c for c in shared_cols if c == name or c.startswith(name + "_")]
            if not cols:
                print(f"    {name:<16}  NOT PRESENT in both headers")
                continue
            ok = all(c in identical_cols for c in cols)
            print(f"    {name:<16}  {'byte-identical' if ok else 'DIFFERS'}"
                  f"   ({len(cols)} column{'s' if len(cols) != 1 else ''}: {', '.join(cols)})")

    # ----------------------------------------------------------- spot check
    if args.show > 0 and matched:
        print()
        print(THIN)
        print(f"SPOT CHECK  first {min(args.show, len(matched))} matched rows, shared fields verbatim")
        print(THIN)
        for rb, ra in matched[: args.show]:
            print(f"  B: {rb[al_b]:<12} {rb[pe_b]:<12} label={rb[lb_b]}")
            print(f"  A: {ra[al_a]:<12} {ra[pe_a]:<12} tgt  ={ra[lb_a]}")
            sample = shared_cols[:4]
            for c in sample:
                same = "==" if rb[c] == ra[c] else "!="
                print(f"       {c:<24} B {rb[c]:<22} {same} A {ra[c]}")
            print()

    # -------------------------------------------------------------- verdict
    contaminated = (
        n_match > 0
        and n_match == n_b
        and label_agree == n_match
        and len(identical_cols) == len(shared_cols)
        and len(shared_cols) > 0
    )
    partial = n_match > 0 and not contaminated

    print(RULE)
    print("VERDICT")
    print(RULE)
    print(f"  naive join                 {naive_hits} / {n_b}")
    print(f"  canonicalised join         {n_match} / {n_b}")
    print(f"  label agreement            {label_agree} / {n_match}")
    print(f"  byte-identical predictors  "
          f"{len(identical_cols) if n_match else 0} / {len(shared_cols)}"
          + ("" if n_match else "  (vacuous: no matched rows)"))
    print(f"  extra rows in A            {len(index_a) - n_match}")
    print()
    if contaminated:
        print("  Every row of B is in A, with the same label, and every shared predictor")
        print("  column agrees to the last printed character. These are not two datasets")
        print("  that overlap. B is a subset of A with extra columns appended.")
        print()
        print("  WHAT THIS SUPPORTS: a score on B is not evidence of capability on data")
        print("  the scorer has not seen, because the table is public.")
        print("  WHAT THIS DOES NOT ESTABLISH: anything about how any file came to exist,")
        print("  or about anyone's intent. This script compares files. It cannot and does")
        print("  not measure what any model was trained on.")
        print()
        print("  CHECK FAILED (exit 1): the contamination is present.")
        return 1
    if partial:
        print("  Partial overlap. The two files share rows but are not the same table on")
        print("  every axis checked. Read the per-column report above before describing")
        print("  this as contamination; the strong claim needs all four lines to agree.")
        print()
        print("  CHECK FAILED (exit 1): overlap is present.")
        return 1
    print("  No overlap after canonicalisation. These files are independent on the")
    print("  evidence of this check.")
    print()
    print("  CHECK PASSED (exit 0).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
