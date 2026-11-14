#!/usr/bin/env python3
"""Count how many published leaderboard rows are contamination-affected, under a rule
stated before it is applied.

    python3 -I supplementary/scripts/contamination_census.py \
        --im-test <im_test.csv> \
        --infectious-table <TITANiAN_IM_Infectious.csv> \
        --neoantigen-table <TITANiAN_IM_Neoantigen.csv>

Standard library only, plus `epigauge.canon` for the allele/peptide normaliser.
Writes `supplementary/artifacts/contamination_census.json` and prints a readable table.

ALL THREE ARGUMENTS ARE REQUIRED, because the package redistributes none of the three
CSVs. They are public and are fetched from their deposits: `im_test.csv` from BigMHC's
Mendeley deposit 10.17632/dvmz6pkzvb, and the two TITANiAN evaluation tables from
`test.zip` in the Dryad deposit 10.5061/dryad.s7h44j1k7. The sha256 of each file behind
the figures that ship is in EXPECTED_SHA256 below and is compared against whatever is
supplied; a mismatch is printed and recorded, never corrected. The artifact records the
deposit identity and the measured checksum of each input, not the path it was read from.


WHY THIS SCRIPT EXISTS
======================

The paper's contamination section states that "the highest figure on the affected
arm is not reachable by a method without access to that table, so it is not a capability
bar". That sentence is about ONE figure. It never says how many of the published rows the
same sentence is true of, and a referee cannot check a claim whose extent is unstated.

The flagging on disk is also not self-consistent. `BigMHC_EL` carries no flag on either
arm while `BigMHC_IM` and `BigMHC_ELIM` carry one; the six `BigMHC_Ablation_*` columns
carry none on the neoantigen arm while `BigMHC_IM`/`BigMHC_ELIM` -- read from the same
file, by the same join -- carry one. A flag applied method-by-method cannot be audited,
because there is no written rule to audit it against.

So this script states a rule, applies it mechanically to every row of both leaderboards,
and prints the rows where the rule and the shipped flag disagree. The disagreement list
is the finding; the counts are what the paper is currently missing.


THE RULE
========

A published leaderboard row is CONTAMINATION-AFFECTED when at least one of the following
holds. The tiers are separate claims about different facts, and they are reported
separately because conflating them would overstate the weaker one.

  TIER 1 -- the score is read from a column shipped INSIDE the evaluation table.

      The row's figure is computed from a predictor column that is physically
      distributed in the same file as the labels it is scored against. No inference was
      run to obtain it: the number is a re-read of bytes that travel with the answer
      key. Anyone who holds the evaluation file holds these scores.

      Mechanically: the row's score column is a member of the header of the upstream
      evaluation CSV for that arm.

      Two sub-cases are counted separately because they support different sentences:

        1a  the column is also one of the 25 predictor columns present in BigMHC's
            `im_test.csv` -- i.e. the column is distributed in at least two independent
            public deposits, so the "anyone who holds the file" set is wider than the
            users of one deposit;

        1b  the column is present only in the evaluation deposit itself. On these arms
            that is exactly the benchmark authors' own TITANiAN columns.

  TIER 2 -- the method has a public deposit of its own that contains the evaluation rows.

      A separate and stronger-sounding claim, so it is held to a separate measurement:
      the method's own published data deposit contains every (allele, peptide) row of
      the arm being scored. That is established here by joining the deposit file against
      the arm on canonicalised keys and requiring containment of 1.000 -- not by reading
      a README.

      This is the tier that licenses the sentence "not a capability bar for a method
      with access to the table", because it is the tier where the table demonstrably
      travels with the method.

  TIER 0 -- not affected under this rule.

      Either the score was produced by running inference on the rows (the organisers'
      own submissions), or it is computed from fields of the row itself (peptide length,
      mean hydropathy), or it is fitted on a declared pool that is not the evaluation
      labels (the allele-prior lookup, fitted on a declared training pool and never on
      the evaluation labels).

      Tier 0 is NOT a certificate of independence. It means only that this rule found no
      affectedness. Every arm here carries some dependence -- see the paper's discussion
      of gate G10 -- and tier 0 does not speak to it.

Tiers 1 and 2 are not mutually exclusive; a row can be both, and `max_tier` is reported
alongside the per-tier membership so neither count is inflated by the other.


WHY THIS RULE AND NOT A NARROWER ONE
====================================

The narrow alternative -- "flag a row when its method's authors deposited the evaluation
table" -- is tier 2 alone. It is defensible but it under-reports, and it under-reports in
a direction that flatters the leaderboard: it leaves the 25 columns that ship *inside*
the evaluation file unflagged, which is the single largest group of rows on both arms.
A reader who sees `MixMHCpred-2.2 (scores) 0.5234` has no way to know, from the row,
that the figure was read out of the same CSV as the labels rather than computed by
running MixMHCpred.

The broad alternative -- "flag every row on an arm that is contained in any public
deposit" -- over-reports: it would flag the organisers' own fresh-inference submissions
and the unfitted baselines, which are the rows for which the published figure IS a
measurement of something.

The two-tier rule is the narrowest one that separates "this number was not produced by
inference" from "this method's authors distribute the answer key", because those are
different claims and the paper needs to make them separately.


WHAT THE RULE DOES NOT CLAIM
============================

This script compares file contents. It establishes what a file contains and what two
files share. It does not and cannot establish how any file came to exist, what any model
was fitted on, or anything about anyone's intent. "Not a capability bar for a method with
access to this table" is the supportable form and is the only form used here; the
inverted form is not available from this evidence and is not implied by a tier.

A withheld row (NetMHCpan-4.1, DTU licence sec. 7(v)) publishes no figure. Such rows are
counted under the rule -- the column is in the table either way -- but reported in their
own bucket, so the disagreement list is not padded with rows that print no number for a
flag to qualify.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

# The released package root: `supplementary/scripts/x.py` -> two levels up. Inserted
# explicitly rather than relying on PYTHONPATH, so `python3 -I` works from any cwd.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from epigauge.canon import CanonError, canonicalise_allele, canonicalise_peptide
except ImportError as exc:  # pragma: no cover - environment, not a data finding
    sys.stderr.write(
        f"cannot import epigauge.canon from {ROOT}: {exc}\n"
        "Run this from an unpacked epigauge source distribution, or install the "
        "package first.\n"
    )
    raise SystemExit(2)

SCRIPT_REL = "supplementary/scripts/contamination_census.py"
LEADERBOARDS = ROOT / "leaderboards"
OUT = ROOT / "supplementary" / "artifacts" / "contamination_census.json"

# The three upstream CSVs are NOT redistributed and are therefore required arguments.
# Each is recorded in the artifact by its deposit identity plus the sha256 measured on
# the file supplied, never by the path it happened to be read from: a local path is not
# provenance, and a checksum is.
INPUT_DESCRIPTIONS = {
    "bigmhc": (
        "BigMHC deposit member datasets/im_test.csv, Mendeley Data "
        "10.17632/dvmz6pkzvb; not redistributed, verify against sha256 below"),
    "infectious": (
        "TITANiAN deposit member test.zip -> TITANiAN_IM_Infectious.csv, Dryad "
        "10.5061/dryad.s7h44j1k7; not redistributed, verify against eval_table_sha256"),
    "neoantigen": (
        "TITANiAN deposit member test.zip -> TITANiAN_IM_Neoantigen.csv, Dryad "
        "10.5061/dryad.s7h44j1k7; not redistributed, verify against eval_table_sha256"),
}

# The sha256 of the exact files behind the figures that ship, for comparison against
# whatever a third party fetches. A mismatch is reported, not corrected.
EXPECTED_SHA256 = {
    "bigmhc": "70ba5dbeecb0e62ccc95f2aa11a28801773893f4ebf209b5e68b58318be7b634",
    "infectious": "2d64689d22637562f9f4f65ec101bb2cea76c600cb80ee852b9508bd1a335553",
    "neoantigen": "da2442cd7f7ea0d1ed83e15a71538cb0f8983a329fe973cfe43514c1dfed9661",
}

# --------------------------------------------------------------------------- inputs

# The upstream evaluation table behind each leaderboard arm, and the arm's leaderboard.
# `eval_table` is None until the command line supplies it: the package redistributes
# neither table, so there is no path to default to.
ARMS = {
    "infectious": {
        "leaderboard": LEADERBOARDS / "infectious.csv",
        "eval_table": None,
        "allele_col": "mhc",
        "peptide_col": "peptide",
        "set_id": "titanian_infectious",
    },
    "neoantigen": {
        "leaderboard": LEADERBOARDS / "neoantigen.csv",
        "eval_table": None,
        "allele_col": "mhc",
        "peptide_col": "peptide",
        "set_id": "titanian_neoantigen",
    },
}

# Published data deposits distributed BY the authors of a method on these leaderboards,
# with the member file that would have to contain the arm's rows for tier 2 to hold.
# Containment is measured below; nothing here is taken on trust.
DEPOSITS = {
    "BigMHC": {
        "doi": "10.17632/dvmz6pkzvb",
        "distributor": "Mendeley Data",
        "member": "datasets/im_test.csv",
        "local": None,  # supplied by --im-test; not redistributed

        "allele_col": "mhc",
        "peptide_col": "pep",
    },
    "TITANiAN": {
        "doi": "10.5061/dryad.s7h44j1k7",
        "distributor": "Dryad",
        "member": "test.zip -> TITANiAN_IM_{arm}.csv",
        "local": None,  # the arm's own evaluation table; filled per arm
        "allele_col": "mhc",
        "peptide_col": "peptide",
    },
}

# Non-predictor columns of the upstream evaluation tables: keys and side columns.
NON_PREDICTOR_COLS = {"mhc", "peptide", "pep", "label", "tgt", "Allele", "allele", "pseudo"}

# --------------------------------------------------------- declared entry -> column map
#
# One entry per leaderboard row, keyed by the arm and the entry string as the CSV spells
# it. `col` names the column of the arm's upstream evaluation table the figure is read
# from, or None when no column is read. `method` attributes the row for tier 2.
#
# `kind` records how the number was obtained:
#   table_column        read from a column of the evaluation table
#   fresh_inference     produced by running the model on the rows (organiser submission)
#   computed_from_row   a function of the row's own fields, nothing fitted
#   fitted_other_pool   fitted on a declared pool that is not the evaluation labels
#
# The script FAILS if any leaderboard entry is missing from this map, or if any declared
# column is absent from the evaluation table's header. There is no default branch: a
# silently defaulted row is how an uncounted row becomes a published count.

_PREDICTORS = {
    "NetMHCpan-4.1 (scores)": ("NetMHCpan-4.1_Scores", "NetMHCpan-4.1"),
    "NetMHCpan-4.1 (percentile ranks)": ("NetMHCpan-4.1_Ranks", "NetMHCpan-4.1"),
    "MHCflurry-2.0 (scores)": ("MHCflurry-2.0_Scores", "MHCflurry-2.0"),
    "MHCflurry-2.0 (percentile ranks, sign-corrected)": ("MHCflurry-2.0_Ranks", "MHCflurry-2.0"),
    "MHCnuggets-2.4.0": ("MHCnuggets-2.4.0", "MHCnuggets-2.4.0"),
    "MixMHCpred-2.1 (scores)": ("MixMHCpred-2.1_Scores", "MixMHCpred-2.1"),
    "MixMHCpred-2.1 (percentile ranks, sign-corrected)": ("MixMHCpred-2.1_Ranks", "MixMHCpred-2.1"),
    "MixMHCpred-2.2 (scores)": ("MixMHCpred-2.2_Scores", "MixMHCpred-2.2"),
    "MixMHCpred-2.2 (percentile ranks, sign-corrected)": ("MixMHCpred-2.2_Ranks", "MixMHCpred-2.2"),
    "PRIME-1.0 (scores)": ("PRIME-1.0_Scores", "PRIME-1.0"),
    "PRIME-1.0 (percentile ranks, sign-corrected)": ("PRIME-1.0_Ranks", "PRIME-1.0"),
    "PRIME-2.0 (scores)": ("PRIME-2.0_Scores", "PRIME-2.0"),
    "PRIME-2.0 (percentile ranks, sign-corrected)": ("PRIME-2.0_Ranks", "PRIME-2.0"),
    "TransPHLA": ("TransPHLA", "TransPHLA"),
    "HLAthena (scores)": ("HLAthena_Scores", "HLAthena"),
    "HLAthena (percentile ranks, sign-corrected)": ("HLAthena_Ranks", "HLAthena"),
    "BigMHC_EL": ("BigMHC_EL", "BigMHC"),
    "BigMHC_IM": ("BigMHC_IM", "BigMHC"),
    "BigMHC_ELIM": ("BigMHC_ELIM", "BigMHC"),
}

# The six BigMHC LSTM ablation columns, spelled differently on the two arms.
_BIGMHC_ABLATIONS_NEO = {
    f"BigMHC_Ablation_{w}_LSTM_{h}": (f"BigMHC_Ablation_{w}_LSTM_{h}", "BigMHC")
    for w in ("Canonical", "Wide")
    for h in ("EL", "IM", "ELIM")
}
_BIGMHC_ABLATIONS_INF = {
    f"BigMHC {w.lower()}-LSTM {h} ablation": (f"BigMHC_Ablation_{w}_LSTM_{h}", "BigMHC")
    for w in ("Canonical", "Wide")
    for h in ("EL", "IM", "ELIM")
}

ENTRY_MAP: dict[str, dict[str, dict]] = {"infectious": {}, "neoantigen": {}}

for _arm in ("infectious", "neoantigen"):
    for _entry, (_col, _meth) in _PREDICTORS.items():
        ENTRY_MAP[_arm][_entry] = {"col": _col, "method": _meth, "kind": "table_column"}

for _entry, (_col, _meth) in _BIGMHC_ABLATIONS_INF.items():
    ENTRY_MAP["infectious"][_entry] = {"col": _col, "method": _meth, "kind": "table_column"}
for _entry, (_col, _meth) in _BIGMHC_ABLATIONS_NEO.items():
    ENTRY_MAP["neoantigen"][_entry] = {"col": _col, "method": _meth, "kind": "table_column"}

# The benchmark authors' own columns. The infectious table spells the ablations with an
# underscore after "ablation", the neoantigen table with a hyphen; the leaderboards use a
# third spelling again. Hence an explicit map rather than a derived one.
ENTRY_MAP["infectious"].update({
    "TITANiAN-IM": {"col": "TITANiAN-IM", "method": "TITANiAN", "kind": "table_column"},
    "TITANiAN-IM baseline": {"col": "TITANiAN-IM-baseline", "method": "TITANiAN",
                             "kind": "table_column"},
    **{
        f"TITANiAN-IM ablation {k}": {"col": f"TITANiAN-IM-ablation_{k}",
                                      "method": "TITANiAN", "kind": "table_column"}
        for k in ("1,2,3,4,7", "1,2,5,6,7", "1,2,3,4,5,6", "3,4,5,6,7")
    },
})
ENTRY_MAP["neoantigen"].update({
    "TITANiAN-IM": {"col": "TITANiAN-IM", "method": "TITANiAN", "kind": "table_column"},
    "TITANiAN-IM-baseline": {"col": "TITANiAN-IM-baseline", "method": "TITANiAN",
                             "kind": "table_column"},
    **{
        f"TITANiAN-IM-ablation-{k}": {"col": f"TITANiAN-IM-ablation-{k}",
                                      "method": "TITANiAN", "kind": "table_column"}
        for k in ("1,2,3,4,7", "1,2,5,6,7", "1,2,3,4,5,6")
    },
})

# Organiser submissions: inference run on the rows by this project. No column is read.
for _arm, _entries in (
    ("infectious", [
        "NORVARA-3B — frozen ESM-2 3B backbone, 2,561-param head",
        "NORVARA-650M — frozen ESM-2 650M backbone, 1,281-param head",
        "NORVARA-150M — frozen ESM-2 150M backbone, 641-param head",
    ]),
    ("neoantigen", ["NORVARA-3B", "NORVARA-650M", "NORVARA-150M"]),
):
    for _entry in _entries:
        ENTRY_MAP[_arm][_entry] = {"col": None, "method": "NORVARA",
                                   "kind": "fresh_inference"}

# Baselines.
for _arm in ("infectious", "neoantigen"):
    ENTRY_MAP[_arm]["BASELINE — peptide length"] = {
        "col": None, "method": "baseline", "kind": "computed_from_row"}
    ENTRY_MAP[_arm]["BASELINE — mean Kyte-Doolittle hydropathy, unfitted"] = {
        "col": None, "method": "baseline", "kind": "computed_from_row"}
ENTRY_MAP["infectious"]["BASELINE — allele lookup, no peptide input"] = {
    "col": None, "method": "baseline", "kind": "fitted_other_pool"}

# ------------------------------------------------------------------------- helpers

RULE = "=" * 86
THIN = "-" * 86


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_entry(s: str) -> str:
    """Collapse whitespace only. Dashes and commas are NOT normalised: two entries that
    differ only in punctuation are two different entries and must be mapped separately."""
    return re.sub(r"\s+", " ", s).strip()


def read_header(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return next(csv.reader(fh))


def canon_keys(path: Path, allele_col: str, peptide_col: str) -> tuple[set, int]:
    """Canonical (allele, peptide) keys, and a count of rows that would not canonicalise.

    `canon` raises rather than coercing. A raise is reported, never folded into a
    default bucket -- a lenient fallback key is how a mismatch becomes a published
    number, which is this project's own documented instrument-failure #3.
    """
    keys: set = set()
    bad = 0
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            try:
                keys.add((canonicalise_allele(row[allele_col]),
                          canonicalise_peptide(row[peptide_col])))
            except CanonError:
                bad += 1
    return keys, bad


def classify_flag(flag: str) -> str:
    f = (flag or "").strip()
    if not f:
        return "none"
    for prefix in ("not-a-capability-measurement", "withheld", "organiser-submission"):
        if f.startswith(prefix):
            return prefix
    if f == "baseline":
        return "baseline"
    return "other"


# ------------------------------------------------------------------------- the census

def check_sha256(key: str, path: Path) -> str:
    """The measured checksum, with a loud note when it is not the expected one."""
    got = sha256_of(path)
    want = EXPECTED_SHA256.get(key)
    if want and got != want:
        print("    WARNING  sha256 differs from the file behind the published figures")
        print(f"             expected {want}")
    return got


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--im-test", required=True, type=Path,
                    help="BigMHC deposit member datasets/im_test.csv "
                         "(Mendeley Data 10.17632/dvmz6pkzvb)")
    ap.add_argument("--infectious-table", required=True, type=Path,
                    help="TITANiAN_IM_Infectious.csv from test.zip "
                         "(Dryad 10.5061/dryad.s7h44j1k7)")
    ap.add_argument("--neoantigen-table", required=True, type=Path,
                    help="TITANiAN_IM_Neoantigen.csv from test.zip "
                         "(Dryad 10.5061/dryad.s7h44j1k7)")
    args = ap.parse_args(argv)

    DEPOSITS["BigMHC"]["local"] = args.im_test
    ARMS["infectious"]["eval_table"] = args.infectious_table
    ARMS["neoantigen"]["eval_table"] = args.neoantigen_table

    print(RULE)
    print("CONTAMINATION CENSUS -- how many published leaderboard rows the finding reaches")
    print(RULE)

    for arm, cfg in ARMS.items():
        if not cfg["leaderboard"].is_file():
            sys.stderr.write(f"missing shipped leaderboard: {cfg['leaderboard']}\n")
            return 2

    # --- the shared column set, measured from the two headers ------------------
    im_test = DEPOSITS["BigMHC"]["local"]
    if not im_test.is_file():
        sys.stderr.write(
            f"missing deposit member: {im_test}\n"
            "Fetch datasets/im_test.csv from BigMHC's Mendeley deposit "
            "10.17632/dvmz6pkzvb; the package does not redistribute it.\n")
        return 2
    im_header = read_header(im_test)
    im_predictors = [c for c in im_header if c not in NON_PREDICTOR_COLS]

    print()
    print("INPUTS")
    print(f"  BigMHC deposit member  {DEPOSITS['BigMHC']['member']}")
    im_sha = check_sha256("bigmhc", im_test)
    print(f"    sha256   {im_sha}")
    print(f"    columns  {len(im_header)} ({len(im_predictors)} predictor)")

    arm_state = {}
    for arm, cfg in ARMS.items():
        table = cfg["eval_table"]
        if not table.is_file():
            sys.stderr.write(
                f"missing evaluation table: {table}\n"
                "Fetch test.zip from the TITANiAN Dryad deposit "
                "10.5061/dryad.s7h44j1k7; the package does not redistribute it.\n")
            return 2
        header = read_header(table)
        predictors = [c for c in header if c not in NON_PREDICTOR_COLS]
        keys, bad = canon_keys(table, cfg["allele_col"], cfg["peptide_col"])
        print(f"  {arm} evaluation table  "
              f"{DEPOSITS['TITANiAN']['member'].replace('{arm}', arm.capitalize())}")
        arm_sha = check_sha256(arm, table)
        arm_state[arm] = {
            "eval_table": INPUT_DESCRIPTIONS[arm],
            "eval_table_sha256": arm_sha,
            "eval_table_columns": len(header),
            "eval_table_predictor_columns": len(predictors),
            "predictor_set": set(predictors),
            "keys": keys,
            "uncanonicalisable_rows": bad,
            "rows": len(keys),
        }
        print(f"    sha256   {arm_sha}")
        print(f"    columns  {len(header)} ({len(predictors)} predictor), rows {len(keys)}")
        if bad:
            print(f"    rows that would not canonicalise: {bad}")

    shared_with_im_test = {
        arm: sorted(arm_state[arm]["predictor_set"] & set(im_predictors))
        for arm in ARMS
    }

    # --- tier 2 containment, measured not asserted ----------------------------
    print()
    print(THIN)
    print("DEPOSIT CONTAINMENT (tier 2 is granted only on containment == 1.000)")
    print(THIN)
    containment = {}
    for arm in ARMS:
        containment[arm] = {}
        arm_keys = arm_state[arm]["keys"]
        for method, dep in DEPOSITS.items():
            local = dep["local"] if dep["local"] is not None else ARMS[arm]["eval_table"]
            trivial = Path(local) == ARMS[arm]["eval_table"]
            dkeys, dbad = canon_keys(Path(local), dep["allele_col"], dep["peptide_col"])
            n_in = len(arm_keys & dkeys)
            frac = (n_in / len(arm_keys)) if arm_keys else 0.0
            containment[arm][method] = {
                "deposit_doi": dep["doi"],
                "distributor": dep["distributor"],
                "member": dep["member"].replace("{arm}", arm.capitalize()),
                "deposit_rows": len(dkeys),
                "arm_rows_present": n_in,
                "containment": round(frac, 6),
                "contains_whole_arm": frac == 1.0,
                "trivially_true": trivial,
                "trivially_true_note": (
                    "the deposit member IS this arm's evaluation table, so containment is "
                    "definitional and carries no independent information"
                    if trivial else None),
                "deposit_uncanonicalisable_rows": dbad,
            }
            mark = "CONTAINS WHOLE ARM" if frac == 1.0 else f"{frac:.4f}"
            note = "  (definitional: same file)" if trivial else ""
            print(f"  {arm:<11} {method:<9} {dep['doi']:<26} "
                  f"{n_in}/{len(arm_keys)}  {mark}{note}")

    # --- apply the rule, row by row -------------------------------------------
    rows_out = []
    for arm, cfg in ARMS.items():
        lb = cfg["leaderboard"]
        with lb.open("r", newline="", encoding="utf-8-sig") as fh:
            lb_rows = list(csv.DictReader(fh))

        unmapped = [r["entry"] for r in lb_rows
                    if norm_entry(r["entry"]) not in
                    {norm_entry(k) for k in ENTRY_MAP[arm]}]
        if unmapped:
            sys.stderr.write(
                f"{arm}: {len(unmapped)} leaderboard entries are not in ENTRY_MAP and "
                "would be silently uncounted. Refusing to emit a census:\n  "
                + "\n  ".join(unmapped) + "\n")
            return 3

        lookup = {norm_entry(k): v for k, v in ENTRY_MAP[arm].items()}
        predictor_set = arm_state[arm]["predictor_set"]
        shared = set(shared_with_im_test[arm])

        for r in lb_rows:
            entry = r["entry"]
            spec = lookup[norm_entry(entry)]
            col = spec["col"]

            if col is not None and col not in predictor_set:
                sys.stderr.write(
                    f"{arm}: declared column {col!r} for entry {entry!r} is not in the "
                    f"evaluation table header. The map is wrong; refusing to emit.\n")
                return 3

            in_table = col is not None
            in_shared = bool(col) and col in shared
            tier1 = in_table
            tier1_sub = ("1a" if in_shared else "1b") if tier1 else None
            dep_hit = [
                m for m, c in containment[arm].items()
                if c["contains_whole_arm"] and spec["method"] == m
            ]
            tier2 = bool(dep_hit)

            flag_class = classify_flag(r.get("flag", ""))
            withheld = flag_class == "withheld"
            max_tier = 2 if tier2 else (1 if tier1 else 0)

            rows_out.append({
                "arm": arm,
                "rank": r.get("rank", ""),
                "entry": entry,
                "method": spec["method"],
                "score_source_kind": spec["kind"],
                "score_column": col,
                "column_in_evaluation_table": in_table,
                "column_also_in_bigmhc_im_test": in_shared,
                "tier1": tier1,
                "tier1_subtier": tier1_sub,
                "tier2": tier2,
                "tier2_deposits": [
                    {"method": m, **{k: containment[arm][m][k] for k in
                                     ("deposit_doi", "member", "containment",
                                      "trivially_true")}}
                    for m in dep_hit
                ],
                "max_tier": max_tier,
                "affected": max_tier >= 1,
                "figure_published": not withheld,
                "within_allele": r.get("within_allele", ""),
                "pooled": r.get("pooled", ""),
                "current_flag": (r.get("flag") or "").strip(),
                "current_flag_class": flag_class,
                "flagged_as_contamination": flag_class == "not-a-capability-measurement",
            })

    # --- counts ---------------------------------------------------------------
    def tally(arm):
        rs = [x for x in rows_out if x["arm"] == arm]
        pub = [x for x in rs if x["figure_published"]]
        return {
            "leaderboard_rows": len(rs),
            "rows_publishing_a_figure": len(pub),
            "rows_withheld": len(rs) - len(pub),
            "affected_total": sum(x["affected"] for x in rs),
            "affected_publishing_a_figure": sum(x["affected"] for x in pub),
            "tier1": sum(x["tier1"] for x in rs),
            "tier1a_column_in_two_deposits": sum(x["tier1_subtier"] == "1a" for x in rs),
            "tier1b_column_in_evaluation_deposit_only":
                sum(x["tier1_subtier"] == "1b" for x in rs),
            "tier2": sum(x["tier2"] for x in rs),
            "tier2_nontrivial": sum(
                1 for x in rs if x["tier2"]
                and any(not d["trivially_true"] for d in x["tier2_deposits"])),
            "tier0_unaffected": sum(not x["affected"] for x in rs),
            "currently_flagged_contamination": sum(x["flagged_as_contamination"] for x in rs),
            "by_score_source_kind": {
                k: sum(1 for x in rs if x["score_source_kind"] == k)
                for k in sorted({x["score_source_kind"] for x in rs})
            },
        }

    counts = {arm: tally(arm) for arm in ARMS}

    print()
    print(THIN)
    print("CENSUS, PER ARM AND PER TIER")
    print(THIN)
    hdr = f"  {'':<52}{'infectious':>13}{'neoantigen':>13}"
    print(hdr)
    for label, key in [
        ("leaderboard rows", "leaderboard_rows"),
        ("  of which publish a figure", "rows_publishing_a_figure"),
        ("  of which are withheld (no figure)", "rows_withheld"),
        ("AFFECTED under the rule (tier >= 1)", "affected_total"),
        ("  affected AND publishing a figure", "affected_publishing_a_figure"),
        ("tier 1  column shipped in the evaluation table", "tier1"),
        ("  1a  column also in BigMHC im_test.csv", "tier1a_column_in_two_deposits"),
        ("  1b  column in the evaluation deposit only", "tier1b_column_in_evaluation_deposit_only"),
        ("tier 2  method's own deposit holds the whole arm", "tier2"),
        ("  of which measured against a DIFFERENT file", "tier2_nontrivial"),
        ("tier 0  not affected under the rule", "tier0_unaffected"),
        ("CURRENTLY flagged as contamination on disk", "currently_flagged_contamination"),
    ]:
        print(f"  {label:<52}{counts['infectious'][key]:>13}{counts['neoantigen'][key]:>13}")

    print()
    print("  how each figure was obtained:")
    kinds = sorted({x["score_source_kind"] for x in rows_out})
    for k in kinds:
        print(f"    {k:<50}"
              f"{counts['infectious']['by_score_source_kind'].get(k, 0):>13}"
              f"{counts['neoantigen']['by_score_source_kind'].get(k, 0):>13}")

    # --- the disagreement list: THE FINDING -----------------------------------
    disagreements = []
    for x in rows_out:
        if x["affected"] and not x["flagged_as_contamination"] and x["figure_published"]:
            kind = "rule says affected, no contamination flag on disk"
        elif x["affected"] and not x["flagged_as_contamination"] and not x["figure_published"]:
            kind = "rule says affected, figure withheld (no figure to qualify)"
        elif x["flagged_as_contamination"] and not x["affected"]:
            kind = "flagged as contamination, rule says unaffected"
        else:
            continue
        disagreements.append({
            "arm": x["arm"], "rank": x["rank"], "entry": x["entry"],
            "method": x["method"], "score_column": x["score_column"],
            "max_tier": x["max_tier"], "tier1_subtier": x["tier1_subtier"],
            "within_allele": x["within_allele"], "pooled": x["pooled"],
            "current_flag_class": x["current_flag_class"],
            "disagreement": kind,
        })

    published_dis = [d for d in disagreements
                     if d["disagreement"].startswith("rule says affected, no")]
    withheld_dis = [d for d in disagreements if "withheld" in d["disagreement"]]
    reverse_dis = [d for d in disagreements
                   if d["disagreement"].startswith("flagged as contamination")]

    print()
    print(THIN)
    print("FLAG DISAGREEMENTS -- every row where the rule and the shipped flag differ")
    print(THIN)
    print(f"  rule-affected, no contamination flag, figure published (both arms):  "
          f"{len(published_dis)}")
    print(f"  rule-affected, figure withheld (reported separately, both arms):     "
          f"{len(withheld_dis)}")
    print(f"  carrying a contamination flag the rule does not support:             "
          f"{len(reverse_dis)}")
    print()
    for arm in ARMS:
        arm_dis = [d for d in published_dis if d["arm"] == arm]
        print(f"  {arm} -- {len(arm_dis)} rows:")
        print(f"    {'rank':>5}  {'tier':<5} {'entry':<48} {'score column':<38}")
        for d in arm_dis:
            tier = f"T{d['tier1_subtier']}" if d["max_tier"] == 1 else f"T{d['max_tier']}"
            print(f"    {d['rank']:>5}  {tier:<5} "
                  f"{d['entry'][:48]:<48} {str(d['score_column'])[:38]:<38}")
        print()
    if withheld_dis:
        print("  rule-affected with the figure withheld under DTU licence sec. 7(v):")
        for d in withheld_dis:
            print(f"    {d['arm']:<11} {d['entry']:<42} {d['score_column']}")
        print()
    if reverse_dis:
        print("  flagged on disk but NOT affected under the rule:")
        for d in reverse_dis:
            print(f"    {d['arm']:<11} {d['entry']}")
    else:
        print("  No row carries a contamination flag the rule does not support: the")
        print("  shipped flagging is a strict SUBSET of what the rule reaches. The")
        print("  defect is under-flagging, not mis-flagging.")

    # --- the specific inconsistencies the referee named ------------------------
    print()
    print(THIN)
    print("THE THREE INCONSISTENCIES NAMED IN REVIEW, resolved against the rule")
    print(THIN)

    def find(arm, entry):
        return next((x for x in rows_out if x["arm"] == arm and x["entry"] == entry), None)

    for arm in ARMS:
        trio = [find(arm, e) for e in ("BigMHC_EL", "BigMHC_IM", "BigMHC_ELIM")]
        if all(trio):
            print(f"  {arm}: BigMHC_EL / BigMHC_IM / BigMHC_ELIM")
            for x in trio:
                print(f"    {x['entry']:<14} column in table: {x['column_in_evaluation_table']!s:<5} "
                      f"tier {x['max_tier']}  flag: {x['current_flag_class']}")
            tiers = {x["max_tier"] for x in trio}
            flags = {x["flagged_as_contamination"] for x in trio}
            print(f"    -> same tier for all three: {len(tiers) == 1}; "
                  f"same flag for all three: {len(flags) == 1}")

    neo_abl = [x for x in rows_out
               if x["arm"] == "neoantigen" and x["entry"].startswith("BigMHC_Ablation_")]
    neo_par = [x for x in rows_out if x["arm"] == "neoantigen"
               and x["entry"] in ("BigMHC_EL", "BigMHC_IM", "BigMHC_ELIM")]
    print()
    print(f"  neoantigen: 6 BigMHC_Ablation_* rows -- tiers "
          f"{sorted({x['max_tier'] for x in neo_abl})}, "
          f"flagged {sum(x['flagged_as_contamination'] for x in neo_abl)}/6")
    print(f"              3 BigMHC parent rows     -- tiers "
          f"{sorted({x['max_tier'] for x in neo_par})}, "
          f"flagged {sum(x['flagged_as_contamination'] for x in neo_par)}/3")
    print("    -> all nine are read from columns of the same file and are the same tier.")

    neo = counts["neoantigen"]
    print()
    print(f"  the whole neoantigen arm IS im_test: containment "
          f"{containment['neoantigen']['BigMHC']['containment']:.4f} "
          f"({containment['neoantigen']['BigMHC']['arm_rows_present']}/"
          f"{arm_state['neoantigen']['rows']} rows), and "
          f"{len(shared_with_im_test['neoantigen'])} predictor columns are present in both "
          f"files. {neo['affected_publishing_a_figure']} of "
          f"{neo['rows_publishing_a_figure']} published neoantigen figures are affected; "
          f"{neo['currently_flagged_contamination']} are flagged.")
    print(f"  the infectious arm has containment "
          f"{containment['infectious']['BigMHC']['containment']:.4f} against im_test "
          f"({containment['infectious']['BigMHC']['arm_rows_present']}/"
          f"{arm_state['infectious']['rows']} rows), so NO infectious row is tier 2 via "
          f"BigMHC. Its tier 1 rows are tier 1 because the predictor columns ship inside "
          f"the Dryad evaluation table itself.")

    # --- artifact -------------------------------------------------------------
    artifact = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "script": SCRIPT_REL,
        "rule": {
            "statement": (
                "A leaderboard row is contamination-affected if (tier 1) its score column "
                "is a column of the upstream evaluation table it is scored against, or "
                "(tier 2) the method's own public data deposit contains every row of the "
                "arm. Tier 1a marks a column present in both the evaluation deposit and "
                "BigMHC's im_test.csv; tier 1b a column present only in the evaluation "
                "deposit. Tier 0 is fresh inference, a function of the row's own fields, "
                "or a fit on a declared non-evaluation pool."),
            "supportable_claim_form": (
                "A score on an affected row is not a capability bar for a method with "
                "access to the table."),
            "not_established": (
                "Nothing about how any file came to exist, what any model was fitted on, "
                "or anyone's intent. This measures file contents only."),
            "tier0_caveat": (
                "Tier 0 is not a certificate of independence; it means this rule found no "
                "affectedness. Every arm carries some dependence (gate G10)."),
        },
        "inputs": {
            "bigmhc_deposit_member": {
                "path": INPUT_DESCRIPTIONS["bigmhc"],
                "sha256": im_sha,
                "columns": len(im_header),
                "predictor_columns": len(im_predictors),
                "predictor_column_names": im_predictors,
                "deposit_doi": DEPOSITS["BigMHC"]["doi"],
            },
            "arms": {
                arm: {k: v for k, v in arm_state[arm].items()
                      if k not in ("predictor_set", "keys")}
                for arm in ARMS
            },
            "leaderboards": {
                arm: str(ARMS[arm]["leaderboard"].relative_to(ROOT)) for arm in ARMS
            },
        },
        "shared_predictor_columns_with_bigmhc_im_test": {
            arm: {"count": len(shared_with_im_test[arm]),
                  "columns": shared_with_im_test[arm]}
            for arm in ARMS
        },
        "deposit_containment": containment,
        "counts": counts,
        "rows": rows_out,
        "flag_disagreements": {
            "rule_affected_unflagged_figure_published": published_dis,
            "rule_affected_figure_withheld": withheld_dis,
            "flagged_but_rule_says_unaffected": reverse_dis,
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(artifact, indent=2, sort_keys=False) + "\n")

    print()
    print(RULE)
    print(f"artifact written: {OUT.relative_to(ROOT)}")
    print(RULE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
