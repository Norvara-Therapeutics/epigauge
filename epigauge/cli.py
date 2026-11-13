#!/usr/bin/env python3
"""The `epigauge` command: list, fetch, score, leaderboard, validate, fixtures.

    epigauge list
    epigauge fetch infectious
    epigauge score PREDICTIONS.csv --set infectious
    epigauge leaderboard --set infectious
    epigauge validate PREDICTIONS.csv
    epigauge fixtures

Exit codes: 0 every check passed, 1 a check FAILED, 2 could NOT RUN. 2 is never
conflated with 0, because the failures here are absences that read as answers -- a
missed allele join returning exactly 0.5000, `0 of 2,093` meaning "nothing looked".
`bootstrap`, `baselines` and `data` are imported defensively; a command needing a
missing one exits 2 naming it, never a figure computed without it. Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from . import report
from .metrics import DegenerateGroupingError, auroc, stratified_auroc
from .report import FAIL, PASS, UNRUNNABLE, MethodRow, Table

try:  # intervals, paired deltas, minimum detectable difference
    from . import bootstrap
except ImportError:  # pragma: no cover - exercised only before bootstrap.py lands
    bootstrap = None

try:  # the baseline roster a submission has to beat
    from . import baselines
except ImportError:  # pragma: no cover
    baselines = None

try:  # row set loading and fetching
    from . import data
except ImportError:  # pragma: no cover
    data = None

__all__ = [
    "EXIT_OK",
    "EXIT_CHECK_FAILED",
    "EXIT_COULD_NOT_RUN",
    "CouldNotRun",
    "SubmissionError",
    "Submission",
    "read_submission",
    "validate_submission",
    "exit_code_for",
    "build_parser",
    "main",
]

EXIT_OK = 0
EXIT_CHECK_FAILED = 1
EXIT_COULD_NOT_RUN = 2


class CouldNotRun(RuntimeError):
    """The command could not be attempted. Always exit 2, never exit 0."""


class SubmissionError(ValueError):
    """The submission file is present but not readable as a submission. Exit 1."""


# ---------------------------------------------------------------------------
# submissions
#
# `row_id,score` plus a DECLARED score direction, mandatory and with no default:
# inferring it from the data means inferring it from the labels.

_ID_COLUMNS = ("row_id", "record_id", "rowid", "id")
_SCORE_COLUMNS = ("score", "prediction", "pred", "value")
_DIRECTION_KEYS = ("direction", "score_direction", "sign")

HIGHER, LOWER = "higher-is-better", "lower-is-better"
_DIRECTIONS = {
    HIGHER: HIGHER, "higher": HIGHER, "high": HIGHER, "descending": HIGHER,
    "max": HIGHER, "maximise": HIGHER, "maximize": HIGHER, "score": HIGHER,
    LOWER: LOWER, "lower": LOWER, "low": LOWER, "ascending": LOWER,
    "min": LOWER, "minimise": LOWER, "minimize": LOWER, "rank": LOWER,
    "percentile-rank": LOWER,
}


def _normalise_direction(raw: str | None) -> str | None:
    if raw is None:
        return None
    key = raw.strip().lower().replace("_", "-").replace(" ", "-")
    return _DIRECTIONS.get(key)


@dataclass
class Submission:
    path: str
    ids: list[str] = field(default_factory=list)
    scores: dict[str, float] = field(default_factory=dict)
    direction: str | None = None
    direction_raw: str | None = None
    unparseable: list[tuple[str, str]] = field(default_factory=list)
    nonfinite: list[tuple[str, float]] = field(default_factory=list)
    label: str = "submission"

    @property
    def duplicates(self) -> list[str]:
        seen, dupes = set(), []
        for i in self.ids:
            if i in seen and i not in dupes:
                dupes.append(i)
            seen.add(i)
        return dupes


def read_submission(path, *, direction: str | None = None,
                    label: str | None = None) -> Submission:
    """Parse a submission. Tolerant by design: the CHECKS live in `validate_submission`.

    A missing id or score column raises SubmissionError; anything merely wrong is
    recorded on the Submission so `validate` prints every problem in one pass.
    """
    source = Path(path)
    if not source.is_file():
        raise CouldNotRun(
            f"no such submission file: {source}. Nothing was read, so nothing passed."
        )
    text = source.read_text(encoding="utf-8", errors="replace")

    meta: dict[str, str] = {}
    body: list[str] = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            key, _, value = line.lstrip().lstrip("#").strip().partition(":")
            if value:
                meta[key.strip().lower().replace(" ", "_")] = value.strip()
        else:
            body.append(line)

    reader = csv.DictReader(body)
    fields = [f.strip() for f in (reader.fieldnames or []) if f]
    lower = {f.lower(): f for f in fields}
    id_col = next((lower[c] for c in _ID_COLUMNS if c in lower), None)
    score_col = next((lower[c] for c in _SCORE_COLUMNS if c in lower), None)
    if id_col is None or score_col is None:
        raise SubmissionError(
            f"{source}: need a `row_id` column and a `score` column; found {fields}. "
            "The join key is row_id and nothing else -- a (peptide, allele) join across "
            "two spellings of the same allele returned 0 of 833 matches here and hid a "
            "whole contaminated benchmark."
        )

    raw_direction = next((meta[k] for k in _DIRECTION_KEYS if k in meta), None)
    sub = Submission(path=str(source), label=label or source.stem)

    dir_col = lower.get("direction")
    for record in reader:
        row_id = (record.get(id_col) or "").strip()
        if not row_id:
            continue
        sub.ids.append(row_id)
        if dir_col and raw_direction is None:
            raw_direction = (record.get(dir_col) or "").strip() or None
        raw = (record.get(score_col) or "").strip()
        try:
            value = float(raw)
        except (TypeError, ValueError):
            sub.unparseable.append((row_id, raw))
            continue
        if not math.isfinite(value):
            sub.nonfinite.append((row_id, value))
            continue
        sub.scores.setdefault(row_id, value)

    sub.direction_raw = direction or raw_direction
    sub.direction = _normalise_direction(sub.direction_raw)
    return sub


def validate_submission(sub: Submission, expected_ids=None) -> list[tuple[str, str, str]]:
    """Format check only. Returns `[(check, status, detail), ...]`; computes no AUROC.

    `expected_ids=None` means the row set could not be loaded, so completeness reports
    `could-not-run`, NOT a pass: a coverage check that did not look is `0 of 2,093`.
    """
    checks: list[tuple[str, str, str]] = []

    dupes = sub.duplicates
    if expected_ids is None:
        checks.append((
            "every benchmark row_id present exactly once", UNRUNNABLE,
            "the benchmark row set could not be loaded, so coverage was never compared. "
            "This is not a pass.",
        ))
    else:
        expected = list(expected_ids)
        expected_set = set(expected)
        got = set(sub.ids)
        missing = [i for i in expected if i not in got]
        extra = sorted(got - expected_set)
        if missing or extra or dupes:
            bits = []
            if missing:
                bits.append(f"{len(missing)} of {len(expected)} benchmark rows missing "
                            f"(e.g. {', '.join(missing[:5])})")
            if extra:
                bits.append(f"{len(extra)} row_id(s) not in the benchmark "
                            f"(e.g. {', '.join(extra[:5])})")
            if dupes:
                bits.append(f"{len(dupes)} duplicated row_id(s) "
                            f"(e.g. {', '.join(dupes[:5])})")
            checks.append(("every benchmark row_id present exactly once", FAIL,
                           "; ".join(bits) + ". A partial submission is scored on a "
                           "different population, which is not the benchmark."))
        else:
            checks.append(("every benchmark row_id present exactly once", PASS,
                           f"{len(expected)} rows, each exactly once"))

    if sub.unparseable or sub.nonfinite:
        bits = []
        if sub.unparseable:
            bits.append(f"{len(sub.unparseable)} unparseable "
                        f"(e.g. {sub.unparseable[0][0]}={sub.unparseable[0][1]!r})")
        if sub.nonfinite:
            bits.append(f"{len(sub.nonfinite)} non-finite "
                        f"(e.g. {sub.nonfinite[0][0]}={sub.nonfinite[0][1]})")
        checks.append(("every score is a finite number", FAIL,
                       "; ".join(bits) + ". A NaN is ordered arbitrarily by the "
                       "mid-rank, so it becomes a silent per-row coin flip."))
    elif not sub.scores:
        checks.append(("every score is a finite number", FAIL,
                       "no scores were read at all"))
    else:
        checks.append(("every score is a finite number", PASS,
                       f"{len(sub.scores)} scores"))

    if sub.direction is None:
        detail = (f"declared {sub.direction_raw!r}, which is not one of "
                  f"{HIGHER} / {LOWER}" if sub.direction_raw else
                  "no direction declared. Add `# direction: higher-is-better` as the "
                  "first line, a `direction` column, or pass --direction.")
        checks.append(("a score direction is declared", FAIL,
                       detail + " The direction is not inferable from the scores without "
                       "reading the labels; a family of percentile-rank columns read "
                       "with the wrong convention produced a 0.4650 'anti-predictive "
                       "floor' that vanished when flipped."))
    else:
        checks.append(("a score direction is declared", PASS, sub.direction))

    distinct = len(set(sub.scores.values()))
    if distinct <= 1:
        checks.append(("the submission is not all-constant", FAIL,
                       f"{distinct} distinct score value(s) across "
                       f"{len(sub.scores)} rows. A constant score is exactly 0.5 under "
                       "the mid-rank correction, which reads as a measured null rather "
                       "than as an empty column."))
    else:
        checks.append(("the submission is not all-constant", PASS,
                       f"{distinct} distinct values"))
    return checks


def exit_code_for(checks) -> int:
    """FAIL outranks could-not-run, and could-not-run never reads as a pass."""
    if any(status == FAIL for _, status, _ in checks):
        return EXIT_CHECK_FAILED
    if any(status == UNRUNNABLE for _, status, _ in checks):
        return EXIT_COULD_NOT_RUN
    return EXIT_OK


# ---------------------------------------------------------------------------
# adapters over the sibling modules
#
# Every call into `data`, `baselines` and `bootstrap` goes through an accessor below,
# each raising CouldNotRun naming the module and the attempt, and never substituting a
# figure. A licence refusal is also a could-not-run: `score --set tesla608` exits 2.

DEFAULT_REPLICATES = 1000


def _require(module, name: str, why: str):
    if module is None:
        raise CouldNotRun(
            f"epigauge.{name} is not available, so {why} could not be attempted. "
            f"Exiting {EXIT_COULD_NOT_RUN} rather than reporting a figure computed "
            "without it."
        )
    return module


def _field(item, *names, default=None):
    for name in names:
        if hasattr(item, name):
            return getattr(item, name)
        try:
            if name in item:
                return item[name]
        except TypeError:
            pass
    return default


def _data_call(what: str, fn, *args, **kwargs):
    """Every entry into `data`, with its refusals turned into exit 2, reason kept."""
    try:
        return fn(*args, **kwargs)
    except getattr(data, "DataError", ()) as error:
        raise CouldNotRun(f"{what}: {type(error).__name__}: {error}") from None
    except (OSError, ValueError) as error:
        raise CouldNotRun(f"{what}: {type(error).__name__}: {error}") from None


def _resolve_set(name: str) -> str:
    """`--set infectious` -> `titanian_infectious`, or a refusal naming the candidates.

    A short name matching two sets is an ERROR, not a pick: a key that resolves to the
    wrong thing gives a clean-looking number on the wrong rows.
    """
    _require(data, "data", f"resolving the row set {name!r}")
    known = [s.set_id for s in _data_call("listing row sets", data.available)]
    if name in known:
        return name
    hits = [s for s in known if s.endswith(f"_{name}") or name in s]
    if len(hits) == 1:
        return hits[0]
    raise CouldNotRun(
        f"{'ambiguous' if hits else 'unknown'} row set {name!r}"
        + (f"; it matches {', '.join(hits)}" if hits else "")
        + f". Known sets: {', '.join(sorted(known))}. `epigauge list` prints them with "
        "their licence and whether they ship."
    )


def _row_sets() -> list[dict]:
    _require(data, "data", "listing the row sets")
    specs = []
    for info in _data_call("listing row sets", data.available):
        entry = _data_call(f"describing {info.set_id}", data.describe, info.set_id)
        specs.append({
            "name": info.set_id,
            "n_rows": info.rows,
            "n_pos": entry.get("positives"),
            "licence": info.licence,
            "status": info.status,
            "note": info.title,
        })
    return specs


def _load_rows(set_id: str, *, scores: bool = False) -> list[dict]:
    _require(data, "data", f"loading the {set_id!r} row set")
    rows = _data_call(f"loading {set_id}", data.load, set_id, scores=scores)
    if not rows:
        raise CouldNotRun(
            f"the {set_id!r} row set loaded zero rows. An empty set scores nothing, and "
            "a printed 0 that can only ever be 0 is the report that read as 'nothing "
            "found'."
        )
    return list(rows)


def _row_ids(rows) -> list[str]:
    return [str(_field(r, "record_id", "row_id", "id")) for r in rows]


def _labels(rows) -> list[int]:
    values = [_field(r, "label", "y") for r in rows]
    if any(v is None for v in values):
        raise CouldNotRun(
            f"{sum(v is None for v in values)} of {len(values)} rows carry no label. An "
            "absent field read as a category is how all 2,093 rows came to be declared "
            "self-derived."
        )
    return [int(v) for v in values]


def _groups(rows) -> list[str]:
    values = [_field(r, "allele", "mhc_allele", "mhc", "stratum") for r in rows]
    if any(v in (None, "") for v in values):
        raise CouldNotRun(
            f"{sum(1 for v in values if v in (None, ''))} of {len(values)} rows carry no "
            "allele. A uniformly-missing grouping turns the within-allele estimate back "
            "into pooled silently: 0.584138 where the correct figure is 0.492317."
        )
    return [str(v) for v in values]


def _peptides(rows) -> list[str]:
    return [str(_field(r, "peptide")) for r in rows]


def _baseline_rows(rows) -> list[dict]:
    """`baselines` spells the allele `mhc_allele`, `data` rows spell it `allele`.

    Aliased here, where the two meet. Dropping rows that fail to join would be the
    0-of-833 bug again.
    """
    out = []
    for row in rows:
        if _field(row, "mhc_allele") is not None:
            out.append(row)
        else:
            copy = dict(row)
            copy["mhc_allele"] = _field(row, "allele", "mhc")
            out.append(copy)
    return out


def _baseline_columns(rows, fit_rows=None) -> tuple[dict[str, list[float]], dict[str, str]]:
    """The baseline roster. Returns `(columns, refusals)`; a refusal is never a gap.

    `eval_ids=()` is typed by hand only for baselines declaring `uses_labels = False`;
    the loop refuses it for any other, since a baseline fitted on the evaluation labels
    is indistinguishable from a good one. The label-reading baselines (allele prior,
    PWM) need a separate `--fit-set` pool, and then run under the live leak guard.
    """
    _require(baselines, "baselines", "adding the baseline columns")
    evaluation = _baseline_rows(rows)
    columns: dict[str, list[float]] = {}
    refusals: dict[str, str] = {}

    for cls in (baselines.Hydropathy, baselines.Length):
        if getattr(cls, "uses_labels", True):
            refusals[cls.name] = ("declares uses_labels=True; it may not be fitted on "
                                  "the evaluation rows")
            continue
        try:
            fitted = cls(eval_ids=()).fit(evaluation)
            columns[fitted.name] = [float(v) for v in fitted.score(evaluation)]
        except baselines.BaselineError as error:
            refusals[cls.name] = str(error)

    label_reading = (baselines.AlleleOnly, baselines.PWM)
    if fit_rows is None:
        for cls in label_reading:
            refusals[cls.name] = (
                "reads labels, so it needs a fitting pool disjoint from the evaluation "
                "set. Pass --fit-set to include it; it is ABSENT here, not zero"
            )
    else:
        eval_ids = set(_row_ids(evaluation))
        pool = _baseline_rows(fit_rows)
        for cls in label_reading:
            try:
                fitted = cls(eval_ids=eval_ids).fit(pool)
                columns[fitted.name] = [float(v) for v in fitted.score(evaluation)]
            except baselines.BaselineError as error:
                refusals[cls.name] = f"{type(error).__name__}: {error}"
    return columns, refusals


def _reference_columns(set_id: str, rows) -> tuple[dict[str, list[float]], list[str]]:
    """The shipped competitor columns, plus the withheld ones named as withheld.

    "Missing" and "withheld" are different findings, and only one is about the method,
    so a licence-withheld column goes in the footnotes rather than being dropped.
    """
    _require(data, "data", "reading the shipped reference columns")
    scored = _load_rows(set_id, scores=True)
    raw: dict[str, list] = {}
    for row in scored:
        for column, value in (_field(row, "scores", default={}) or {}).items():
            raw.setdefault(column, []).append(value)

    columns: dict[str, list[float]] = {}
    for column, values in raw.items():
        if len(values) != len(scored):
            continue
        try:
            floats = [float(v) for v in values]
        except (TypeError, ValueError):
            continue
        if len(set(floats)) <= 1:
            continue  # a constant column is exactly 0.5 by construction, not a result
        columns[column] = floats

    notes: list[str] = []
    info = _data_call(f"reading {set_id} score columns", data.score_columns, set_id)
    withheld = sorted(info.get("withheld") or {})
    if withheld:
        notes.append(f"{len(withheld)} predictor column(s) are WITHHELD, not missing: "
                     f"{', '.join(withheld[:6])}"
                     f"{', ...' if len(withheld) > 6 else ''}.")
        notes.append("NetMHCpan columns in particular are absent because DTU 7(v) "
                     "forbids publishing results")
        notes.append("run on the product without written consent, which does not exist. "
                     "See data.score_columns().")
    return columns, notes


def _mde(set_id, labels, groups, peptides, arms, args) -> tuple[float, str]:
    """The minimum detectable paired difference: measured if it can be, assumed if not.

    Clusters need source proteomes, so it is measured only with `--proteome`. Otherwise
    the documented 0.059-0.095 band prints "assumed, NOT computed" in the header.
    """
    basis = report.DOCUMENTED_MDE_BASIS
    if "--proteome" not in basis:
        basis += "; pass --proteome FASTA to measure it for YOUR comparison"
    fallback = (report.DOCUMENTED_MDE, basis)
    proteomes = list(getattr(args, "proteome", None) or [])
    if bootstrap is None or not proteomes or len(arms) < 2:
        if proteomes and bootstrap is None:
            print("note: --proteome given but epigauge.bootstrap is not importable; "
                  "the MDE below is the documented band, not a measurement.",
                  file=sys.stderr)
        return fallback

    (a_name, a_scores), (b_name, b_scores) = arms[0], arms[1]
    try:
        clusters = bootstrap.assign_clusters(peptides, proteomes)
        recipe = bootstrap.SeedRecipe(version="epigauge-cli-1", set_id=set_id,
                                      comparison_id=f"{a_name} vs {b_name}")
        matrix = bootstrap.draws(clusters.n_clusters, args.replicates,
                                 seed_recipe=recipe)
        delta = bootstrap.paired(stratified_auroc, labels, a_scores, b_scores, groups,
                                 clusters, matrix)
        measured = bootstrap.mde(delta)
    except Exception as error:  # noqa: BLE001 - any refusal degrades to the band
        print(f"note: the MDE could not be measured ({type(error).__name__}: {error}); "
              "the table falls back to the documented band.", file=sys.stderr)
        return fallback
    return float(measured.mde), (
        f"MEASURED on {a_name} vs {b_name}: se {measured.se_delta:.4f}, n_eff "
        f"{measured.n_eff:.0f} after design effect {measured.design_effect:.3f}, "
        f"basis {measured.basis}"
    )


# ---------------------------------------------------------------------------
# scoring


def _method_row(name: str, labels, scores, groups, *, sign_flipped=False,
                direction_declared=False, flag=None) -> MethodRow:
    """One table row. Stratified and pooled are computed together or not at all."""
    common = {"name": name, "sign_flipped": sign_flipped, "flag": flag,
              "direction_declared": direction_declared}
    try:
        strat = stratified_auroc(labels, scores, groups)
    except DegenerateGroupingError as error:
        return MethodRow(note=f"degenerate grouping: {error}", **common)
    return MethodRow(stratified=strat, pooled=auroc(labels, scores), **common)


def _oriented(name: str, direction: str, values):
    """Apply the declared direction, and make any flip visible in the printed name."""
    if direction == LOWER:
        return f"{name} [declared lower-is-better, SIGN-FLIPPED]", [-v for v in values]
    return name, list(values)


def _join(sub: Submission, rows) -> list[float]:
    ids = _row_ids(rows)
    missing = [i for i in ids if i not in sub.scores]
    if missing:
        raise SubmissionError(
            f"{len(missing)} of {len(ids)} benchmark rows have no score "
            f"(e.g. {', '.join(missing[:5])}). Run `epigauge validate` first; a partial "
            "submission would be scored on a different population than the benchmark."
        )
    return [sub.scores[i] for i in ids]


def _rank_arms(labels, groups, arms):
    """Scorable arms, strongest first. A degenerate arm is dropped, never ranked at 0.5."""
    scored = []
    for name, values in arms:
        try:
            scored.append((stratified_auroc(labels, values, groups).value, name, values))
        except DegenerateGroupingError:
            continue
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(name, values) for _value, name, values in scored]


def _set_header(set_id: str, rows, labels) -> str:
    n_pos = sum(1 for y in labels if y == 1)
    return (f"{set_id}  n={len(rows)}  {n_pos} positive / {len(rows) - n_pos} negative "
            f"({n_pos / len(rows):.1%})  {len(set(_groups(rows)))} strata")


def _baseline_rows_for_table(rows, fit_rows, table_rows, arms) -> None:
    """Append the baseline roster to a table, refusals included as named absences."""
    try:
        columns, refusals = _baseline_columns(rows, fit_rows)
    except CouldNotRun as error:
        table_rows.append(MethodRow(name="BASELINE roster", note=str(error)))
        return
    labels, groups = _labels(rows), _groups(rows)
    for name, values in sorted(columns.items()):
        table_rows.append(_method_row(f"BASELINE {name}", labels, values, groups,
                                      flag=report.capability_flag(name)))
        arms.append((f"BASELINE {name}", values))
    for name, reason in sorted(refusals.items()):
        table_rows.append(MethodRow(name=f"BASELINE {name}", note=reason,
                                    flag=report.capability_flag(name)))


# ---------------------------------------------------------------------------
# commands


def cmd_list(args) -> int:
    sys.stdout.write(report.render_row_sets(_row_sets()))
    return EXIT_OK


def cmd_fetch(args) -> int:
    _require(data, "data", "fetching a row set")
    if not args.set:
        sys.stdout.write(report.render_row_sets(_row_sets()))
        print("name one of the sets above: `epigauge fetch titanian_ada`.")
        return EXIT_COULD_NOT_RUN
    set_id = _resolve_set(args.set)
    path = _data_call(f"fetching {set_id}", data.fetch, set_id)
    print(f"{set_id}: {path}")
    return EXIT_OK


def cmd_validate(args) -> int:
    sub = read_submission(args.predictions, direction=args.direction)
    try:
        expected = _row_ids(_load_rows(_resolve_set(args.set)))
    except CouldNotRun as error:
        print(f"note: {error}", file=sys.stderr)
        expected = None
    checks = validate_submission(sub, expected)
    sys.stdout.write(report.render_validation(sub.path, checks))
    print("this is a format check only; no AUROC was computed.")
    return exit_code_for(checks)


def cmd_score(args) -> int:
    set_id = _resolve_set(args.set)
    rows = _load_rows(set_id)
    sub = read_submission(args.predictions, direction=args.direction)
    checks = validate_submission(sub, _row_ids(rows))
    code = exit_code_for(checks)
    if code != EXIT_OK:
        sys.stdout.write(report.render_validation(sub.path, checks))
        print("refusing to score a submission that failed the format check.")
        return code

    labels, groups = _labels(rows), _groups(rows)
    name, scores = _oriented(sub.label, sub.direction, _join(sub, rows))

    arms = [(name, scores)]
    table_rows = [_method_row(name, labels, scores, groups,
                              sign_flipped=sub.direction == LOWER,
                              direction_declared=True)]
    fit_rows = _load_rows(_resolve_set(args.fit_set)) if args.fit_set else None
    _baseline_rows_for_table(rows, fit_rows, table_rows, arms)

    # The submission stays arm A: the MDE wanted is the one for THEIR comparison,
    # against the strongest arm rather than whichever baseline was built first.
    mde, basis = _mde(set_id, labels, groups, _peptides(rows),
                      arms[:1] + _rank_arms(labels, groups, arms[1:]), args)
    sys.stdout.write(report.render_table(Table(
        title=f"SCORECARD  {sub.path}",
        subtitle=_set_header(set_id, rows, labels),
        rows=table_rows,
        mde=mde,
        mde_basis=basis,
        footnotes=[
            "baselines are rows in this table, not footnotes: an unfitted hydropathy "
            "scale places",
            "above most published methods here, and an allele lookup with no peptide "
            "input reaches",
            "pooled 0.7122 at within-allele exactly 0.5000.",
            "no interval is printed. Intervals require the cluster bootstrap and a "
            "declared basis;",
            "row-level intervals on clustered rows understate the width by the design "
            "effect.",
        ],
    )))
    return EXIT_OK


def cmd_leaderboard(args) -> int:
    set_id = _resolve_set(args.set)
    rows = _load_rows(set_id)
    labels, groups = _labels(rows), _groups(rows)

    columns, notes = _reference_columns(set_id, rows)
    table_rows: list[MethodRow] = []
    arms: list[tuple[str, list[float]]] = []
    for column, values in sorted(columns.items()):
        flipped = report.looks_like_percentile_rank(column)
        name, oriented = (report.flip_percentile_rank(column, values) if flipped
                          else (column, list(values)))
        table_rows.append(_method_row(name, labels, oriented, groups,
                                      sign_flipped=flipped,
                                      flag=report.capability_flag(column)))
        arms.append((name, oriented))

    fit_rows = _load_rows(_resolve_set(args.fit_set)) if args.fit_set else None
    _baseline_rows_for_table(rows, fit_rows, table_rows, arms)

    if not any(r.stratified is not None for r in table_rows):
        raise CouldNotRun(f"no scorable column on the {set_id!r} row set")

    mde, basis = _mde(set_id, labels, groups, _peptides(rows),
                      _rank_arms(labels, groups, arms), args)
    sys.stdout.write(report.render_table(Table(
        title=f"LEADERBOARD  {set_id}",
        subtitle=_set_header(set_id, rows, labels),
        rows=table_rows,
        mde=mde,
        mde_basis=basis,
        footnotes=notes + [
            "a flagged row is not a weak method: it is a column measuring something "
            "other than",
            "capability, and the flag says which. Read the reason, not the position.",
        ],
    )))
    return EXIT_OK


def cmd_fixtures(args) -> int:
    """Run every fixture file through its standalone runner, so pytest is not required."""
    tests = sorted((Path(__file__).resolve().parent / "tests").glob("test_*.py"))
    if not tests:
        raise CouldNotRun("no fixture files found under epigauge/tests/")
    failed = []
    for path in tests:
        print(f"\n=== {path.name} " + "=" * max(0, 60 - len(path.name)))
        completed = subprocess.run([sys.executable, "-I", str(path)], check=False)
        if completed.returncode != 0:
            failed.append(path.name)
    print(f"\n{len(tests)} fixture file(s), {len(failed)} failing: "
          f"{', '.join(failed) if failed else 'none'}")
    if failed:
        print("a check that cannot be demonstrated to fail is a failed check; these "
              "demonstrate it.")
    return EXIT_CHECK_FAILED if failed else EXIT_OK


# ---------------------------------------------------------------------------
# argument parsing


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="epigauge",
        description="A measurement instrument for class-I T-cell epitope immunogenicity "
                    "prediction. Reports the pair-weighted stratified estimator with "
                    "pooled beside it, never one alone, and never ranks on pooled.",
        epilog="exit codes: 0 clean, 1 a check failed, 2 could not run. 2 is never "
               "conflated with 0.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subs = parser.add_subparsers(dest="command", required=True)

    def with_set(sub):
        sub.add_argument("--set", default="titanian_infectious", metavar="NAME",
                         help="row set, as printed by `epigauge list`. A unique short "
                              "name resolves (`infectious` -> `titanian_infectious`); "
                              "an ambiguous one is refused rather than picked")

    def with_mde(sub):
        sub.add_argument("--fit-set", default=None, metavar="NAME",
                         help="row set to fit the label-reading baselines on. It must be "
                              "disjoint from --set by record_id, and the leak guard "
                              "raises if it is not. Without it the allele prior and the "
                              "PWM are reported ABSENT with the reason, not omitted")
        sub.add_argument("--proteome", action="append", metavar="FASTA",
                         help="source proteome FASTA, repeatable. Supplying them lets "
                              "the MDE be MEASURED from cluster draws instead of taken "
                              "from the documented 0.059-0.095 band")
        sub.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES,
                         metavar="B",
                         help=f"bootstrap replicates behind a measured MDE (default "
                              f"{DEFAULT_REPLICATES}); the honest p-value floor is "
                              "1/(B+1)")

    listing = subs.add_parser(
        "list", help="available row sets",
        description="Row sets with size, class balance, licence, and whether they ship "
                    "or must be fetched. `fetch required` is a licence consequence: "
                    "TESLA, VDJdb, McPAS-TCR and the NetMHCpan score columns may not be "
                    "redistributed from here.")
    listing.set_defaults(func=cmd_list)

    fetch = subs.add_parser(
        "fetch", help="download a row set that does not ship",
        description="Fetch and build a row set locally. Every upstream file is pinned by "
                    "URL, SHA-256 and byte count, so a changed upstream fails loudly "
                    "rather than silently changing a published number.")
    fetch.add_argument("set", nargs="?", default=None, metavar="SET",
                       help="row set name; omit to list what is fetchable")
    fetch.set_defaults(func=cmd_fetch)

    score = subs.add_parser(
        "score", help="score a submission against a row set",
        description="Join a submission on row_id, compute the pair-weighted stratified "
                    "estimator with pooled printed beside it, add the baseline rows, and "
                    "print a scorecard. The submission is format-checked first and is "
                    "not scored if that check fails.")
    score.add_argument("predictions", metavar="PREDICTIONS.csv",
                       help="CSV with row_id and score columns, plus a declared "
                            "direction (`# direction: higher-is-better`)")
    with_set(score)
    with_mde(score)
    score.add_argument("--direction", default=None, choices=sorted({HIGHER, LOWER}),
                       help="declare the score direction here instead of in the file")
    score.set_defaults(func=cmd_score)

    board = subs.add_parser(
        "leaderboard", help="the full ranked table for a row set",
        description="Every shipped reference column plus the baselines, ranked on the "
                    "stratified estimator. Percentile-rank columns are sign-flipped and "
                    "say so in the printed name; columns that are not capability "
                    "measurements carry the specific reason.")
    with_set(board)
    with_mde(board)
    board.set_defaults(func=cmd_leaderboard)

    validate = subs.add_parser(
        "validate", help="format check only, no scoring",
        description="Check that every benchmark row_id is present exactly once, that "
                    "every score is finite, that a score direction is declared, and that "
                    "the submission is not all-constant. Computes no AUROC. Exits 2, not "
                    "0, if the benchmark row set could not be loaded to compare against.")
    validate.add_argument("predictions", metavar="PREDICTIONS.csv")
    with_set(validate)
    validate.add_argument("--direction", default=None, choices=sorted({HIGHER, LOWER}))
    validate.set_defaults(func=cmd_validate)

    fixtures = subs.add_parser(
        "fixtures", help="run every regression fixture",
        description="Run each tests/test_*.py through its standalone runner, so a "
                    "reviewer without pytest can reproduce every documented failure. "
                    "Exits non-zero if any fixture fails.")
    fixtures.set_defaults(func=cmd_fixtures)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except CouldNotRun as error:
        print(f"epigauge {args.command}: COULD NOT RUN: {error}", file=sys.stderr)
        return EXIT_COULD_NOT_RUN
    except SubmissionError as error:
        print(f"epigauge {args.command}: CHECK FAILED: {error}", file=sys.stderr)
        return EXIT_CHECK_FAILED
    except report.ReportError as error:
        print(f"epigauge {args.command}: REFUSED TO PRINT: {error}", file=sys.stderr)
        return EXIT_CHECK_FAILED


if __name__ == "__main__":
    raise SystemExit(main())
