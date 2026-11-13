"""Formatting. Separated from `cli.py` so that every reporting rule has a fixture.

The rules live here as code and `render_table` raises rather than emit a table that
breaks one: pooled and stratified ship together, ranking reads `stratified` only, every
figure carries its stratum and evaluable-row counts, the MDE heads the table, a
non-capability column carries a specific reason, and a sign flip shows in the name.
Standard library only.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "ReportError",
    "MethodRow",
    "Table",
    "DOCUMENTED_MDE",
    "DOCUMENTED_MDE_BASIS",
    "NOT_A_CAPABILITY",
    "SIGN_FLIP_TAG",
    "capability_flag",
    "looks_like_percentile_rank",
    "flip_percentile_rank",
    "render_table",
    "render_row_sets",
    "render_validation",
]


class ReportError(ValueError):
    """A table violates a reporting rule. Raised instead of printing the table."""


#: Fallback minimum detectable paired difference, used only when the cluster bootstrap
#: cannot run. The measured band is 0.059-0.095, roughly double the unmeasured 0.03-0.05
#: band it supersedes; the fallback is the conservative end, because an MDE half the real
#: one lets a margin be called resolvable when it is not.
DOCUMENTED_MDE = 0.095
DOCUMENTED_MDE_BASIS = (
    "assumed, NOT computed: the conservative end of the MEASURED 0.059-0.095 band "
    "(clustered paired SE 0.021-0.034 at z-sum 2.801585, design effect 1.10-1.72, "
    "n=2,093 infectious). Supersedes the unmeasured 0.03-0.05 band that circulated "
    "earlier. A measured MDE for YOUR comparison needs the cluster bootstrap: "
    "pass --proteome"
)

#: Shipped reference columns that measure something other than capability, each with
#: the reason specific to it.
NOT_A_CAPABILITY: dict[str, str] = {
    "titanian": "label recoverable from the IEDB assay table at 0.9101, 100% coverage",
    "titanian_no_pretrain": "same assay-table provenance as the other TITANiAN columns",
    "titanian_ablation_": "same assay-table provenance as the other TITANiAN columns",
    "TITANiAN-IM": "label recoverable from the IEDB assay table at 0.9101, 100% coverage",
    "TITANiAN-IM-baseline": "same assay-table provenance as TITANiAN-IM",
    "bigmhc_im": "im_test.csv contains this benchmark's neoantigen arm entirely "
                 "(833/833 pairs, 833/833 labels, all 25 shared columns byte-identical)",
    "bigmhc_elim": "im_test.csv contains this benchmark's neoantigen arm entirely",
    "BigMHC_IM": "im_test.csv contains this benchmark's neoantigen arm entirely "
                 "(833/833 pairs, 833/833 labels, all 25 shared columns byte-identical)",
    "BigMHC_ELIM": "im_test.csv contains this benchmark's neoantigen arm entirely",
    "AlleleOnly": "an allele positive-rate lookup with NO PEPTIDE INPUT: within-allele "
                  "0.5000 by construction, pooled 0.7122",
}

#: Part of the printed name, not a legend entry.
SIGN_FLIP_TAG = "[percentile rank, SIGN-FLIPPED: lower raw = better binder]"

#: A flag shorter or vaguer than this is an asterisk wearing words.
_MIN_FLAG_CHARS = 12
_GENERIC_FLAGS = frozenset(
    {"*", "**", "***", "see note", "note", "caveat", "careful", "n/a", "tbd", "see below"}
)

_RANK_HINTS = ("rank", "percentile", "pctile")


def capability_flag(name: str) -> str | None:
    """The specific reason `name` is not a capability measurement, or None.

    Longest prefix wins: `titanian_ablation_*` must get the ablation's reason, not the
    one written for the headline `titanian` column.
    """
    for key in sorted(NOT_A_CAPABILITY, key=len, reverse=True):
        if name.startswith(key):
            return NOT_A_CAPABILITY[key]
    return None


def looks_like_percentile_rank(name: str) -> bool:
    """Whether `name` names a percentile-rank column, which must be sign-flipped.

    Over-inclusive on purpose: a false positive costs one explicit `sign_flipped=True`
    at the call site, a false negative costs the 0.4650 anti-predictive floor.
    """
    low = name.lower()
    return any(h in low for h in _RANK_HINTS)


def flip_percentile_rank(name: str, values) -> tuple[str, list[float]]:
    """Negate a lower-is-better column and make the negation visible in its name.

    >>> flip_percentile_rank("PRIME_Ranks", [0.1, 0.9])
    ('PRIME_Ranks [percentile rank, SIGN-FLIPPED: lower raw = better binder]', [-0.1, -0.9])
    """
    return f"{name} {SIGN_FLIP_TAG}", [-float(v) for v in values]


@dataclass
class MethodRow:
    """One line of a table.

    `stratified` is a `metrics.StratifiedResult` (or None), `pooled` a float (or None);
    one without the other raises at render time.
    """

    name: str
    stratified: object | None = None
    pooled: float | None = None
    sign_flipped: bool = False
    #: True only when the submitter DECLARED a direction. A declaration is the one thing
    #: that excuses a rank-shaped name carrying unflipped scores.
    direction_declared: bool = False
    flag: str | None = None
    note: str | None = None


@dataclass
class Table:
    """A renderable table. `mde` is mandatory: there is no unqualified ranking here."""

    title: str
    rows: list[MethodRow]
    mde: float | None = DOCUMENTED_MDE
    mde_basis: str = DOCUMENTED_MDE_BASIS
    subtitle: str = ""
    stratum_label: str = "allele"
    footnotes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# rule enforcement


def _check_row(row: MethodRow) -> None:
    has_strat = row.stratified is not None
    has_pool = row.pooled is not None

    if has_strat != has_pool:
        missing = "pooled" if has_strat else "stratified"
        raise ReportError(
            f"{row.name!r} carries "
            f"{'stratified' if has_strat else 'pooled'} without {missing}. Pooled 0.6304 "
            "and within-allele 0.5306 have been produced from the same scores; either "
            "figure alone is read as 'the' AUROC. They ship together or neither ships."
        )
    if not has_strat and not row.note:
        raise ReportError(
            f"{row.name!r} has no figures and no note. A blank cell must say why it is "
            "blank -- 'could not run' is not 'did not pass'."
        )

    if (looks_like_percentile_rank(row.name) and not row.sign_flipped
            and not row.direction_declared):
        raise ReportError(
            f"{row.name!r} reads as a percentile-rank column but is not marked "
            "sign_flipped. Reading *_Ranks with the Scores sign convention produced the "
            "0.4650 'anti-predictive floor'; flipped, those columns return ~0.535."
        )
    if row.sign_flipped and "SIGN-FLIPPED" not in row.name:
        raise ReportError(
            f"{row.name!r} is sign-flipped but the printed name does not say so. An "
            "invisible flip is indistinguishable from a different model."
        )

    if row.flag is not None:
        reason = row.flag.strip()
        if reason.lower() in _GENERIC_FLAGS or len(reason) < _MIN_FLAG_CHARS:
            raise ReportError(
                f"{row.name!r} carries the generic flag {row.flag!r}. A flag states the "
                "specific reason the column is not a capability measurement; a generic "
                "asterisk is read as a minor caveat."
            )


def _rank_key(row: MethodRow) -> float:
    """Ranking reads `stratified` only. An allele lookup with no peptide input scores
    pooled 0.7122 at within-allele exactly 0.5000."""
    return row.stratified.value if row.stratified is not None else float("-inf")


# ---------------------------------------------------------------------------
# rendering


def _wrap(text: str, width: int) -> list[str]:
    """Greedy wrap. Not textwrap: its `break_long_words` would split a sha256 or a DOI."""
    lines, current = [], ""
    for word in str(text).split():
        if current and len(current) + 1 + len(word) > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    if current:
        lines.append(current)
    return lines or [""]


# Counts, not footnotes: on the infectious set only 36 of 78 alleles are admissible,
# covering 1,787 of 2,093 rows, so a bare 0.55 reads as measured on all of them.
def _strata(result) -> str:
    return f"{result.n_strata_contributing}/{result.n_strata_total}"


def _evaluable(result) -> str:
    return (f"{result.n_rows_evaluable}/{result.n_rows_total} "
            f"({result.evaluable_fraction:.1%})")


def render_table(table: Table) -> str:
    """Render, or raise. Every rule is checked before a single line is emitted."""
    if table.mde is None:
        raise ReportError(
            "no minimum detectable difference supplied. A table without an MDE presents "
            "every gap as a result, including gaps of 0.003 at an MDE of 0.059-0.095."
        )
    for row in table.rows:
        _check_row(row)

    scored = sorted((r for r in table.rows if r.stratified is not None),
                    key=_rank_key, reverse=True)
    unscored = [r for r in table.rows if r.stratified is None]

    head = f"within-{table.stratum_label}"
    widths = (4, max(len(head), 13), 9, max(22, max((len(_evaluable(r.stratified))
                                                     for r in scored), default=22)))
    out: list[str] = ["", table.title]
    if table.subtitle:
        out.append(table.subtitle)
    out.append(f"minimum detectable difference: {table.mde:.4f}")
    out += [f"  {line}" for line in _wrap(table.mde_basis, 94)]
    out.append("  gaps smaller than it print '=' and the measured gap, NOT a ranking.")
    out.append(f"estimator: pair-weighted within-{table.stratum_label} AUROC over "
               "admissible strata only.")
    out.append("  Pooled is printed beside it and is never ranked on.")

    # Concentration is a property of the design, not of any one method: hoist it when
    # every scored row shares it, so a reader does not learn to skip 14 copies.
    shared = {(bool(getattr(r.stratified, "provisional", False)),
               round(getattr(r.stratified, "top3_pair_weight", 0.0), 6))
              for r in scored}
    hoisted = len(shared) == 1 and scored and next(iter(shared))[0]
    if hoisted:
        out.append(f"PROVISIONAL, every row: the top three strata hold "
                   f"{scored[0].stratified.top3_pair_weight:.1%} of the pair weight on "
                   "this design.")
    out.append("")
    out.append(f"{'rank':>{widths[0]}}  {head:>{widths[1]}}  {'pooled':>{widths[2]}}  "
               f"{'strata':>6}  {'evaluable rows':<{widths[3]}}  method")
    out.append(f"{'-' * widths[0]:>{widths[0]}}  {'-' * widths[1]}  {'-' * widths[2]}  "
               f"{'-' * 6}  {'-' * widths[3]}  {'-' * 40}")

    previous: float | None = None
    for position, row in enumerate(scored, 1):
        result = row.stratified
        gap = None if previous is None else previous - result.value
        resolvable = gap is None or gap >= table.mde
        label = str(position) if resolvable else "="
        out.append(
            f"{label:>{widths[0]}}  {result.value:>{widths[1]}.4f}  "
            f"{row.pooled:>{widths[2]}.4f}  {_strata(result):>6}  "
            f"{_evaluable(result):<{widths[3]}}  {row.name}"
        )
        if not resolvable:
            out.append(f"{'':>{widths[0]}}  NOT RESOLVABLE against the row above: gap "
                       f"{gap:.4f} < MDE {table.mde:.4f}")
        if row.flag:
            out.append(f"{'':>{widths[0]}}  ^-- NOT A CAPABILITY MEASUREMENT: {row.flag}")
        if getattr(result, "provisional", False) and not hoisted:
            out.append(f"{'':>{widths[0]}}  ^-- PROVISIONAL: the top three strata hold "
                       f"{result.top3_pair_weight:.1%} of the pair weight")
        if row.note:
            out.append(f"{'':>{widths[0]}}  ^-- {row.note}")
        previous = result.value

    for row in unscored:
        out.append(
            f"{'--':>{widths[0]}}  {'n/a':>{widths[1]}}  {'n/a':>{widths[2]}}  "
            f"{'n/a':>6}  {'n/a':<{widths[3]}}  {row.name}"
        )
        out.append(f"{'':>{widths[0]}}  ^-- NOT SCORED (this is not a low score): "
                   f"{row.note}")
        if row.flag:
            out.append(f"{'':>{widths[0]}}  ^-- NOT A CAPABILITY MEASUREMENT: {row.flag}")

    out.append("")
    out.append("pooled is reported, never ranked on: an allele positive-rate lookup with "
               "no peptide")
    out.append("input reaches pooled 0.7122 at within-allele exactly 0.5000 by "
               "construction.")
    for note in table.footnotes:
        out.append(note)
    out.append("")
    return "\n".join(out)


def _or_unknown(value, fmt: str = "{}") -> str:
    """An absent count prints `unknown`: `0 of 2,093` must never mean "nothing looked"."""
    return "unknown" if value is None else fmt.format(value)


#: `excluded` is a licence forbidding the rows to leave their source, listed rather than
#: hidden so the reason is visible before an evaluation is written against it.
_STATUS = {"ship": "ships", "fetch": "fetch required",
           "excluded": "excluded (licence)"}


def _clip(text: str, width: int) -> str:
    """Truncate to `width` with an unmistakable marker: a clipped licence must not read
    as the whole of one."""
    text = str(text)
    return text if len(text) <= width else text[:width - 4] + " ..."


_LICENCE_WIDTH = 34


def render_row_sets(specs: list[dict]) -> str:
    """The `list` table: size, class balance, licence, and ships-or-fetch."""
    name_w = max([len("name")] + [len(str(s.get("name", "?"))) for s in specs])
    out = ["", "ROW SETS", ""]
    out.append(f"{'name':<{name_w}}  {'rows':>7}  {'positive':>13}  "
               f"{'licence':<{_LICENCE_WIDTH}}  {'availability':<18}  note")
    out.append(f"{'-' * name_w}  {'-' * 7}  {'-' * 13}  {'-' * _LICENCE_WIDTH}  "
               f"{'-' * 18}  {'-' * 24}")
    for spec in specs:
        n = spec.get("n_rows")
        n_pos = spec.get("n_pos")
        balance = (f"{n_pos} ({n_pos / n:.1%})" if n and n_pos is not None else
                   _or_unknown(n_pos))
        availability = spec.get("availability") or _STATUS.get(spec.get("status"))
        if availability is None:
            ships = spec.get("ships")
            availability = ("unknown" if ships is None else
                            "ships" if ships else "fetch required")
        out.append(f"{str(spec.get('name', '?')):<{name_w}}  "
                   f"{_or_unknown(n, '{:>7}'):>7}  {balance:>13}  "
                   f"{_clip(spec.get('licence') or 'unknown', _LICENCE_WIDTH):<{_LICENCE_WIDTH}}  "
                   f"{availability:<18}  {spec.get('note', '') or ''}")
    out.append("")
    out.append("`fetch required` is a licence consequence, not a packaging oversight: "
               "TESLA, VDJdb,")
    out.append("McPAS-TCR and the NetMHCpan score columns may not be redistributed from "
               "here.")
    out.append("")
    return "\n".join(out)


#: The three statuses a check can end in. `could-not-run` is NOT a pass; it exits 2.
PASS, FAIL, UNRUNNABLE = "pass", "FAIL", "could-not-run"


def render_validation(path: str, checks: list[tuple[str, str, str]]) -> str:
    """Format the `validate` result. `checks` is [(name, status, detail), ...]."""
    out = ["", f"FORMAT CHECK  {path}", ""]
    for name, status, detail in checks:
        out.append(f"  [{status:>13}]  {name}")
        if detail:
            out.append(f"                   {detail}")
    n_fail = sum(1 for _, s, _ in checks if s == FAIL)
    n_unrun = sum(1 for _, s, _ in checks if s == UNRUNNABLE)
    out.append("")
    out.append(f"{len(checks)} checks, {n_fail} failed, {n_unrun} could not run")
    if n_unrun:
        out.append("a check that could not run is NOT a check that passed; this exits 2.")
    out.append("")
    return "\n".join(out)
