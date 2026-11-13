"""Regression fixtures for the command-line interface and the formatter.

Governing rule in `test_canon_and_metrics.py`. Nothing here touches the network, the
data tree, or `epigauge fixtures` (the fixture runner launches this file). Standalone:
`python3 -I epigauge/tests/test_cli.py`.
"""

from __future__ import annotations

import contextlib
import io
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from epigauge import cli, report  # noqa: E402
from epigauge.metrics import StratifiedResult  # noqa: E402
from epigauge.report import FAIL, PASS, UNRUNNABLE, MethodRow, ReportError, Table  # noqa: E402

FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def expect_raises(exc, fn, *a, **kw) -> bool:
    try:
        fn(*a, **kw)
    except exc:
        return True
    except Exception as e:  # noqa: BLE001
        FAILURES.append(f"expected {exc.__name__}, got {type(e).__name__}: {e}")
        return False
    FAILURES.append(f"expected {exc.__name__}, nothing raised: {getattr(fn, '__name__', fn)}")
    return False


# The measured shape of the infectious set: 36 of 78 alleles admissible, 1,787 of
# 2,093 rows. Verbatim, so the printed counts are checkable by eye.
def _strat(value: float) -> StratifiedResult:
    return StratifiedResult(value=value, n_strata_contributing=36, n_strata_total=78,
                            n_rows_evaluable=1787, n_rows_total=2093,
                            herfindahl=0.0912, top3_pair_weight=0.49, provisional=False)


def _status(checks, name):
    for check_name, status, _detail in checks:
        if check_name.startswith(name):
            return status
    FAILURES.append(f"no check named {name!r} in {[c[0] for c in checks]}")
    return None


def _detail(checks, name):
    for check_name, _status, detail in checks:
        if check_name.startswith(name):
            return detail
    return ""


def _write(directory, text: str) -> str:
    path = Path(directory) / "submission.csv"
    path.write_text(text, encoding="utf-8")
    return str(path)


_GOOD = ("# direction: higher-is-better\n"
         "row_id,score\nr1,0.9\nr2,0.1\nr3,0.5\n")
_IDS = ["r1", "r2", "r3"]


def test_submission_missing_rows_is_rejected() -> None:
    """A partial join is not a worse score, it is a different population -- the
    (peptide, allele) join that returned 0 of 833 is the extreme case."""
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "# direction: higher-is-better\nrow_id,score\nr1,0.9\nr2,0.1\n")
        sub = cli.read_submission(path)
        checks = cli.validate_submission(sub, _IDS)
    check(_status(checks, "every benchmark row_id") == FAIL,
          "a submission missing r3 must FAIL the coverage check")
    check("missing" in _detail(checks, "every benchmark row_id"),
          f"the detail must name the missing rows: {_detail(checks, 'every benchmark row_id')!r}")
    check(cli.exit_code_for(checks) == cli.EXIT_CHECK_FAILED,
          f"missing rows must exit 1, got {cli.exit_code_for(checks)}")


def test_duplicated_and_unknown_row_ids_are_rejected() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "# direction: higher-is-better\n"
                           "row_id,score\nr1,0.9\nr1,0.8\nr2,0.1\nr3,0.5\nr9,0.4\n")
        checks = cli.validate_submission(cli.read_submission(path), _IDS)
    detail = _detail(checks, "every benchmark row_id")
    check(_status(checks, "every benchmark row_id") == FAIL,
          "a duplicated row_id plus an unknown row_id must FAIL")
    check("duplicated" in detail and "not in the benchmark" in detail,
          f"both problems must be named, got {detail!r}")


def test_submission_without_declared_direction_is_rejected() -> None:
    """Percentile-rank columns read with the Scores sign convention made a whole family
    of tools look anti-predictive at 0.4650; flipped they return ~0.535."""
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "row_id,score\nr1,0.9\nr2,0.1\nr3,0.5\n")
        sub = cli.read_submission(path)
        checks = cli.validate_submission(sub, _IDS)
    check(sub.direction is None, f"no direction should have been inferred, got {sub.direction!r}")
    check(_status(checks, "a score direction") == FAIL,
          "an undeclared score direction must FAIL")
    check(cli.exit_code_for(checks) == cli.EXIT_CHECK_FAILED,
          "an undeclared direction must exit 1")


def test_unrecognised_direction_is_rejected_rather_than_defaulted() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "# direction: whichever\nrow_id,score\nr1,0.9\nr2,0.1\nr3,0.5\n")
        sub = cli.read_submission(path)
        checks = cli.validate_submission(sub, _IDS)
    check(sub.direction is None, "an unrecognised direction must not resolve")
    check(_status(checks, "a score direction") == FAIL,
          "an unrecognised direction must FAIL rather than default to higher-is-better")
    check("whichever" in _detail(checks, "a score direction"),
          "the rejected value must be echoed back")


def test_declared_lower_is_better_flips_and_says_so() -> None:
    """The direction is declared, never inferred, and a flip is visible in the name."""
    name, values = cli._oriented("mytool", cli.LOWER, [0.1, 0.9])
    check(values == [-0.1, -0.9], f"lower-is-better must be negated, got {values}")
    check("SIGN-FLIPPED" in name, f"the flip must be visible in the name, got {name!r}")
    same, unchanged = cli._oriented("mytool", cli.HIGHER, [0.1, 0.9])
    check(unchanged == [0.1, 0.9] and same == "mytool",
          "higher-is-better must pass through untouched")


def test_rank_shaped_column_without_a_flip_is_refused() -> None:
    row = MethodRow(name="PRIME-2.0_Ranks", stratified=_strat(0.4650), pooled=0.4700)
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[row], mde=0.04))
    # Flipped and tagged, it renders.
    flipped, _ = report.flip_percentile_rank("PRIME-2.0_Ranks", [0.1, 0.9])
    ok = MethodRow(name=flipped, stratified=_strat(0.5350), pooled=0.5400,
                   sign_flipped=True)
    text = report.render_table(Table(title="t", rows=[ok], mde=0.04))
    check("SIGN-FLIPPED" in text, "the printed table must carry the flip tag")


def test_invisible_flip_is_refused() -> None:
    row = MethodRow(name="quietly flipped", stratified=_strat(0.55), pooled=0.60,
                    sign_flipped=True)
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[row], mde=0.04))


def test_all_constant_submission_is_rejected() -> None:
    """A constant column scores exactly 0.5 under the mid-rank correction, which reads
    as a measured null rather than as an empty submission."""
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "# direction: higher-is-better\n"
                           "row_id,score\nr1,0.5\nr2,0.5\nr3,0.5\n")
        checks = cli.validate_submission(cli.read_submission(path), _IDS)
    check(_status(checks, "the submission is not all-constant") == FAIL,
          "an all-constant submission must FAIL")
    check("0.5" in _detail(checks, "the submission is not all-constant"),
          "the detail must say what a constant column scores")
    check(cli.exit_code_for(checks) == cli.EXIT_CHECK_FAILED,
          "an all-constant submission must exit 1")


def test_non_finite_scores_are_rejected() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "# direction: higher-is-better\n"
                           "row_id,score\nr1,0.9\nr2,nan\nr3,oops\n")
        checks = cli.validate_submission(cli.read_submission(path), _IDS)
    check(_status(checks, "every score is a finite number") == FAIL,
          "NaN and unparseable scores must FAIL")


def test_a_clean_submission_passes() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        checks = cli.validate_submission(cli.read_submission(_write(tmp, _GOOD)), _IDS)
    check(all(s == PASS for _, s, _ in checks),
          f"a clean submission must pass everything: {[(c, s) for c, s, _ in checks if s != PASS]}")
    check(cli.exit_code_for(checks) == cli.EXIT_OK, "a clean submission must exit 0")


def test_formatter_raises_when_pooled_is_omitted() -> None:
    """Pooled 0.6304 and within-allele 0.5306 come from the same scores; either alone
    reads as "the" AUROC, and the gap between them is the finding."""
    row = MethodRow(name="a method", stratified=_strat(0.5306), pooled=None)
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[row], mde=0.04))


def test_formatter_raises_when_stratified_is_omitted() -> None:
    row = MethodRow(name="a method", stratified=None, pooled=0.6304)
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[row], mde=0.04))


def test_a_blank_row_must_say_why_it_is_blank() -> None:
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[MethodRow(name="a method")], mde=0.04))
    text = report.render_table(Table(
        title="t", mde=0.04,
        rows=[MethodRow(name="a method", note="all rows share one allele")]))
    check("NOT SCORED (this is not a low score)" in text,
          "an unscored row must not read as a bottom-ranked row")


def test_the_table_is_ranked_on_stratified_and_never_on_pooled() -> None:
    """Ranking on pooled puts an allele lookup with no peptide input first: pooled ~0.71
    at within-allele exactly 0.5000 by construction."""
    lookup = MethodRow(name="allele prior (no peptide input)",
                       stratified=_strat(0.5000), pooled=0.7100)
    model = MethodRow(name="a real model", stratified=_strat(0.5540), pooled=0.6169)
    text = report.render_table(Table(title="t", rows=[lookup, model], mde=0.04))
    check(text.index("a real model") < text.index("allele prior"),
          "the higher pooled figure was ranked first; the sort must read stratified")
    check("0.7100" in text, "pooled must still be printed, just never ranked on")


def test_every_stratified_figure_prints_its_strata_and_row_counts() -> None:
    """A bare 0.55 read as measured on 78 alleles and 2,093 rows, when only 36 alleles
    are admissible and they cover 1,787 rows."""
    text = report.render_table(Table(
        title="t", mde=0.04,
        rows=[MethodRow(name="a method", stratified=_strat(0.5540), pooled=0.6169)]))
    check("36/78" in text, "the admissible-stratum count must be printed beside the figure")
    check("1787/2093" in text, "the evaluable-row count must be printed beside the figure")
    check("85.4%" in text, f"the evaluable fraction must be printed; got:\n{text}")


def test_a_table_without_an_mde_is_refused() -> None:
    """Ordering two methods 0.0036 apart when the MDE is 0.059-0.095 reports resampling
    noise as a ranking."""
    row = MethodRow(name="a method", stratified=_strat(0.5540), pooled=0.6169)
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[row], mde=None))


def test_a_gap_below_the_mde_prints_not_resolvable_instead_of_a_rank() -> None:
    top = MethodRow(name="method alpha", stratified=_strat(0.5531), pooled=0.6304)
    near = MethodRow(name="method beta", stratified=_strat(0.5495), pooled=0.5912)
    far = MethodRow(name="method gamma", stratified=_strat(0.4800), pooled=0.5100)
    text = report.render_table(Table(title="t", rows=[top, near, far], mde=0.04))
    check("NOT RESOLVABLE" in text,
          "a 0.0036 gap at an MDE of 0.04 must print NOT RESOLVABLE")
    check("0.0036" in text, "the measured gap must be printed, not just the verdict")
    beta_line = [ln for ln in text.splitlines() if "method beta" in ln][0]
    check(beta_line.strip().startswith("="),
          f"an unresolvable row must not carry a rank number: {beta_line!r}")
    gamma_line = [ln for ln in text.splitlines() if "method gamma" in ln][0]
    check(gamma_line.strip().startswith("3"),
          f"a resolvable row keeps its rank number: {gamma_line!r}")
    check("minimum detectable difference" in text.splitlines()[3].lower()
          or any("minimum detectable difference" in ln.lower()
                 for ln in text.splitlines()[:6]),
          "the MDE must be printed at the TOP of the table, not in a footnote")


def test_a_generic_flag_is_refused_and_a_specific_reason_prints() -> None:
    """A reference column marked with an asterisk reads as a minor caveat; the label is
    recoverable from the IEDB assay table at 0.9101, 100% coverage."""
    vague = MethodRow(name="TITANiAN-IM", stratified=_strat(0.84), pooled=0.88, flag="*")
    expect_raises(ReportError, report.render_table,
                  Table(title="t", rows=[vague], mde=0.04))

    reason = report.capability_flag("TITANiAN-IM")
    check(reason is not None and "0.9101" in reason,
          f"the registry must carry the measured reason, got {reason!r}")
    row = MethodRow(name="TITANiAN-IM", stratified=_strat(0.84), pooled=0.88, flag=reason)
    text = report.render_table(Table(title="t", rows=[row], mde=0.04))
    check("NOT A CAPABILITY MEASUREMENT" in text and "0.9101" in text,
          "the specific reason must appear in the table")


def test_could_not_run_is_not_a_pass() -> None:
    """"could not run" reported as "passed": `0 of 2,093` must not be able to mean
    "nothing looked"."""
    with tempfile.TemporaryDirectory() as tmp:
        sub = cli.read_submission(_write(tmp, _GOOD))
    checks = cli.validate_submission(sub, expected_ids=None)
    check(_status(checks, "every benchmark row_id") == UNRUNNABLE,
          "with no benchmark row set loaded, coverage must be could-not-run")
    check(cli.exit_code_for(checks) == cli.EXIT_COULD_NOT_RUN,
          f"could-not-run must exit 2, got {cli.exit_code_for(checks)}")
    text = report.render_validation(sub.path, checks)
    check("NOT a check that passed" in text,
          "the printed report must say that a could-not-run is not a pass")


def test_a_failure_outranks_a_could_not_run() -> None:
    checks = [("a", FAIL, ""), ("b", UNRUNNABLE, ""), ("c", PASS, "")]
    check(cli.exit_code_for(checks) == cli.EXIT_CHECK_FAILED,
          "a FAIL must outrank a could-not-run")


def test_missing_input_exits_two_not_zero() -> None:
    missing = str(Path(tempfile.gettempdir()) / "epigauge-no-such-submission-7f3a.csv")
    check(not Path(missing).exists(), "the fixture path must not exist")
    for argv in (["validate", missing], ["score", missing, "--set", "titanian_infectious"]):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = cli.main(argv)
        check(code == cli.EXIT_COULD_NOT_RUN,
              f"{argv[0]} on a missing file returned {code}, expected "
              f"{cli.EXIT_COULD_NOT_RUN}")
        check(code != cli.EXIT_OK, f"{argv[0]} must never exit 0 on an absent input")
        check("COULD NOT RUN" in err.getvalue(),
              f"{argv[0]} must say it could not run: {err.getvalue()!r}")


def test_a_malformed_file_is_a_failed_check_not_a_could_not_run() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = _write(tmp, "peptide,allele,score\nSIINFEKL,HLA-A*02:01,0.9\n")
        expect_raises(cli.SubmissionError, cli.read_submission, path)
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = cli.main(["validate", path])
        check(code == cli.EXIT_CHECK_FAILED,
              f"a file with no row_id column must exit 1, got {code}")


SUBCOMMANDS = ("list", "fetch", "score", "leaderboard", "validate", "fixtures")


def test_help_works_for_every_subcommand() -> None:
    parser = cli.build_parser()
    for command in SUBCOMMANDS:
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                parser.parse_args([command, "--help"])
        except SystemExit as exit_code:
            check(exit_code.code == 0,
                  f"`epigauge {command} --help` exited {exit_code.code}")
        else:
            FAILURES.append(f"`epigauge {command} --help` did not exit")
        check(command in out.getvalue(),
              f"`epigauge {command} --help` printed nothing about {command}")

    out = io.StringIO()
    try:
        with contextlib.redirect_stdout(out):
            parser.parse_args(["--help"])
    except SystemExit:
        pass
    for command in SUBCOMMANDS:
        check(command in out.getvalue(), f"top-level --help omits {command}")


def test_every_subcommand_degrades_rather_than_crashing_without_data() -> None:
    """data.py, baselines.py and bootstrap.py are written alongside this file."""
    if cli.data is not None:
        return  # data.py has landed; this fixture is about the window before it does
    for argv in (["list"], ["fetch", "titanian_infectious"], ["leaderboard", "--set", "titanian_infectious"]):
        with contextlib.redirect_stderr(io.StringIO()) as err, \
                contextlib.redirect_stdout(io.StringIO()):
            code = cli.main(argv)
        check(code == cli.EXIT_COULD_NOT_RUN,
              f"{argv[0]} without data.py returned {code}, expected could-not-run")
        check("epigauge.data" in err.getvalue(),
              f"{argv[0]} must name the module it needed: {err.getvalue()!r}")


def test_the_mde_falls_back_to_the_documented_band_and_says_so() -> None:
    """The fallback MDE must be the conservative end of the measured 0.059-0.095 band:
    an assumed MDE below the real one lets an unresolvable margin be declared resolvable."""
    class _Args:
        proteome = None
        replicates = 100

    # arms=[] and no --proteome, so the measured path cannot run and must fall back.
    value, basis = cli._mde("titanian_infectious", [1, 0, 1, 0], ["a", "a", "b", "b"],
                            ["AAAAAAAAA"] * 4, [], _Args())
    check(0.0 < value <= 0.095,
          f"the fallback MDE must sit inside the measured 0.059-0.095 band, got {value}")
    check(value >= 0.059,
          f"the fallback must be the CONSERVATIVE end of the measured band, got {value}. "
          "An assumed MDE smaller than the real one lets an unresolvable margin be "
          "declared resolvable, which is the failure this constant prevents.")
    if cli.bootstrap is None:
        check("assumed" in basis and "NOT computed" in basis,
              f"an assumed MDE must say it was not computed: {basis!r}")


def test_a_licence_refusal_is_a_could_not_run_not_an_empty_table() -> None:
    """A set nobody may redistribute must not come back as zero rows and a clean exit."""
    if cli.data is None:
        return  # exercised once data.py is importable
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = cli.main(["leaderboard", "--set", "vdjdb"])
    check(code == cli.EXIT_COULD_NOT_RUN,
          f"an excluded set must exit {cli.EXIT_COULD_NOT_RUN}, got {code}")
    check("Licence" in err.getvalue() or "licence" in err.getvalue(),
          f"the refusal must carry the licence reason: {err.getvalue()[:200]!r}")


def test_an_ambiguous_short_set_name_is_refused_rather_than_picked() -> None:
    if cli.data is None:
        return
    check(cli._resolve_set("infectious") == "titanian_infectious",
          "a unique short name must resolve")
    expect_raises(cli.CouldNotRun, cli._resolve_set, "titanian")
    expect_raises(cli.CouldNotRun, cli._resolve_set, "no-such-set")



_STANDALONE = False


def test_zz_every_fixture_above_reported_no_failures() -> None:
    """Turns failures recorded by `check()` into a PYTEST failure; without it pytest
    stays green whatever was recorded. Last in the file so it runs last."""
    if _STANDALONE:
        return  # main() prints and returns the same information
    assert not FAILURES, "\n".join(["recorded failures:"] + [f"  - {f}" for f in FAILURES])

def main() -> int:
    global _STANDALONE
    _STANDALONE = True
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        before = len(FAILURES)
        try:
            t()
        except Exception as e:  # noqa: BLE001
            FAILURES.append(f"{t.__name__} raised {type(e).__name__}: {e}")
        status = "ok" if len(FAILURES) == before else "FAIL"
        print(f"  [{status:>4}] {t.__name__}")
    print(f"\n{len(tests)} fixtures, {len(FAILURES)} failure(s)")
    for f in FAILURES:
        print(f"  - {f}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
