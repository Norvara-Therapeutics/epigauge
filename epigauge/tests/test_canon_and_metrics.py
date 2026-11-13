"""Regression fixtures for documented instrument failures.

Each test reproduces a measurement this project got confidently wrong and asserts the
current code refuses it. The governing rule for every test module here: a guard that
cannot be DEMONSTRATED to fire is not a guard. Run under pytest, or standalone with
`python3 -I epigauge/tests/test_canon_and_metrics.py`.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from epigauge.canon import (  # noqa: E402
    CanonError,
    canonicalise_allele,
    canonicalise_organism,
    canonicalise_peptide,
)
from epigauge.metrics import (  # noqa: E402
    DegenerateGroupingError,
    auroc,
    mixture_decomposition,
    pooled_rank_auroc,
    stratified_auroc,
)

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
    FAILURES.append(f"expected {exc.__name__}, nothing raised: {fn.__name__}")
    return False


def test_allele_join_reproduces_contamination_trap() -> None:
    """The allele-format join that returned 0 of 833 and hid a contaminated benchmark:
    im_test.csv spells HLA-A*01:01, the TITANiAN CSV spells A0101, same table."""
    bigmhc_style = ["HLA-A*01:01", "HLA-A*02:01", "HLA-B*07:02"]
    titanian_style = ["A0101", "A0201", "B0702"]

    naive = set(bigmhc_style) & set(titanian_style)
    check(len(naive) == 0,
          "the trap no longer reproduces; the raw spellings should not intersect")

    canon_a = {canonicalise_allele(x) for x in bigmhc_style}
    canon_b = {canonicalise_allele(x) for x in titanian_style}
    check(canon_a == canon_b,
          f"canonicalisation failed to recover the join: {canon_a} vs {canon_b}")
    check(len(canon_a & canon_b) == 3,
          "all three pairs must match after canonicalisation")

    for variant in ("hla-a0101", "A*01:01", "HLA-A01:01", "HLA-A*01:01:01:02N", " A0101 "):
        check(canonicalise_allele(variant) == "HLA-A*01:01",
              f"variant {variant!r} -> {canonicalise_allele(variant)!r}, expected HLA-A*01:01")


def test_allele_raises_rather_than_falling_back() -> None:
    """A missed allele join produced exactly 0.5000, caught only because 0.5000 looks
    impossible; a plausible fallback value would have been published."""
    for bad in ("", "   ", "not-an-allele", "HLA-", "A", "???"):
        expect_raises(CanonError, canonicalise_allele, bad)
    expect_raises(CanonError, canonicalise_allele, None)


def test_organism_folds_one_parenthetical_only() -> None:
    """One organism audited as two blocks: 8,857 rows at 48.4% positive and 7,492 at
    31.8%."""
    check(canonicalise_organism("Homo sapiens") == canonicalise_organism("Homo sapiens (human)"),
          "the two human spellings must fold to one key")
    check(canonicalise_organism("Homo sapiens") == "homo sapiens",
          f"unexpected human key: {canonicalise_organism('Homo sapiens')!r}")


def test_organism_does_not_fold_viruses_named_human() -> None:
    """The substring fix for that split would fold 13 viruses named "human" into the
    human key."""
    human = canonicalise_organism("Homo sapiens")
    viruses = [
        "Human betaherpesvirus 5",
        "Human immunodeficiency virus 1",
        "Human papillomavirus type 16",
        "Human T-cell leukemia virus type I",
        "Human adenovirus 5",
    ]
    for v in viruses:
        got = canonicalise_organism(v)
        check(got != human,
              f"{v!r} folded into the human key -- substring matching has crept back in")
        check("human" in got,
              f"{v!r} -> {got!r}; the word should survive, only the key must stay distinct")
    check(len({canonicalise_organism(v) for v in viruses}) == len(viruses),
          "the five viruses must remain five distinct keys")


def test_organism_does_not_merge_strain_spellings() -> None:
    """Merging strains is a judgement, so it is reported, not applied."""
    a = canonicalise_organism("Lymphocytic choriomeningitis virus (Armstrong)")
    b = canonicalise_organism("Lymphocytic choriomeningitis virus (WE)")
    check(a == b, "a trailing strain parenthetical is the one thing that IS folded")
    c = canonicalise_organism("Influenza A virus")
    d = canonicalise_organism("Influenza A virus H1N1")
    check(c != d, "distinct spellings must NOT be silently merged")


_Y = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
_S = [0.9, 0.1, 0.8, 0.7, 0.6, 0.2, 0.95, 0.3, 0.55, 0.4, 0.85, 0.15]
_G = ["A", "A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"]


def test_missing_grouping_raises_and_would_have_equalled_pooled() -> None:
    """A uniformly-missing allele column makes within-allele equal pooled exactly
    (0.584138 vs the correct 0.492317). Both look publishable, so this raises."""
    one_group = ["A"] * len(_Y)
    expect_raises(DegenerateGroupingError, stratified_auroc, _Y, _S, one_group)

    silent = pooled_rank_auroc(_Y, _S, one_group)
    pooled = auroc(_Y, _S)
    check(abs(silent - pooled) < 1e-12,
          f"expected the degenerate case to equal pooled exactly; got {silent} vs {pooled}")

    correct = stratified_auroc(_Y, _S, _G)
    check(abs(correct.value - pooled) > 1e-6,
          "the stratified and pooled estimates must differ on this fixture, or it proves nothing")


def test_all_singleton_grouping_raises_instead_of_returning_half() -> None:
    """Bypassed, all-singleton grouping returns exactly 0.500000 -- which reads as
    "no signal" rather than "no grouping"."""
    singles = [f"g{i}" for i in range(len(_Y))]
    expect_raises(DegenerateGroupingError, stratified_auroc, _Y, _S, singles)
    check(abs(pooled_rank_auroc(_Y, _S, singles) - 0.5) < 1e-12,
          "the bypassed path should return exactly 0.5 on all-singleton grouping")


def test_empty_group_key_raises() -> None:
    g = list(_G)
    g[3] = ""
    expect_raises(DegenerateGroupingError, stratified_auroc, _Y, _S, g)
    g[3] = None
    expect_raises(DegenerateGroupingError, stratified_auroc, _Y, _S, g)


def test_single_class_strata_are_excluded_not_scored_half() -> None:
    y = [1, 0, 1, 0, 1, 1, 1]
    s = [0.9, 0.1, 0.8, 0.2, 0.5, 0.6, 0.7]
    g = ["A", "A", "A", "A", "Z", "Z", "Z"]  # Z is all-positive
    r = stratified_auroc(y, s, g)
    check(r.n_strata_contributing == 1, f"expected 1 contributing stratum, got {r.n_strata_contributing}")
    check(r.n_strata_total == 2, f"expected 2 total strata, got {r.n_strata_total}")
    check(r.per_stratum["Z"]["pair_weight"] == 0.0, "single-class stratum must carry weight 0")
    check(r.per_stratum["Z"]["auroc"] is None, "single-class stratum must not be scored at 0.5")
    check(r.n_rows_evaluable == 4, f"expected 4 evaluable rows, got {r.n_rows_evaluable}")


def test_no_admissible_stratum_raises() -> None:
    y = [1, 1, 0, 0]
    s = [0.9, 0.8, 0.2, 0.1]
    g = ["A", "A", "B", "B"]  # each stratum single-class
    expect_raises(DegenerateGroupingError, stratified_auroc, y, s, g)


def test_constant_score_is_exactly_half_not_zero_or_one() -> None:
    check(abs(auroc([1, 0, 1, 0], [0.5] * 4) - 0.5) < 1e-12,
          "a constant score must give exactly 0.5 under mid-ranks")


def test_mixture_flags_pooled_exceeding_every_constituent() -> None:
    """A pooled figure above BOTH constituents -- measured pathogen 0.5806, self 0.5727,
    pooled 0.5837. Here blocks are 0.375 and 0.625, pooled 24/36, excess +0.041667."""
    y = [1, 1, 1, 1, 0, 0] + [1, 1, 0, 0, 0, 0]
    s = [0.90, 0.80, 0.70, 0.60, 0.85, 0.75] + [0.40, 0.20, 0.35, 0.30, 0.25, 0.15]
    b = ["P"] * 6 + ["S"] * 6

    m = mixture_decomposition(y, s, b)
    check(abs(m.blocks["P"]["auroc"] - 0.375) < 1e-12,
          f"block P within-AUROC {m.blocks['P']['auroc']}, expected 0.375")
    check(abs(m.blocks["S"]["auroc"] - 0.625) < 1e-12,
          f"block S within-AUROC {m.blocks['S']['auroc']}, expected 0.625")
    check(abs(m.pooled - 24 / 36) < 1e-12, f"pooled {m.pooled}, expected {24/36}")
    check(abs(m.max_constituent - 0.625) < 1e-12,
          f"max constituent {m.max_constituent}, expected 0.625")
    check(abs(m.excess - (24 / 36 - 0.625)) < 1e-12,
          f"excess {m.excess:+.6f}, expected {24/36 - 0.625:+.6f}")
    check(m.excess > 0, f"expected positive excess, got {m.excess:+.6f}")
    check(m.artifact, f"excess {m.excess:+.6f} should have been flagged as an artifact")
    check(m.pooled > m.max_constituent,
          f"pooled {m.pooled:.6f} should exceed max constituent {m.max_constituent:.6f}")

    # The stratified estimator is diagonal-only, so it cannot show the same inflation.
    strat = stratified_auroc(y, s, b)
    check(strat.value <= m.max_constituent + 1e-12,
          f"stratified {strat.value:.6f} exceeded max constituent {m.max_constituent:.6f}; "
          "the estimator is admitting cross-block pairs")


def test_mixture_reports_no_artifact_when_blocks_are_exchangeable() -> None:
    y = [1, 0] * 10
    s = [0.9, 0.1] * 10
    b = (["P"] * 10) + (["S"] * 10)
    m = mixture_decomposition(y, s, b)
    check(not m.artifact,
          f"exchangeable blocks must not be flagged; excess was {m.excess:+.6f}")


def test_peptide_raises_on_non_canonical_rather_than_masking() -> None:
    """An absent field read as a real category: all 2,093 rows declared "self-derived"
    because the default was a category rather than an error."""
    for bad in ("SIINFEKX", "SIIN FEK*", "siinfekl1"):
        expect_raises(CanonError, canonicalise_peptide, bad)
    expect_raises(CanonError, canonicalise_peptide, "SIINFE")        # too short
    expect_raises(CanonError, canonicalise_peptide, "S" * 30)        # too long
    check(canonicalise_peptide("siinfekl") == "SIINFEKL", "happy path broken")



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
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
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
