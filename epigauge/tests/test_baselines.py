"""Regression fixtures for the baseline roster.

Governing rule in `test_canon_and_metrics.py`. Run with pytest, or standalone:

    python3 -I epigauge/tests/test_baselines.py
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from epigauge.baselines import (  # noqa: E402
    AlleleOnly,
    BaselineError,
    EvaluationLeakError,
    Hydropathy,
    HydropathyLengthAllele,
    Length,
    PWM,
    ShuffledSentinel,
    SingularDesignError,
    mean_kd,
    ols,
    residualise,
    sentinel_distribution,
    within_allele_rank,
)
from epigauge.metrics import auroc, stratified_auroc  # noqa: E402

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
    name = getattr(fn, "__name__", repr(fn))
    FAILURES.append(f"expected {exc.__name__}, nothing raised: {name}")
    return False


# Alleles with very different positive rates: the condition the infectious benchmark
# has (137 alleles, pooled allele-lookup ~0.70) and the one the ordering defect needs.
ALLELES = {
    "HLA-A*02:01": 0.85, "HLA-A*01:01": 0.70, "HLA-B*07:02": 0.55,
    "HLA-B*44:02": 0.40, "HLA-C*07:01": 0.25, "HLA-C*04:01": 0.12,
}
AA = "ACDEFGHIKLMNPQRSTVWY"


def make_rows(n: int, *, seed: int, prefix: str, kd_effect: float = 0.0,
              lengths=(8, 9, 10, 11)) -> list[dict]:
    """Rows with an allele-label association; `kd_effect` tilts the positive probability
    by mean Kyte-Doolittle hydropathy for the residualisation fixtures."""
    rng = random.Random(seed)
    names = sorted(ALLELES)
    rows = []
    for i in range(n):
        allele = names[i % len(names)]
        L = lengths[i % len(lengths)]
        pep = "".join(rng.choice(AA) for _ in range(L))
        p = ALLELES[allele] + kd_effect * mean_kd(pep)
        p = min(0.98, max(0.02, p))
        rows.append({
            "record_id": f"{prefix}-{i:05d}",
            "peptide": pep,
            "mhc_allele": allele,
            "label": 1 if rng.random() < p else 0,
        })
    return rows


FIT_ROWS = make_rows(900, seed=11, prefix="fit")
EVAL_ROWS = make_rows(450, seed=97, prefix="eval")
EVAL_IDS = frozenset(r["record_id"] for r in EVAL_ROWS)
Y = [r["label"] for r in EVAL_ROWS]
G = [r["mhc_allele"] for r in EVAL_ROWS]


def pooled(scores) -> float:
    return auroc(Y, scores)


def within(scores) -> float:
    return stratified_auroc(Y, scores, G).value


def test_shuffle_ordering_after_the_encoding_fakes_a_result() -> None:
    """Permuting a sentinel after the allele target encoding fakes a model-shaped number:
    permute-after 0.6083 pooled / 0.5120 within, permute-before 0.4448 / 0.4751."""
    before = sentinel_distribution(
        HydropathyLengthAllele, FIT_ROWS, EVAL_ROWS, pooled,
        eval_ids=EVAL_IDS, seeds=range(30), order="before")
    after = sentinel_distribution(
        HydropathyLengthAllele, FIT_ROWS, EVAL_ROWS, pooled,
        eval_ids=EVAL_IDS, seeds=range(30), order="after")

    # Assertions are on means: the correct sentinel has sd ~0.09, so a single draw of
    # 0.62 or 0.34 is available from the same null.
    check(before.n == 30 and after.n == 30,
          f"expected 30 draws each, got {before.n} and {after.n}")
    check(before.sd > after.sd,
          f"the correctly ordered null should be the wide one (sd {before.sd:.4f} "
          f"against {after.sd:.4f}); the defect is narrow because it is tracking a "
          "real allele prior rather than a permutation")

    check(abs(before.mean - 0.5) < 0.05,
          f"correctly ordered sentinel is not a null: pooled mean {before.mean:.4f} "
          f"(sd {before.sd:.4f}); it must sit at chance")
    check(after.mean - before.mean > 0.08,
          f"the ordering defect no longer reproduces: permute-after {after.mean:.4f} "
          f"against permute-before {before.mean:.4f}. If this gap has closed, either "
          "the target encoding stopped surviving the permutation or the fixture lost "
          "its allele-label association -- check the fixture before deleting the test.")
    check(after.mean > 0.58,
          f"permute-after pooled mean {after.mean:.4f}; the measured defect is 0.6083, "
          "so anything near chance means the encoding is no longer leaking through")

    # The within-allele axis is where the defect is quiet: 0.5120 reads as a clean null.
    after_within = sentinel_distribution(
        HydropathyLengthAllele, FIT_ROWS, EVAL_ROWS, within,
        eval_ids=EVAL_IDS, seeds=range(30), order="after")
    check(abs(after_within.mean - 0.5) < 0.06,
          f"permute-after within-allele mean {after_within.mean:.4f}; the measured "
          "value is 0.5120 -- near chance, which is exactly what makes the pooled "
          "figure beside it convincing")

    real = HydropathyLengthAllele(eval_ids=EVAL_IDS).fit(FIT_ROWS)
    check(pooled(real.score(EVAL_ROWS)) > after.mean - 0.02,
          "the unpermuted model should be at least as good as its wrongly-permuted "
          "sentinel; if it is not, the sentinel is not a sentinel")


def test_sentinel_default_order_is_the_safe_one() -> None:
    s = ShuffledSentinel(HydropathyLengthAllele, eval_ids=EVAL_IDS, seed=0)
    check(s.order == "before",
          f"default order is {s.order!r}; the default must be the reportable one")
    expect_raises(BaselineError, ShuffledSentinel, HydropathyLengthAllele,
                  eval_ids=EVAL_IDS, order="afterwards")
    # A sentinel over a label-free baseline is the baseline with a misleading name.
    expect_raises(BaselineError, ShuffledSentinel, Hydropathy, eval_ids=EVAL_IDS)


def test_sentinel_refuses_a_single_draw() -> None:
    expect_raises(BaselineError, sentinel_distribution, HydropathyLengthAllele,
                  FIT_ROWS, EVAL_ROWS, pooled, eval_ids=EVAL_IDS, seeds=[0])


def test_sentinel_does_not_mutate_the_fitting_pool() -> None:
    labels = [r["label"] for r in FIT_ROWS]
    ShuffledSentinel(HydropathyLengthAllele, eval_ids=EVAL_IDS, seed=3).fit(FIT_ROWS)
    check(labels == [r["label"] for r in FIT_ROWS],
          "the sentinel permuted the caller's rows in place; every later column in the "
          "table would have been fitted on shuffled labels")


def test_fit_raises_on_an_evaluation_record_id() -> None:
    """A baseline fitted on the evaluation labels is indistinguishable from a good one
    until someone checks, so fit() checks, by record_id."""
    contaminated = FIT_ROWS[:50] + EVAL_ROWS[:1]
    for factory in (AlleleOnly, Hydropathy, Length, HydropathyLengthAllele):
        b = factory(eval_ids=EVAL_IDS)
        expect_raises(EvaluationLeakError, b.fit, contaminated)
    pwm = PWM(eval_ids=EVAL_IDS, min_rows=1)
    expect_raises(EvaluationLeakError, pwm.fit, contaminated)
    sent = ShuffledSentinel(HydropathyLengthAllele, eval_ids=EVAL_IDS)
    expect_raises(EvaluationLeakError, sent.fit, contaminated)

    # The clean pool must still fit, or the guard proves nothing.
    AlleleOnly(eval_ids=EVAL_IDS).fit(FIT_ROWS)

    # The guard is on record_id: a re-spelled evaluation row keeps its id.
    respelled = dict(EVAL_ROWS[0])
    respelled["mhc_allele"] = "A0201"
    respelled["peptide"] = respelled["peptide"].lower()
    expect_raises(EvaluationLeakError, AlleleOnly(eval_ids=EVAL_IDS).fit,
                  FIT_ROWS[:20] + [respelled])


def test_eval_ids_must_be_declared_explicitly() -> None:
    expect_raises(BaselineError, AlleleOnly, eval_ids=None)
    expect_raises(TypeError, AlleleOnly)
    # Declaring "no held-out set" is allowed, but it has to be typed.
    AlleleOnly(eval_ids=()).fit(FIT_ROWS)


def test_allele_only_is_exactly_half_within_allele() -> None:
    """Pooled ~0.70 read as capability: AlleleOnly has no peptide input and is EXACTLY
    0.5000 within-allele, by construction."""
    b =AlleleOnly(eval_ids=EVAL_IDS).fit(FIT_ROWS)
    s = b.score(EVAL_ROWS)
    w = stratified_auroc(Y, s, G)
    check(abs(w.value - 0.5) < 1e-12,
          f"within-allele {w.value!r}; it must be EXACTLY 0.5 -- the score is constant "
          "inside a stratum, so any other value means the strata are not the alleles")
    for g, e in w.per_stratum.items():
        if e["admissible"]:
            check(abs(e["auroc"] - 0.5) < 1e-12,
                  f"stratum {g} scored {e['auroc']!r}, expected exactly 0.5")

    check(len(set(s)) == len(ALLELES),
          f"expected one distinct score per allele, got {len(set(s))}")
    check(pooled(s) > 0.60,
          f"pooled {pooled(s):.4f}; a lookup with no peptide input should reach far "
          "above chance on this fixture, which is the entire reason it is a row")

    scrambled = [dict(r, peptide="AAAAAAAAA") for r in EVAL_ROWS]
    check(b.score(scrambled) == s,
          "AlleleOnly changed its answer when the peptides changed; it has acquired a "
          "peptide input and is no longer the shortcut it is shipped as")


def test_unseen_allele_falls_back_to_the_global_rate_and_is_counted() -> None:
    b = AlleleOnly(eval_ids=EVAL_IDS).fit(FIT_ROWS)
    check("HLA-A*24:02" not in ALLELES, "pick an allele the fixture does not contain")
    unseen = [dict(EVAL_ROWS[0], mhc_allele="HLA-A*24:02")]
    got = b.score(unseen)[0]
    check(abs(got - b.encoder.global_rate) < 1e-12,
          f"unseen allele scored {got!r}, expected the global rate "
          f"{b.encoder.global_rate!r} -- never 0.5, which is the value a missed join "
          "produces and the value that got caught only because it looks impossible")
    check(b.encoder.n_fallback == 1,
          f"fallback counter is {b.encoder.n_fallback}, expected 1; an uncounted "
          "fallback is a silent join failure")


def test_hydropathy_is_deterministic_and_label_free() -> None:
    """Hydropathy reads no labels, so fitting must be incapable of changing it."""
    unfitted = Hydropathy(eval_ids=EVAL_IDS).score(EVAL_ROWS)
    fitted = Hydropathy(eval_ids=EVAL_IDS).fit(FIT_ROWS).score(EVAL_ROWS)
    check(unfitted == fitted, "fitting changed Hydropathy's scores; it reads no labels")

    flipped = [dict(r, label=1 - r["label"]) for r in FIT_ROWS]
    check(Hydropathy(eval_ids=EVAL_IDS).fit(flipped).score(EVAL_ROWS) == fitted,
          "inverting every fitting label changed Hydropathy's scores")

    other_pool = make_rows(300, seed=4242, prefix="other")
    check(Hydropathy(eval_ids=EVAL_IDS).fit(other_pool).score(EVAL_ROWS) == fitted,
          "a different fitting pool changed Hydropathy's scores")

    check(fitted == [mean_kd(r["peptide"]) for r in EVAL_ROWS],
          "Hydropathy is not mean Kyte-Doolittle hydropathy")
    check(len(set(fitted)) > 1, "degenerate fixture: hydropathy is constant")

    # Length is the other unfitted row, and it is the one that is anti-predictive.
    ls = Length(eval_ids=EVAL_IDS).fit(FIT_ROWS).score(EVAL_ROWS)
    check(ls == [float(len(r["peptide"])) for r in EVAL_ROWS], "Length is not length")


KD_ROWS = make_rows(1200, seed=5, prefix="kd", kd_effect=0.18)
KD_Y = [r["label"] for r in KD_ROWS]
KD_G = [r["mhc_allele"] for r in KD_ROWS]


def test_residualising_hydropathy_against_itself_drives_it_to_chance() -> None:
    """A score that IS the basis must residualise to chance; this self-check is the one
    that FAILS on the neoantigen set at 0.5748, n=833, so no figure there is reportable."""
    raw = [mean_kd(r["peptide"]) for r in KD_ROWS]
    raw_w = stratified_auroc(KD_Y, raw, KD_G).value
    check(raw_w > 0.68,
          f"fixture is not hydropathy-driven: raw within-allele {raw_w:.4f}; the test "
          "proves nothing unless the signal is there to remove")

    res = residualise(raw, KD_ROWS)
    res_w = stratified_auroc(KD_Y, res, KD_G).value
    check(abs(res_w - 0.5) < 0.06,
          f"residualised within-allele {res_w:.4f}; a score that IS the basis must go "
          "to chance. This self-check is the one that FAILS on the neoantigen set "
          "(0.5748 at n=833), which is why no residualised neoantigen figure is "
          "reportable -- and why it is asserted here rather than assumed")
    check(abs(res_w - 0.5) < abs(raw_w - 0.5) - 0.15,
          f"residualisation removed too little: {raw_w:.4f} -> {res_w:.4f}")
    # Not exactly 0.5: a quadratic basis against a sigmoid within-allele rank leaves a
    # cubic remnant, ~+0.03 at n=1,200 here and +0.0748 at n=833 on the neoantigen set.

    # OLS orthogonality holds to machine precision, or the solve is wrong.
    kd = [mean_kd(r["peptide"]) for r in KD_ROWS]
    length = [float(len(r["peptide"])) for r in KD_ROWS]
    for nm, col in (("1", [1.0] * len(res)), ("kd", kd),
                    ("kd2", [v * v for v in kd]), ("length", length)):
        dot = sum(a * b for a, b in zip(res, col))
        check(abs(dot) < 1e-6 * max(1.0, sum(abs(v) for v in col)),
              f"residual is not orthogonal to {nm}: dot product {dot:.3e}")


def test_residualisation_leaves_a_column_the_basis_cannot_explain() -> None:
    """An artifact shrinks, a leak does not move: ESM-2 0.5540 -> 0.5336 and PRIME-2.0
    0.5348 -> 0.5183, while the leaked TITANiAN column went 0.7040 -> 0.7113."""
    rng = random.Random(7)
    leaked = [r["label"] + 0.30 * rng.random() for r in KD_ROWS]
    before = stratified_auroc(KD_Y, leaked, KD_G).value
    after = stratified_auroc(KD_Y, residualise(leaked, KD_ROWS), KD_G).value
    check(before > 0.95, f"stand-in leaked column is only {before:.4f}; rebuild it")
    check(after > before - 0.05,
          f"a column the basis cannot explain lost {before - after:+.4f} "
          f"({before:.4f} -> {after:.4f}); residualisation is removing more than "
          "hydropathy and length, so the artifact/leak contrast is not measurable")


def test_residualise_raises_on_a_design_it_cannot_fit() -> None:
    nine_mers = make_rows(120, seed=8, prefix="nine", lengths=(9,))
    raw = [mean_kd(r["peptide"]) for r in nine_mers]
    # All 9-mers, so `length` is collinear with the intercept; dropping it quietly would
    # claim length was removed when it was not.
    expect_raises(SingularDesignError, residualise, raw, nine_mers)
    # Naming the admissible basis is allowed; having it chosen for you is not.
    out = residualise(raw, nine_mers, predictors=("kd", "kd2"))
    check(len(out) == len(nine_mers), "residualise returned the wrong length")
    expect_raises(BaselineError, residualise, raw, nine_mers, predictors=())
    expect_raises(BaselineError, residualise, raw, nine_mers, predictors=("charge",))


def test_within_allele_rank_puts_a_singleton_at_exactly_half() -> None:
    rows = [dict(EVAL_ROWS[0], record_id="x", mhc_allele="HLA-A*01:01")]
    check(within_allele_rank([3.0], rows) == [0.5],
          "a singleton stratum must rank to exactly 0.5")
    r = within_allele_rank([1.0, 2.0, 3.0, 4.0],
                           [dict(EVAL_ROWS[0], mhc_allele="HLA-A*02:01")] * 4)
    check(r == [0.125, 0.375, 0.625, 0.875], f"mid-rank mapping is wrong: {r}")


def test_pwm_has_760_parameters_and_refuses_unfitted_lengths() -> None:
    p = PWM(eval_ids=EVAL_IDS, min_rows=1).fit(FIT_ROWS)
    check(sorted(p.matrices) == [8, 9, 10, 11],
          f"fitted lengths {sorted(p.matrices)}, expected 8-11")
    check(p.n_parameters == 760,
          f"{p.n_parameters} parameters, expected 760 = 20 * (8+9+10+11)")
    s = p.score(EVAL_ROWS)
    check(len(set(s)) > 1, "PWM produced a constant column")

    long_pep = [dict(EVAL_ROWS[0], peptide="ACDEFGHIKLMNA")]  # never fitted
    expect_raises(BaselineError, p.score, long_pep)
    lenient = PWM(eval_ids=EVAL_IDS, min_rows=1, strict=False).fit(FIT_ROWS)
    check(lenient.score(long_pep)[0] != lenient.score(long_pep)[0],
          "non-strict PWM should return NaN, which is not rankable, rather than 0.0")
    expect_raises(BaselineError, PWM(eval_ids=EVAL_IDS, min_rows=10**6).fit, FIT_ROWS)


def test_four_input_model_beats_its_own_one_input_parts() -> None:
    m = HydropathyLengthAllele(eval_ids=EVAL_IDS).fit(FIT_ROWS)
    check(m.coefficients is not None and len(m.coefficients) == 4,
          f"expected 4 OLS coefficients, got {m.coefficients}")
    s = m.score(EVAL_ROWS)
    allele_only = AlleleOnly(eval_ids=EVAL_IDS).fit(FIT_ROWS).score(EVAL_ROWS)
    check(pooled(s) >= pooled(allele_only) - 0.03,
          f"four-input pooled {pooled(s):.4f} against allele-only "
          f"{pooled(allele_only):.4f}; the measured pair is 0.7084 against ~0.70")
    check(stratified_auroc(Y, s, G).value != 0.5,
          "the four-input model is constant within allele; its peptide terms are dead")
    expect_raises(BaselineError, HydropathyLengthAllele(eval_ids=EVAL_IDS).score,
                  EVAL_ROWS)


def test_ols_recovers_a_known_fit_and_refuses_a_singular_one() -> None:
    X = [[1.0, x, x * x] for x in (-2, -1, 0, 1, 2, 3)]
    y = [2.0 + 3.0 * r[1] - 0.5 * r[2] for r in X]
    b = ols(X, y)
    for got, want in zip(b, (2.0, 3.0, -0.5)):
        check(abs(got - want) < 1e-9, f"OLS recovered {b}, expected [2, 3, -0.5]")
    expect_raises(SingularDesignError, ols, [[1.0, 1.0], [1.0, 1.0]], [0.0, 1.0])
    expect_raises(SingularDesignError, ols, [[1.0, 2.0]], [1.0])


def test_rows_with_bad_fields_raise_instead_of_scoring() -> None:
    b = AlleleOnly(eval_ids=EVAL_IDS).fit(FIT_ROWS)
    expect_raises(BaselineError, b.score, [{"peptide": "SIINFEKL"}])
    expect_raises(BaselineError, AlleleOnly(eval_ids=EVAL_IDS).fit,
                  [dict(FIT_ROWS[0], label=2)])
    expect_raises(BaselineError, AlleleOnly(eval_ids=EVAL_IDS).fit,
                  [dict(FIT_ROWS[0], record_id="")])
    expect_raises(BaselineError, AlleleOnly(eval_ids=EVAL_IDS).fit, [])



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
