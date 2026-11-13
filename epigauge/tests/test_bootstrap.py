"""Regression fixtures for the resampling layer.

Governing rule in `test_canon_and_metrics.py`. Each fixture either reproduces a wrong
measurement or computes a quantity by hand. Run under pytest, or standalone with
`python3 -I epigauge/tests/test_bootstrap.py`.
"""

from __future__ import annotations

import random
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from epigauge.bootstrap import (  # noqa: E402
    ClusterMap,
    Delta,
    Estimate,
    MixedBasisError,
    PreRegistrationUnderpowered,
    SeedRecipe,
    Unmeasurable,
    UnmeasuredBasis,
    assign_clusters,
    draws,
    format_table,
    interval,
    mde,
    paired,
    study_level_basis,
    z_sum,
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
    FAILURES.append(f"expected {exc.__name__}, nothing raised: {getattr(fn,'__name__',fn)}")
    return False


PROT_A = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"
)
PROT_B = (
    "MGSSVPRGSHMASMTGGQQMGRGSEFEKLAVYNTWQRDVAIAGCSVNPYTQLKWMFDTAKGLSPEQCIDVYNRVAEG"
)
PROT_C = "MDPNTVSSFQVDCFLWHVRKQVADQELGDAPFLDRLRRDQRQISFVKSEWTNHIFYNGICKEYPDLV"
SHARED_9 = "QRQISFVKS"  # occurs in PROT_A and in PROT_C -> must merge the two proteins


def _tmpdir() -> Path:
    return Path(tempfile.mkdtemp(prefix="epigauge_boot_"))


def _write_fasta(dirpath: Path, name: str, entries) -> Path:
    p = dirpath / name
    with open(p, "w", encoding="utf-8") as fh:
        for acc, seq in entries:
            fh.write(f">sp|{acc}|{acc}_TEST test protein\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    return p


def mean_score(labels, scores):
    """A mean has a closed-form bootstrap variance, which is what makes the design
    effect below hand-computable."""
    return sum(scores) / len(scores)


RECIPE = SeedRecipe(version="epigauge-1", set_id="infectious", comparison_id="arm-vs-prime")


def test_assigned_fraction_reproduces_the_one_proteome_collapse() -> None:
    """One proteome assigns 96 of 2,093 rows (4.6%) and was reported as a clustered
    design; all nine, with an 8-mer window, reach 85.9% over 323 proteins."""
    d = _tmpdir()
    human = _write_fasta(d, "human_reviewed.fasta", [("P11111", PROT_A)])
    pathogen = _write_fasta(d, "virus_reviewed.fasta", [("P22222", PROT_B), ("P33333", PROT_C)])

    from_a = [PROT_A[0:9], PROT_A[20:30]]
    from_b = [PROT_B[i : i + 9] for i in range(0, 45, 5)]
    from_c = [PROT_C[i : i + 10] for i in range(0, 40, 8)]
    peptides = from_a + from_b + from_c

    one = assign_clusters(peptides, [human])
    all_nine = assign_clusters(peptides, [human, pathogen])

    check(
        one.assigned_fraction < 0.25,
        f"one-proteome assignment should collapse, got {one.assigned_fraction:.3f}",
    )
    check(
        one.provisional,
        "an assignment below 0.80 must be flagged provisional; that is the only thing "
        "separating a 4.6% run from an 85.9% run on the page",
    )
    check(
        all_nine.assigned_fraction == 1.0,
        f"all-proteome assignment should reach every row, got {all_nine.assigned_fraction:.3f}",
    )
    check(not all_nine.provisional, "a fully assigned map must not be provisional")
    # Cluster count alone cannot tell the two maps apart, which is why the fraction is
    # carried on the map and printed into the basis string.
    check(
        one.n_clusters > 0 and all_nine.n_clusters > 0,
        "both maps must produce clusters at all",
    )


def test_eight_mers_are_matched_whole_not_dropped() -> None:
    d = _tmpdir()
    fa = _write_fasta(d, "p.fasta", [("P11111", PROT_A)])
    eight = PROT_A[10:18]
    check(len(eight) == 8 and eight in PROT_A, "fixture 8-mer must be a real substring")

    cmap = assign_clusters([eight, PROT_A[0:9]], [fa], k_window=9, short_window=8)
    check(
        cmap.assigned_fraction == 1.0,
        f"the 8-mer must assign under the 8-mer window; got {cmap.assigned_fraction}",
    )
    check(
        cmap.keys[cmap.cluster_of_row[0]].startswith("PROTEIN:"),
        "the 8-mer row landed in a SINGLETON despite occurring in the proteome",
    )
    # A 7-mer has no window and must raise: a skipped row is indistinguishable from an
    # unmatched one in the assigned fraction.
    expect_raises(ValueError, assign_clusters, ["MKTAYIA"], [fa])


def test_multi_protein_hits_merge_instead_of_double_counting() -> None:
    d = _tmpdir()
    fa = _write_fasta(d, "p.fasta", [("P11111", PROT_A), ("P33333", PROT_C)])
    check(SHARED_9 in PROT_A and SHARED_9 in PROT_C, "fixture peptide must hit both proteins")

    only_a = PROT_A[0:9]
    only_c = PROT_C[0:9]
    check(only_a not in PROT_C and only_c not in PROT_A, "fixture anchors must be unique")

    cmap = assign_clusters([only_a, SHARED_9, only_c], [fa])
    c = cmap.cluster_of_row
    check(
        c[0] == c[1] == c[2],
        f"an ambiguous peptide must MERGE its proteins by union-find; clusters were {c}. "
        "Assigning it to one winner lets a replicate draw both proteins and count the "
        "shared row twice.",
    )
    check(cmap.n_clusters == 1, f"expected one merged cluster, got {cmap.n_clusters}")

    # Order independence: the reversed FASTA gives the same partition.
    fa_rev = _write_fasta(d, "p_rev.fasta", [("P33333", PROT_C), ("P11111", PROT_A)])
    rev = assign_clusters([only_a, SHARED_9, only_c], [fa_rev])
    check(rev.keys == cmap.keys, f"merge target depends on file order: {rev.keys} vs {cmap.keys}")


def test_unassigned_rows_become_singletons_not_one_giant_cluster() -> None:
    """Pooling the unassigned rows at 85.9% of 2,093 invents a single 295-row cluster --
    the largest in the design, and an artifact of which FASTA files happen to ship."""
    d = _tmpdir()
    fa = _write_fasta(d, "p.fasta", [("P11111", PROT_A)])

    assigned = [PROT_A[0:9], PROT_A[5:14], PROT_A[30:39]]
    # Twelve rows no proteome can place; pooled they would be the largest cluster here.
    unassigned = [
        "WWWWWWWWW", "YYYYYYYYY", "CCCCCCCCC", "HHHHHHHHH",
        "WWWWWWWWC", "YYYYYYYYC", "CCCCCCCCH", "HHHHHHHHW",
        "WWWWWWWCC", "YYYYYYYCC", "CCCCCCCHH", "HHHHHHHWW",
    ]
    for pep in unassigned:
        check(pep not in PROT_A, f"fixture peptide {pep} must not occur in the proteome")

    peptides = assigned + unassigned
    cmap = assign_clusters(peptides, [fa])

    check(
        cmap.n_rows == len(peptides),
        f"no row may be dropped: {cmap.n_rows} of {len(peptides)} survived",
    )
    check(
        cmap.n_assigned == len(assigned),
        f"expected {len(assigned)} assigned rows, got {cmap.n_assigned}",
    )

    unassigned_clusters = {cmap.cluster_of_row[i] for i in range(len(assigned), len(peptides))}
    check(
        len(unassigned_clusters) == len(unassigned),
        f"the {len(unassigned)} unmatched rows collapsed into {len(unassigned_clusters)} "
        "cluster(s); pooling them fabricates the largest cluster in the design",
    )
    sizes = [len(cmap.members[c]) for c in unassigned_clusters]
    check(set(sizes) == {1}, f"unassigned clusters must be singletons, sizes were {sorted(sizes)}")

    biggest = max(len(m) for m in cmap.members)
    check(
        biggest == len(assigned),
        f"the largest cluster has {biggest} rows; the only legitimate large cluster here "
        "is the real protein, not an unassigned pool",
    )
    check(
        cmap.n_singleton_clusters == len(unassigned),
        f"expected {len(unassigned)} singleton clusters, got {cmap.n_singleton_clusters}",
    )
    check(
        abs(cmap.assigned_fraction - len(assigned) / len(peptides)) < 1e-12,
        f"assigned_fraction {cmap.assigned_fraction} disagrees with the row counts",
    )


def test_seed_recipe_is_deterministic_and_order_independent() -> None:
    """A sequentially-seeded bootstrap makes replicate b depend on how many comparisons
    ran before it, while every arm still looks internally consistent."""
    a = draws(25, B=64, seed_recipe=RECIPE)
    b = draws(25, B=64, seed_recipe=RECIPE)
    check(a.matrix == b.matrix, "the same recipe must give identical draws")

    # Replicates do not depend on how many were asked for, so a table recomputed in a
    # different order, or with a cell omitted, reproduces.
    short = draws(25, B=8, seed_recipe=RECIPE)
    check(
        short.matrix == a.matrix[:8],
        "replicate b changed when B changed; the seeding is order-dependent",
    )

    other = draws(25, B=64, seed_recipe=SeedRecipe("epigauge-1", "infectious", "arm-vs-mhcflurry"))
    check(
        other.matrix != a.matrix,
        "two comparisons on one set must draw independently, or their deltas share noise",
    )
    check(
        other.recipe.digest() == RECIPE.digest(),
        "comparison_id must not change the comparability digest, or every multi-arm table "
        "becomes a mixed-basis error for a non-methodological reason",
    )

    newer = draws(25, B=8, seed_recipe=SeedRecipe("epigauge-2", "infectious", "arm-vs-prime"))
    check(newer.recipe.digest() != RECIPE.digest(), "a version bump must change the digest")

    # The documented seed, recomputed independently of the implementation.
    import hashlib

    want = int.from_bytes(
        hashlib.sha256(b"epigauge-1|infectious|arm-vs-prime|7").digest()[:8], "big"
    )
    check(RECIPE.seed(7) == want, f"seed(7) = {RECIPE.seed(7)}, expected {want}")


def test_two_stage_draw_leaves_n_ragged_and_does_not_renormalise() -> None:
    keys = ["BIG"] * 30 + [f"SINGLETON:{i}" for i in range(10)]
    cmap = ClusterMap.from_keys(keys)
    dm = draws(cmap.n_clusters, B=200, seed_recipe=RECIPE)
    counts = {sum(len(cmap.members[c]) for c in drawn) for drawn in dm.matrix}
    check(
        len(counts) > 1,
        "every replicate had the same row count; the two-stage draw was renormalised, "
        "which restores exactly the row-level variance clustering was adopted to avoid",
    )
    check(
        min(counts) < cmap.n_rows < max(counts),
        f"replicate sizes {min(counts)}-{max(counts)} should straddle n={cmap.n_rows}",
    )


def test_design_effect_matches_the_hand_computed_cluster_size() -> None:
    """K clusters of m identical rows: cluster Var is m times row Var, so widths differ
    by sqrt(m) and DE = m exactly in the large-B limit. Here K=20, m=5."""
    K, m = 20, 5
    values = [k / float(K) for k in range(K)]
    labels, a_scores, b_scores, keys = [], [], [], []
    for k in range(K):
        for _ in range(m):
            labels.append(1)
            a_scores.append(values[k])
            b_scores.append(0.0)
            keys.append(f"PROTEIN:P{k}")
    groups = ["g"] * len(labels)
    cmap = ClusterMap.from_keys(keys)
    check(cmap.n_clusters == K, f"expected {K} clusters, got {cmap.n_clusters}")

    dm = draws(K, B=4000, seed_recipe=RECIPE)
    dl = paired(mean_score, labels, a_scores, b_scores, groups, cmap, dm, study_ids=["S"] * len(labels))

    check(
        abs(dl.d - sum(values) / K) < 1e-12,
        f"point delta {dl.d} should be the mean of the cluster values {sum(values)/K}",
    )
    check(
        abs(dl.design_effect - m) < 1.5,
        f"design effect {dl.design_effect:.3f} is not within 1.5 of the hand-computed "
        f"value m = {m} (cluster width {dl.cluster_width:.5f}, row width "
        f"{dl.row_width:.5f}, ratio {dl.cluster_width/dl.row_width:.4f} vs sqrt(m) "
        f"{m ** 0.5:.4f})",
    )
    check(
        abs(dl.n_eff - len(labels) / dl.design_effect) < 1e-9,
        "n_eff must be exactly n / DE",
    )
    check(
        dl.n_eff < len(labels),
        f"n_eff {dl.n_eff:.1f} should fall below n = {len(labels)} when DE > 1",
    )
    check(dl.basis.startswith("cluster/"), f"delta basis {dl.basis!r} must be clustered")
    check(dl.row_basis.startswith("row/"), f"row basis {dl.row_basis!r} must be row-level")
    check(
        isinstance(dl.study_level, Unmeasurable),
        "a constant study_id must come back as the sentinel, not a basis string",
    )


def test_paired_uses_one_shared_index_so_common_noise_cancels() -> None:
    """Arm A is arm B plus 0.05: shared-index pairing must return a zero-width interval,
    which two independently drawn marginal intervals would not."""
    rnd = random.Random(11)
    n = 60
    labels = [1] * n
    b_scores = [rnd.random() for _ in range(n)]
    a_scores = [x + 0.05 for x in b_scores]
    groups = ["g"] * n
    cmap = ClusterMap.from_keys([f"PROTEIN:P{i//3}" for i in range(n)])
    dm = draws(cmap.n_clusters, B=300, seed_recipe=RECIPE)
    dl = paired(mean_score, labels, a_scores, b_scores, groups, cmap, dm)
    check(abs(dl.d - 0.05) < 1e-12, f"point delta {dl.d}, expected 0.05")
    check(
        dl.cluster_width < 1e-12,
        f"shared-index pairing left width {dl.cluster_width:.3e}; the arms are not being "
        "scored on the same resampled rows",
    )
    check(
        abs(dl.p - 1.0 / (dm.B + 1)) < 1e-12,
        f"p = {dl.p}, expected the honest floor 1/(B+1) = {1.0/(dm.B+1)}; a bare "
        "#{delta<=0}/B prints 0.0000, which B replicates cannot support",
    )


def test_z_sum_is_the_filed_constant() -> None:
    """"+0.0246 over PRIME-2.0, p = 0.078" recomputes paired and clustered to +0.0187
    [-0.0111, +0.0466], p = 0.215, SE 0.0144 -- below even MDE_row of 0.0403."""
    zs = z_sum(alpha=0.025, power=0.80)
    check(abs(zs - 2.801585) < 5e-7, f"z_sum = {zs:.6f}, expected 2.801585")
    mde_row = zs * 0.0144
    check(
        abs(mde_row - 0.0403) < 5e-4,
        f"MDE_row from the record's SE 0.0144 is {mde_row:.4f}, expected about 0.0403",
    )
    check(
        mde_row > 0.0187,
        f"the record's +0.0187 margin must sit BELOW MDE_row {mde_row:.4f}, or this "
        "fixture no longer reproduces the failure it exists for",
    )


def _prime_like_delta(B: int = 1200) -> Delta:
    """K=100 protein clusters of 4 rows, calibrated to the record: margin near +0.0187,
    clustered SE near the measured 0.0144."""
    rnd = random.Random(2026)
    K, m = 100, 4
    labels, a_scores, b_scores, keys = [], [], [], []
    for k in range(K):
        per_cluster = rnd.gauss(0.0187, 0.144)
        for _ in range(m):
            labels.append(1)
            a_scores.append(per_cluster)
            b_scores.append(0.0)
            keys.append(f"PROTEIN:P{k}")
    cmap = ClusterMap.from_keys(keys)
    dm = draws(K, B=B, seed_recipe=RECIPE)
    return paired(
        mean_score,
        labels,
        a_scores,
        b_scores,
        ["g"] * len(labels),
        cmap,
        dm,
        study_ids=["32060144"] * len(labels),
    )


def test_underpowered_pre_registration_raises_rather_than_warns() -> None:
    dl = _prime_like_delta()
    check(
        0.008 < dl.se < 0.025,
        f"fixture SE {dl.se:.5f} should sit near the record's 0.0144, or the acceptance "
        "test below is not testing what it claims",
    )

    report = mde(dl, alpha=0.025, power=0.80)
    check(
        report.mde > 0.0187,
        f"fixture MDE {report.mde:.4f} must exceed the +0.0187 margin for the refusal to "
        "be the correct outcome",
    )
    check(abs(report.z_sum - 2.801585) < 5e-7, f"report z_sum {report.z_sum}")
    check(
        abs(report.mde - report.z_sum * dl.se) < 1e-12,
        "MDE must be exactly (z_{1-alpha} + z_{1-power}) * se_delta",
    )

    # The acceptance test from the record: a declared threshold of 0.0187 is REFUSED.
    expect_raises(
        PreRegistrationUnderpowered,
        mde,
        dl,
        alpha=0.025,
        power=0.80,
        declared_threshold=0.0187,
    )
    # A raise, not a flag on a returned object: there is nothing to ignore.
    try:
        mde(dl, declared_threshold=0.0187)
    except PreRegistrationUnderpowered as e:
        msg = str(e)
        check("not pre-registerable" in msg, f"the refusal must say why: {msg!r}")
        check("n_ef" in msg and "design effect" in msg,
              f"the refusal must quote the effective n and design effect: {msg!r}")

    # A resolvable threshold passes, so the gate is not a blanket no.
    ok = mde(dl, declared_threshold=0.25)
    check(ok.declared_threshold == 0.25, "a resolvable threshold must be accepted")
    check("clears the MDE" in ok.verdict, f"unexpected verdict {ok.verdict!r}")


def test_mde_refuses_when_cluster_draws_do_not_exist() -> None:
    dl = _prime_like_delta(B=200)

    orphan = Delta(**{**dl.__dict__, "draws": None})
    expect_raises(UnmeasuredBasis, mde, orphan, declared_threshold=0.5)

    row_only = Delta(
        **{
            **dl.__dict__,
            "draws": draws(8, B=4, seed_recipe=RECIPE, unit="row"),
        }
    )
    expect_raises(UnmeasuredBasis, mde, row_only, declared_threshold=0.5)


def _two_arm_estimates(B: int = 400):
    rnd = random.Random(5)
    n = 80
    labels = [1] * n
    groups = ["g"] * n
    keys = [f"PROTEIN:P{i // 4}" for i in range(n)]
    cmap = ClusterMap.from_keys(keys)
    rowmap = ClusterMap.rows_as_clusters(n)
    arm_a = [rnd.random() for _ in range(n)]
    arm_b = [rnd.random() for _ in range(n)]

    cl_draws = draws(cmap.n_clusters, B=B, seed_recipe=RECIPE)
    rw_draws = draws(n, B=B, seed_recipe=RECIPE.derive("row"), unit="row")
    study = ["32060144"] * n

    clustered = interval(mean_score, labels, arm_a, groups, cmap, cl_draws, study_ids=study)
    clustered_b = interval(mean_score, labels, arm_b, groups, cmap, cl_draws, study_ids=study)
    row_level = interval(mean_score, labels, arm_b, groups, rowmap, rw_draws, study_ids=study)
    return clustered, clustered_b, row_level


def test_mixed_basis_table_raises_instead_of_printing() -> None:
    """A row-level interval beside a clustered one widens by a model-dependent 1.00-1.08x,
    so the table does not look wrong -- it just reorders the ranking."""
    clustered, clustered_b, row_level = _two_arm_estimates()

    check(
        clustered.basis != row_level.basis,
        "the fixture no longer reproduces: the two bases must differ",
    )
    check(clustered.basis.startswith("cluster/"), f"basis {clustered.basis!r}")
    check(row_level.basis.startswith("row/"), f"basis {row_level.basis!r}")

    expect_raises(
        MixedBasisError,
        format_table,
        [("ours", clustered), ("PRIME-2.0", row_level)],
    )
    try:
        format_table({"ours": clustered, "PRIME-2.0": row_level})
    except MixedBasisError as e:
        check("reorders" in str(e), f"the refusal must say what goes wrong: {str(e)!r}")
        check("PRIME-2.0" in str(e) and "ours" in str(e),
              f"the refusal must name both cells: {str(e)!r}")

    # Same unit, different B is ALSO a mixed basis: B moves the width too.
    rnd = random.Random(5)
    other_B = _two_arm_estimates(B=200)[0]
    check(other_B.basis != clustered.basis, "B must appear in the basis string")
    expect_raises(MixedBasisError, format_table, [("a", clustered), ("b", other_B)])

    # Happy path: one basis, and every row prints its study-level sentinel.
    out = format_table({"ours": clustered, "PRIME-2.0": clustered_b}, title="infectious")
    check("unmeasurable" in out, f"the study sentinel must be printed, not omitted:\n{out}")
    check(clustered.basis in out, f"the basis line must be printed:\n{out}")
    check(out.count("\n") >= 3, f"expected a title, a basis line and two rows:\n{out}")


def test_study_level_is_a_required_field_carrying_a_sentinel() -> None:
    """study_id is constant across all 2,093 infectious rows, so study-level clustering
    is unmeasurable -- a value, not a missing key."""
    constant = study_level_basis(["32060144"] * 2093)
    check(isinstance(constant, Unmeasurable), f"expected a sentinel, got {constant!r}")
    check(str(constant) == "unmeasurable", "the sentinel must render as 'unmeasurable'")
    check("constant" in constant.reason, f"the sentinel must carry a reason: {constant.reason!r}")
    check(not constant, "the sentinel must be falsy so it cannot be mistaken for a basis")

    absent = study_level_basis(None)
    check(isinstance(absent, Unmeasurable), "an absent column is also unmeasurable, not None")

    real = study_level_basis(["a"] * 50 + ["b"] * 50 + ["c"] * 50)
    check(isinstance(real, str) and "3-clusters" in real, f"measurable case gave {real!r}")

    # The field cannot be omitted: the Estimate refuses to exist without it.
    kw = dict(
        value=0.55, lo=0.53, hi=0.58, mc_se=(0.001, 0.001), basis="cluster/10/fresh/ab/0.859",
        n_rows=2093, n_clusters=323, assigned_fraction=0.859,
        n_replicates_used=10, n_replicates_failed=0, provisional=False,
    )
    expect_raises(UnmeasuredBasis, Estimate, study_level=None, **kw)
    ok = Estimate(study_level=constant, **kw)
    check(isinstance(ok.study_level, Unmeasurable), "the sentinel must survive construction")


def test_interval_carries_mc_se_and_counts_refused_replicates() -> None:
    clustered, _, _ = _two_arm_estimates(B=400)
    check(len(clustered.mc_se) == 2, "mc_se must cover BOTH bounds")
    check(
        all(x >= 0.0 for x in clustered.mc_se),
        f"mc_se must be non-negative, got {clustered.mc_se}",
    )
    check(
        max(clustered.mc_se) < clustered.width,
        f"mc_se {clustered.mc_se} exceeds the interval width {clustered.width}; the "
        "interval would not be quotable at all",
    )
    check(clustered.n_replicates_failed == 0, "no replicate should fail on a mean metric")
    check(
        clustered.n_replicates_used + clustered.n_replicates_failed == 400,
        "used + failed must account for every replicate",
    )

    # A metric that refuses most resamples must raise, not quietly estimate on the rest.
    n = 40

    def picky(labels, scores):
        # Admissible on the full data, refused on almost every resample.
        if len(set(scores)) < n - 2:
            raise ValueError("refused")
        return sum(scores) / len(scores)

    cmap = ClusterMap.from_keys([f"P{i // 4}" for i in range(n)])
    dm = draws(cmap.n_clusters, B=50, seed_recipe=RECIPE)
    expect_raises(
        UnmeasuredBasis,
        interval,
        picky, [1] * n, [i / n for i in range(n)], ["g"] * n, cmap, dm,
    )


def test_draw_matrix_must_describe_the_cluster_map_it_is_used_with() -> None:
    n = 20
    cmap = ClusterMap.from_keys([f"P{i // 2}" for i in range(n)])  # 10 clusters
    wrong = draws(7, B=10, seed_recipe=RECIPE)
    expect_raises(
        UnmeasuredBasis,
        interval,
        mean_score, [1] * n, [0.5] * n, ["g"] * n, cmap, wrong,
    )



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
