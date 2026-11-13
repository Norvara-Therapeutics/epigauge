# Metrics

## `auroc(labels, scores)`

`AUROC = (sum of mid-ranks of positives - n+ (n+ + 1) / 2) / (n+ n-)`

Ties take mid-ranks (`0.5*(i+j)+1` over the tied block), so a constant score is exactly 0.5000 and partial ties do not depend on input order. Raises `DegenerateGroupingError` when `n+ == 0` or `n- == 0`: undefined, not 0.5 and not the prevalence.

## `stratified_auroc(labels, scores, groups)` — the ranked estimator

```
A_wa = sum_a w_a * C_a    C_a = AUROC within stratum a
                          w_a = (n+_a * n-_a) / sum_b (n+_b * n-_b)
```

Strata are alleles. Admissible when `n+_a >= min_per_class` and `n-_a >= min_per_class` (default 1). Inadmissible strata take weight zero and are **excluded, never imputed to 0.5**; they stay in `per_stratum` and in `n_strata_total`. A single-class stratum is not an error — it returns `auroc: None`, `pair_weight: 0.0`.

Pair weight, not an unweighted macro, because the alleles are unequal. Infectious set (n = 2,093, 78 alleles): 36 strata admissible over 1,787 rows (85.4%); three alleles hold 49% of evaluable rows; 14 of the 36 hold 66 rows — 3.7% of rows but **38.9% of the weight under an unweighted macro**. One arm: macro **0.580835** vs pair-weighted **0.562251**. Not rank-then-pool either (see below).

Returns a `StratifiedResult`, never a bare float: `0.531844  strata 36/78  rows 1787/2093  H=0.0863`

Fields: `value` (`A_wa`), `n_strata_contributing`/`n_strata_total`, `n_rows_evaluable`/`n_rows_total`, `evaluable_fraction`, `herfindahl` (`sum_a w_a^2`, where 1.0 means one stratum is the whole estimate), `top3_pair_weight`, `provisional` (`True` when `top3_pair_weight > 0.5`, printed `[PROVISIONAL]`), and `per_stratum` (per allele `n`, `n_pos`, `n_neg`, `admissible`, `auroc`, `pair_weight`).

### Guards — all raise, none warn

| condition | the number it prevents |
|---|---|
| any empty or `None` group key | every row lands in one group, within-allele silently becomes pooled: **0.584138 reported where the correct figure was 0.492317**. Row index named |
| all rows share one group key | same failure; the estimate would equal pooled AUROC to machine precision |
| every group is a singleton | **exactly 0.500000 regardless of the scores** — reads "no signal", means "no grouping" |
| no admissible stratum | nothing to average; message gives `min_per_class`, group count, row count |
| length mismatch | `ValueError`, before any statistic is computed |

Each guard ships a fixture firing it *and* asserting the unguarded result: `test_missing_grouping_raises_and_would_have_equalled_pooled` checks the bypassed path equals pooled AUROC to within 1e-12.

## Pooled AUROC

`auroc` over all rows. Printed always, **never ranked on**. An allele lookup table with **no peptide input** scores pooled **0.7122** at within-allele exactly **0.5000** (constant inside a stratum), so read pooled above ~0.70 as a statement about allele sampling. The estimators also reverse rankings on identical rows: hydropathy 0.5264 within-allele / 0.5650 pooled; PRIME-2.0 0.5205 / 0.5663.

## `mixture_decomposition(labels, scores, blocks)`

Exact identity over any block partition, not an approximation:

```
A_pool = sum_{j,k} w_jk * C_jk    w_jk = (n+_j * n-_k) / (N+ * N-)
                                  C_jk = P(pos in j > neg in k), ties 0.5
```

Diagonal `C_jj` are within-block AUROCs; off-diagonal terms are cross-block pairs no within-block estimator considers. Raises when no block carries both classes, or when the whole set is single-class. Fields: `pooled`, `diagonal_mass` (`sum_j w_jj`), `diagonal_only` (weight-renormalised diagonal), `off_diagonal_mean`, `max_constituent` (`max_j C_jj`), `excess` (`pooled - max_constituent`), `artifact` (`excess > 0`), `blocks`, `cross`.

**Rule: `excess > 0` means the pooled figure exceeds every block it is composed of, and is flagged a mixture artifact. The flag is on `max`, not on a weighted mean.** A weighted average cannot exceed its largest term; a mixture can, because block offsets order cross pairs for free. Measured: pathogen 0.5806, self 0.5727, pooled 0.5837, excess **+0.0031**. Printed inline:

```
pooled 0.572103  max constituent 0.569881  excess +0.002222  diag mass 0.5140  <-- MIXTURE ARTIFACT
```

`stratified_auroc` is diagonal-only and immune; a fixture asserts it stays at or below 0.625 on a hand-computed case whose pooled figure is 0.6667.

## The chance bar

**Chance is 0.5000, exactly. 0.522 is a 5% critical value, not a chance level** — the point estimate needed to exclude 0.5 at n = 2,093 under a one-sided Hanley–McNeil normal approximation; at n = 833 the same construction gives ~0.543. A model at 0.515 is above chance *and* not distinguishable from chance. The 0.543 once printed in the same table row is half a bootstrap CI width: a different construction, not one series. Printed separately, constructions named:

```
chance                                        0.5000    (exactly; not 0.522)
5% critical value at this n (Hanley-McNeil)   0.522     one-sided, excludes 0.5
```

**MDE for a paired difference at n = 2,093 is ~0.059–0.095** (effective n 1,217–1,903 after the measured source-protein design effect of 1.10–1.72). A margin below it is a tie, not a win — one re-measured margin moved from "+0.0246, p = 0.078" to **+0.0187 [−0.0111, +0.0466], p = 0.215**. A pre-registered threshold below the MDE is refused in advance. Significance is a **paired** bootstrap over **source-protein clusters**: resample once per replicate, score every model on that resample, read percentiles of the per-replicate difference. Measured widening from correct clustering is **1.00–1.08×**; this project's earlier "20–121% too narrow" claim for row-level intervals is **withdrawn** and must not be cited, including as support. Until the cluster bootstrap ships, no interval prints. Any interval carries its basis `{row|cluster}/B/{saved|fresh}/seed/assigned-fraction`; mixed-basis tables are rejected.

## `pooled_rank_auroc(labels, scores, groups)` — diagnostic only

Rank-normalise within group to `(midrank - 0.5) / n_group`, then one pooled AUROC over all rows. This is what much of the literature calls "within-allele AUROC"; kept **only** to reproduce historical figures, and **never reportable as a within-allele number**. Cross-group pairs still contribute (pooled 0.5837 above both constituents, pathogen 0.5806 / self 0.5727), and singleton and single-class groups are injected at 0.5 instead of excluded. Degenerate inputs return publishable-looking numbers: **0.584138** under a uniformly missing grouping where the correct figure is 0.492317, **exactly 0.500000** under an all-singleton grouping.

Same arm, three estimators: **0.554016** rank-then-pool, **0.562251** pair-weighted, **0.580835** macro. On neoantigen the same scores give 0.486 rank-then-pool and 0.513 macro — below and above chance. Both intervals span chance, so the conclusion holds, but **direction may never be stated without naming the estimator**. True under all three: "does not demonstrate within-allele signal; the interval includes chance".

## Pooled and stratified ship together, or neither ships

> Report the stratified and the pooled estimator together, or report neither. A pooled figure that exceeds every constituent is a mixture artifact and must print its decomposition beside it.

The same scores have produced pooled 0.6304 against within-allele 0.5306 here. The scorecard has no mode printing one alone, `epigauge score` has no flag to suppress either, and an artifact flag prints in the same block as the pooled figure, never in a footnote.

## Keys and zeros

`epigauge.canon` is the only thing permitted to normalise a key, and everything in it raises rather than falling back. `canonicalise_allele(raw, *, fields=2)` maps `'A0101'`, `'HLA-A*01:01'` and `'hla-a0101'` to `'HLA-A*01:01'` (two-field by default and it says so; expression suffixes `N/L/S/C/A/Q` dropped, they do not change the groove; mouse `H-2-Kb` → `H2-Kb` passes through). `canonicalise_organism(raw)` is exactly "strip one trailing balanced parenthetical, lowercase, collapse whitespace", nothing else. `canonicalise_peptide(raw, *, min_len=8, max_len=15)` uppercases, strips, and asserts the canonical 20 residues and the length window.

- **No fallback values.** A missed allele join produced exactly 0.5000 and was caught only because 0.5000 looks impossible; a global-rate fallback would have produced a publishable 0.52. Unparseable input raises `CanonError`; a tolerant caller catches it and records the drop.
- **No substring matching for organisms.** `Homo sapiens` and `Homo sapiens (human)` were audited as two blocks — 8,857 rows at 48.4% positive, 7,492 at 31.8%, 180 studies using both spellings. The substring `human` occurs in 16 organism strings here and **13 are viruses** (human betaherpesvirus 5, 85% positive).
- **Non-canonical residues raise**, because every adjacency and clustering metric is defined over canonical 20-letter sequences. Allele-spelling mismatches are also how a contaminated benchmark hid in plain sight — see [CONTAMINATION.md](CONTAMINATION.md).
- **Any printed zero comes from a counter a fixture proves can be non-zero.** `0.00e+00` and `0 of 2,093` must never mean "nothing looked"; every zero prints its denominator.
