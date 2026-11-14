# EpiGauge leaderboard — infectious benchmark

Pair-weighted within-allele AUROC for class-I T-cell epitope immunogenicity prediction, scored on
one specified estimator across every method EpiGauge ships a cleared score column for.

**2,093 rows · 1,488 positive / 605 negative (71.1%) · 78 alleles**
Estimator: pair-weighted within-allele AUROC over **36 of 78 admissible strata** (1,787 / 2,093 rows).
Minimum detectable paired difference at this n: **≈0.059–0.095**.

> **Gaps smaller than the MDE are not rankings.** Nearly every gap in this table is smaller than the
> MDE. Read the column as a band, not an order: entries within ~0.06 of each other are not
> distinguishable on this benchmark, whatever their position here.

## The table

Ranked on within-allele. Pooled is printed but never ranked on.

**Key to the markers**

| marker | meaning |
|---|---|
| ⚠ | **Not a capability measurement.** The label is recoverable from public data, or the benchmark rows are contained in the method's training set. Details [below the table](#why-the--rows-are-not-capability-measurements). |
| † | **Figure transcribed, not recomputed.** The score comes from a predictor column distributed inside the evaluation table itself, so no inference was run here. This is a statement about where the number came from, not about what the method saw. |
| **bold row** | Reference baseline — the bar a method has to clear before its number means anything. |

<!-- BEGIN GENERATED TABLE -->

| # | within-allele | pooled | entry |
|---|---|---|---|
| 1 | 0.7444 | 0.7844 | TITANiAN-IM ⚠ † |
| 2 | 0.5742 | 0.6053 | TITANiAN-IM ablation 1,2,3,4,7 ⚠ † |
| 3 | 0.5712 | 0.5952 | NORVARA-3B — frozen ESM-2 3B backbone, 2,561-param head |
| 4 | 0.5675 | 0.6106 | NORVARA-650M — frozen ESM-2 650M backbone, 1,281-param head |
| 5 | 0.5664 | 0.6169 | NORVARA-150M — frozen ESM-2 150M backbone, 641-param head |
| 6 | 0.5657 | 0.6210 | TITANiAN-IM ablation 3,4,5,6,7 ⚠ † |
| 7 | 0.5445 | 0.5758 | TITANiAN-IM ablation 1,2,5,6,7 ⚠ † |
| 8 | 0.5328 | 0.5687 | BigMHC canonical-LSTM IM ablation † |
| 9 | 0.5319 | 0.5669 | BigMHC canonical-LSTM ELIM ablation † |
| 10 | 0.5267 | 0.5934 | TITANiAN-IM ablation 1,2,3,4,5,6 ⚠ † |
| **11** | **0.5264** | 0.5650 | **BASELINE — mean Kyte-Doolittle hydropathy, unfitted** |
| 12 | 0.5205 | 0.5663 | PRIME-2.0 (scores) † |
| 13 | 0.5205 | 0.5518 | PRIME-2.0 (percentile ranks, sign-corrected) † |
| 14 | 0.5195 | 0.5673 | BigMHC_IM ⚠ † |
| 15 | 0.5189 | 0.5692 | BigMHC_ELIM ⚠ † |
| 16 | 0.5180 | 0.5614 | PRIME-1.0 (scores) † |
| 17 | 0.5175 | 0.5540 | PRIME-1.0 (percentile ranks, sign-corrected) † |
| 18 | 0.5162 | 0.5555 | MixMHCpred-2.1 (percentile ranks, sign-corrected) † |
| 19 | 0.5148 | 0.5518 | MixMHCpred-2.1 (scores) † |
| 20 | 0.5129 | 0.5541 | MHCflurry-2.0 (scores) † |
| 21 | 0.5129 | 0.5541 | MHCflurry-2.0 (percentile ranks, sign-corrected) † |
| 22 | 0.5119 | 0.5608 | MixMHCpred-2.2 (scores) † |
| 23 | 0.5119 | 0.5567 | MixMHCpred-2.2 (percentile ranks, sign-corrected) † |
| 24 | 0.5088 | 0.5652 | BigMHC wide-LSTM IM ablation † |
| 25 | 0.5080 | 0.5632 | BigMHC wide-LSTM ELIM ablation † |
| **26** | **0.5000** | 0.7122 | **BASELINE — allele lookup, no peptide input** |
| 27 | 0.4969 | 0.4540 | MHCnuggets-2.4.0 † |
| 28 | 0.4901 | 0.5599 | TITANiAN-IM baseline ⚠ † |
| 29 | 0.4879 | 0.5465 | TransPHLA † |
| 30 | 0.4813 | 0.5421 | BigMHC_EL † |
| 31 | 0.4803 | 0.5460 | BigMHC canonical-LSTM EL ablation † |
| 32 | 0.4746 | 0.5417 | BigMHC wide-LSTM EL ablation † |
| 33 | 0.4714 | 0.5233 | HLAthena (percentile ranks, sign-corrected) † |
| 34 | 0.4705 | 0.5048 | HLAthena (scores) † |
| **35** | **0.4132** | 0.4512 | **BASELINE — peptide length** |
| 35-36 | withheld | withheld | NetMHCpan-4.1 (scores) † |
| 35-36 | withheld | withheld | NetMHCpan-4.1 (percentile ranks) † |

<!-- END GENERATED TABLE -->

Ranks 35–36 are two rows of the machine-readable table — NetMHCpan-4.1 scores and percentile ranks —
withheld under the same licence term and collapsed into one line here.

## Notes on the marked rows

### Why the ⚠ rows are not capability measurements

**TITANiAN-IM and its ablations:** this column's label is recoverable from a public IEDB assay table
at within-allele 0.9101 with 100% coverage.

**BigMHC_IM / BigMHC_ELIM:** `im_test.csv` contains this benchmark's neoantigen arm in full — 833 of
833 rows, same labels, 25 of 25 shared predictor columns byte-identical. Verify it yourself with
[`scripts/verify_contamination.py`](https://github.com/norvara-therapeutics/epigauge/blob/main/scripts/verify_contamination.py).

---

## Three things this table shows

**An unfitted 1982 physicochemical scale places 11th.** Mean Kyte-Doolittle hydropathy — zero free
parameters, nothing trained — outranks PRIME-2.0, MHCflurry-2.0, MixMHCpred, NetMHCpan-4.1,
TransPHLA, MHCnuggets and six BigMHC variants.

**A lookup with no peptide input scores pooled 0.7122.** Its within-allele figure is exactly 0.5000,
because the score is constant inside every stratum. On a pooled leaderboard it would rank second.
This is why EpiGauge prints pooled but never ranks on it.

**The honest reading of the whole column is "a band near 0.52".** Excluding the flagged rows, every
non-baseline entry sits between 0.47 and 0.54, and the spread is smaller than the minimum detectable
difference. No method here is measurably better than any other on this benchmark.

## Reproduce it

```bash
pip install epigauge
epigauge leaderboard --set titanian_infectious
```

Scores for every shipped column, the row set and the per-file checksums are in
[`benchmarks/`](https://github.com/norvara-therapeutics/epigauge/tree/main/benchmarks) in the
repository. The estimator, the degeneracy guards and the chance-bar construction are specified in
[METRICS](../docs/METRICS.md); the submission contract is in [SUBMITTING](../docs/SUBMITTING.md).
Machine-readable form of this table: [leaderboard.csv](infectious.csv).
