# Neoantigen benchmark — no method performs above chance

**833 rows · 178 positive / 655 negative (21.4%) · 49 alleles**
Estimator: pair-weighted within-allele AUROC over **22 of 49 admissible strata** (605 / 833 rows).

> **This is a null result, not a ranking.** All 35 methods fall between **0.405 and 0.540**, with
> chance at 0.500 — **22 of them below it**. The ordering below is sorted, but the spread is
> narrower than this benchmark can resolve. Nothing here is measurably better than anything else,
> including the entries submitted by this benchmark's own authors, which rank 25th, 34th and 35th.

## Read these three before the table

**Every row is contained in another public dataset.** All 833 are present in `im_test.csv` in the
BigMHC Mendeley deposit — same labels, 25 of 25 shared predictor columns byte-identical. Any method
in that lineage may have trained on them. See [CONTAMINATION.md](../docs/CONTAMINATION.md).

**The labels barely support learning.** A model trained *on these very rows* reaches only 0.5535.
That is the ceiling for anything fitted here, and it is close to chance.

**Less than three quarters of the set carries weight.** 27 of 49 alleles are singleton or
single-class, so 228 of 833 rows are excluded from the estimator rather than imputed.

<!-- BEGIN GENERATED TABLE -->

| # | within-allele | pooled | entry |
|---|---|---|---|
| **1** | **0.5397** | 0.5201 | **BASELINE — peptide length** |
| 2 | 0.5265 | 0.5256 | MixMHCpred-2.1 (percentile ranks, sign-corrected) |
| 3 | 0.5242 | 0.5240 | MixMHCpred-2.1 (scores) |
| 4 | 0.5234 | 0.4788 | MHCflurry-2.0 (scores) |
| 5 | 0.5234 | 0.5383 | MixMHCpred-2.2 (scores) |
| 6 | 0.5234 | 0.5522 | MixMHCpred-2.2 (percentile ranks, sign-corrected) |
| 7 | 0.5230 | 0.4788 | MHCflurry-2.0 (percentile ranks, sign-corrected) |
| 8 | 0.5118 | 0.5215 | HLAthena (percentile ranks, sign-corrected) |
| 9 | 0.5102 | 0.5100 | BigMHC_EL |
| 10 | 0.5081 | 0.5024 | HLAthena (scores) |
| 11 | 0.5075 | 0.5159 | BigMHC_Ablation_Wide_LSTM_EL |
| 12 | 0.4955 | 0.4742 | TransPHLA |
| 13 | 0.4874 | 0.5568 | MHCnuggets-2.4.0 |
| 14 | 0.4835 | 0.5135 | TITANiAN-IM-baseline ⚠ |
| 15 | 0.4827 | 0.5039 | BigMHC_Ablation_Canonical_LSTM_EL |
| 16 | 0.4827 | 0.6048 | TITANiAN-IM ⚠ |
| 17 | 0.4800 | 0.5107 | TITANiAN-IM-ablation-1,2,3,4,5,6 ⚠ |
| 18 | 0.4751 | 0.5352 | PRIME-2.0 (scores) |
| 19 | 0.4751 | 0.5452 | PRIME-2.0 (percentile ranks, sign-corrected) |
| 20 | 0.4724 | 0.5026 | PRIME-1.0 (percentile ranks, sign-corrected) |
| 21 | 0.4722 | 0.4998 | PRIME-1.0 (scores) |
| 22 | 0.4639 | 0.4801 | TITANiAN-IM-ablation-1,2,5,6,7 ⚠ |
| 23 | 0.4556 | 0.4663 | NORVARA-150M |
| **24** | **0.4497** | 0.4681 | **BASELINE — mean Kyte-Doolittle hydropathy, unfitted** |
| 25 | 0.4489 | 0.5646 | TITANiAN-IM-ablation-1,2,3,4,7 ⚠ |
| 26 | 0.4438 | 0.5524 | BigMHC_ELIM ⚠ |
| 27 | 0.4434 | 0.5459 | BigMHC_Ablation_Canonical_LSTM_ELIM |
| 28 | 0.4413 | 0.5094 | BigMHC_IM ⚠ |
| 29 | 0.4383 | 0.5157 | BigMHC_Ablation_Canonical_LSTM_IM |
| 30 | 0.4377 | 0.5804 | BigMHC_Ablation_Wide_LSTM_ELIM |
| 31 | 0.4309 | 0.5180 | BigMHC_Ablation_Wide_LSTM_IM |
| 32 | 0.4203 | 0.4669 | NORVARA-3B |
| 33 | 0.4053 | 0.4253 | NORVARA-650M |
| 33-34 | withheld | withheld | NetMHCpan-4.1 (scores) |
| 33-34 | withheld | withheld | NetMHCpan-4.1 (percentile ranks) |

<!-- END GENERATED TABLE -->

⚠ marks a column that is not a capability measurement; the reasons are
[below](#why-the--rows-are-not-capability-measurements).

---

## What this arm shows that the infectious arm cannot

### Peptide length ranks first here and last on infectious

A single integer — the length of the peptide — tops this table at **0.5397**, and places **37th of
37 at 0.4132** on [the infectious benchmark](infectious.md). A feature cannot be the best predictor
of immunogenicity on one set and the worst on another. It is reading how each set was assembled.

### TITANiAN-IM falls 0.744 → 0.483

| | infectious | neoantigen |
|---|---|---|
| TITANiAN-IM | 0.7444 | **0.4827** |
| NORVARA-150M | 0.5664 | 0.4556 |
| PRIME-2.0 | 0.5205 | 0.4751 |
| MHCflurry-2.0 | 0.5129 | 0.5234 |

On the infectious arm, that column's label is recoverable from a public IEDB assay table at
within-allele **0.9101 with 100% coverage**. Here the same lookup covers about **1.3%** of rows, and
the column drops below chance — a fall of 0.26, far larger than any other method's. This is the
control for the capability flag on the infectious page: it is a measurement, not an assertion.

### The estimator inversion is starker here

TITANiAN-IM scores **pooled 0.6048** against **within-allele 0.4827**. A leaderboard ranked on
pooled would place a below-chance column first. This is why EpiGauge prints pooled and never ranks
on it.

### Why the ⚠ rows are not capability measurements

**TITANiAN-IM and its ablations** — the label is recoverable from a public IEDB assay table, as
above.

**BigMHC_IM / BigMHC_ELIM** — `im_test.csv` contains this entire arm. Reproduce it with
`python3 scripts/verify_contamination.py`.

---

## Reproduce

```bash
pip install -e .
epigauge leaderboard --set titanian_neoantigen
```

Machine-readable: [neoantigen.csv](neoantigen.csv). The estimator and its
guards are specified in [METRICS](../docs/METRICS.md); the infectious arm is
[here](infectious.md).
