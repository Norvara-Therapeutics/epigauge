# A public dataset file contains the whole of a published neoantigen benchmark

`datasets/im_test.csv` from the BigMHC dataset release contains all 833 rows of the TITANiAN neoantigen immunogenicity benchmark, same labels, reference-predictor columns reproduced character for character. One file is the other with 104 extra rows and eight extra columns. The naive check returns a clean-looking zero because one file spells the allele `HLA-A*01:01` and the other `A0101`.

## Reproduce it

```
python3 scripts/verify_contamination.py \
    --a <path>/bigmhc/datasets/im_test.csv \
    --b <path>/TITANiAN_IM_Neoantigen.csv
```

Standard library only; run from the `epigauge/` directory so `epigauge/canon.py` resolves. Both paths are required — the script never guesses where a file lives. Exit 1 = overlap found, 0 = independent on this evidence, 2 = could not run. It never guesses a column name and never substitutes a fallback join key: if it cannot identify the allele, peptide or label column it refuses and prints the header it was given.

| step | measurement |
|---|---|
| naive join on the raw allele strings | **0 / 833** |
| same join after `canonicalise_allele` | **833 / 833** |
| label agreement on matched rows (`label` vs `tgt`) | **833 / 833** (100.000%) |
| shared predictor columns byte-identical across all 833 rows | **25 / 25** |
| shared predictor columns in the same order in both headers | yes |
| rows of `im_test.csv` not in the benchmark | **104** (937 − 833) |
| positives among the matched rows | 178 |

Byte-identical means the CSV text as written, not parsed floats and not a tolerance: `MHCflurry-2.0_Scores` reads `0.9662107380862464` in both files for peptide `ALDKLSSQHLY`.


## Scope

This check compares files. It establishes that one file contains another's rows, with the same
labels and the same predictor columns, and therefore that a score on that table is not evidence of
capability on data the scorer has not seen.

- The claims are about files: a file *contains* rows, a label *is recoverable from* a prior table.
- 0.704 is not reachable without access to that assay table, so it is not a capability bar.

## Files referenced, with hashes

SHA-256 of the file as held locally. Retrieval of record is the download manifest, written **2026-11-10T01:36:15Z**. The dates below are local modification times, not retrieval dates: the Dryad deposit ships a zip and zip preserves member mtimes.

Upstream data hashes are the evidence and are pinned. The two repo-internal rows mark the
version that produced the figures above; they move whenever those files are edited, so treat
git as authoritative for them.

| file | sha256 | size | local mtime | source |
|---|---|---|---|---|
| `bigmhc/datasets/im_test.csv` | `70ba5dbeecb0e62ccc95f2aa11a28801773893f4ebf209b5e68b58318be7b634` | 279,819 B | 2026-11-08 | extracted from the archive below |
| `bigmhc/datasets.zip` (container) | `33a51710aa55ff1af99944ee617fc9f5e2a5224f4dfd2df85dd36d6fe7f44c5d` | 251,406,808 B | 2026-11-08 | Mendeley Data **10.17632/dvmz6pkzvb.4** (v4, 2023-04-14), CC BY 4.0 |
| `test/TITANiAN_IM_Neoantigen.csv` | `da2442cd7f7ea0d1ed83e15a71538cb0f8983a329fe973cfe43514c1dfed9661` | 328,449 B | 2025-05-06 | Dryad **10.5061/dryad.s7h44j1k7**, CC0 1.0 |
| `iedb/tcell_full_v3.zip` | `3863e35be7730f3e594ca4f7db3a916595e05f536439090ebd5db7b23047f1f3` | 45,275,246 B | 2026-11-10 | `iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip`, CC BY 4.0 |
| `data/raw/download_manifest.json` | content digest `eda0e7aaaeb24b471ce62e8883b355e01ac79b443ff452489d5e9d55428ba094` | — | 2026-11-10 | written by this repository |
| `scripts/verify_contamination.py` | `b03f619ecc9f08cb4c2f281791b7a0c933c87b1b252dd4bde1dfb763755d2460` | — | — | this repository |
| `epigauge/canon.py` | `faf9a1fbd96d6c2534c4dbf37b517c0a32678575300fc7f9cbfeeb85340e438a` | — | — | this repository |

The TITANiAN files exist in two licensed copies — the Dryad deposit (CC0 1.0) and a T-SCAPE copy (CC BY-NC 4.0). The hash above is the Dryad copy, and the Dryad copy is the one to cite. Both source datasets are openly licensed and both DOIs resolve.

The authors of both files are given this document, `verify_contamination.py` and the full row-level join output before publication. That does not gate publication, and any correction is published in the same place with the original left visible beside it.

**Scope.** `verify_contamination.py` compares two CSVs on a canonicalised `(allele, peptide)` key and reports pair overlap, label agreement and per-column byte equality. It measures no model's training data, computes no AUROC, evaluates no predictor, and does not adjudicate how either file was assembled.
