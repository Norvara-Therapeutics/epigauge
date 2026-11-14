# A public dataset file contains the whole of a published neoantigen benchmark

`datasets/im_test.csv` from the BigMHC dataset release (Mendeley Data `10.17632/dvmz6pkzvb`, **version 2**, CC BY 4.0) contains all 833 rows of the TITANiAN neoantigen immunogenicity benchmark, same labels, reference-predictor columns reproduced character for character. One file is the other with 104 extra rows and eight extra columns. The naive check returns a clean-looking zero because one file spells the allele `HLA-A*01:01` and the other `A0101`.

## Reproduce it

```
python3 scripts/verify_contamination.py \
    --a <path>/bigmhc/datasets/im_test.csv \
    --b <path>/TITANiAN_IM_Neoantigen.csv
```

Standard library only; run from the `epigauge/` directory so `epigauge/canon.py` resolves. Both paths are required — the script never guesses where a file lives. Exit 1 = overlap found, 0 = independent on this evidence, 2 = could not run. It never guesses a column name and never substitutes a fallback join key: if it cannot identify the allele, peptide or label column it refuses and prints the header it was given.

### The two files

**Neither ships with EpiGauge** — both are third-party releases, and this repository
redistributes neither.

| | file | source | sha256 |
|---|---|---|---|
| `--a` | `datasets/im_test.csv`, 937 rows, 279,819 B | inside `datasets.zip` at Mendeley Data **10.17632/dvmz6pkzvb version 2**, CC BY 4.0 | `70ba5dbeecb0e62ccc95f2aa11a28801773893f4ebf209b5e68b58318be7b634` |
| `--b` | `TITANiAN_IM_Neoantigen.csv`, 833 rows, 328,449 B | inside `test.zip` at Dryad **10.5061/dryad.s7h44j1k7**, CC0 1.0 | `da2442cd…1dfed9661` |

Both DOIs resolve and both licences are open, so anyone can obtain these. Check what you
download against the hashes above; the figures on this page are from exactly those bytes.

Both are the original CSVs. EpiGauge ships the same 833 rows in normalised form at
`benchmarks/cc0/titanian_neoantigen.jsonl`, but the check reads the originals, because it
compares the predictor columns **as bytes** and normalisation would destroy exactly what it
is measuring.

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
- A score on that arm is not a capability bar for a method with access to that table. Whether any particular method had that access is not established here, and nothing about how either file came to exist is measured.

**Scope.** `verify_contamination.py` compares two CSVs on a canonicalised `(allele, peptide)` key and reports pair overlap, label agreement and per-column byte equality. It measures no model's training data, computes no AUROC, evaluates no predictor, and does not adjudicate how either file was assembled.

## Which deposit version these bytes came from

The archive this project extracted `im_test.csv` from is **version 2** of the Mendeley
deposit, not version 4. Established by digest rather than by reading a version number:
the deposit publishes four versions, `datasets.zip` is absent from v1, is 251,406,808
bytes with sha256 `33a51710aa55ff1af99944ee617fc9f5e2a5224f4dfd2df85dd36d6fe7f44c5d` in
v2, and is 144,274,974 bytes with sha256 `0d152a452756cf2e0014ffced6afc25118e7c11cf1f626b26e49f50f79edffaa`
in both v3 and v4. The local archive's recorded size and digest match v2 uniquely.

These docs previously cited version 4. That was wrong for these bytes.

**What this does not establish.** v3 and v4 share one archive digest, so `im_test.csv`
may well be byte-identical across v2, v3 and v4 — that has not been checked, because it
needs the v4 archive downloaded and its member hashed. The correct statement today is
that the local copy is v2's member, **not** that v4 lacks these bytes.
