# ITSNdb — segregated slice, read this before redistributing

Source: `github.com/elmerfer/ITSNdb`, `data/ITSNdb.csv`
(sha256 `f8b4b1622a8e38fec72788f7ed17c7a23208f77f8c5c78ce936e6e7632c66603`).
Nibeyro et al. 2023. 199 rows, 129 positive / 70 negative, 31 alleles.

**The licence position is unresolved, and this notice does not assert a right you can rely
on.** The repository's LICENSE is **GPL-3.0 boilerplate with no named copyright holder**,
and **the LICENSE covers the R package; the data records carry no separate licence.**

- This slice lives in **its own directory**, under its own notice, and is **not** part of
  the CC BY 4.0 bundle. Do not merge it into a CC BY 4.0 redistribution.
- If you redistribute it onward, form your own view of the licence first, or contact the
  ITSNdb authors. Nothing here grants you a licence.
- For a release with no unresolved licence position, use `cc-by-4.0/` and `cc0/` only.

## Known data discrepancy

ITSNdb's `Length` column reports 149 nine-mers and 50 ten-mers; the sequences are 150 and
49 (row `IIGAGPAEV` is 9 residues and carries `Length = 10`). **The sequence is used**, so
this file reports 150/49. The original value is preserved in each row's
`extra.declared_length`.
