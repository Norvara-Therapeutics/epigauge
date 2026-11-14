# CC0 1.0 data

`titanian_infectious` and `titanian_neoantigen` derive from the TITANiAN benchmark tables
deposited at Dryad, doi:10.5061/dryad.s7h44j1k7, which the Dryad API reports as
`spdx: CC0-1.0`.

The local bytes match the project's download manifest. That manifest records a digest and
a manifest-generation date, **not a retrieval from this deposit**: the local CSVs are
dated 2025-05-06 and the Dryad deposit published 2025-11-20, so the shipped rows predate
it and cannot be its published bytes. Dryad exposes a sha256 for `test.zip`
(`cefad8a75a5b7ea9100f08820cb66df24d84d0c39de266b688d678219562ed94`, 518,023 bytes) but
no member-level digest, and its download endpoint refuses unauthenticated requests, so
the archive has not been compared against the shipped rows. Doing so is one manual
download; until then, treat the upstream correspondence as asserted rather than
verified.

CC0 waives copyright; attribution is not a legal condition. **Cite the TITANiAN authors anyway** — a scientific norm, not a licence term.

The `seoklab/T-SCAPE` repository is CC BY-NC 4.0 and holds no copy of these rows — it carries
code, weights and HLA pseudosequence tables. The Dryad deposit is the only public source of this
evaluation set, and these rows come from it.

Licence text: https://creativecommons.org/publicdomain/zero/1.0/

## Two warnings about `titanian_neoantigen`

1. **Contamination.** These 833 rows are the same evaluation table as `im_test.csv` in the
   BigMHC Mendeley deposit: 833/833 `(allele, peptide)` pairs, 833/833 labels agreeing, all
   25 shared predictor columns byte-identical. The two files spell the allele differently, so a
   naive join finds nothing. A score here is not a capability bar for any model with access to
   that table.
2. **Open question.** This arm has exactly 178 positives and Müller et al.'s harmonised
   benchmark exactly 178 immunogenic neo-peptides; no peptide-level intersection has been
   run, so independence is not established.
