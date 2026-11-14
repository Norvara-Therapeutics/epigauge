# CC0 1.0 data

`titanian_infectious` and `titanian_neoantigen` derive from the TITANiAN benchmark tables
deposited at Dryad, doi:10.5061/dryad.s7h44j1k7, which the Dryad API reports as
`spdx: CC0-1.0`.

**Upstream identity is verified, not asserted.** The deposit was downloaded on
2026-11-14 and `test.zip` matches the sha256 Dryad publishes
(`cefad8a75a5b7ea9100f08820cb66df24d84d0c39de266b688d678219562ed94`, 518,023 bytes). All
seven shipped CSVs are **byte-identical** to the archive's members, including
`TITANiAN_IM_Infectious.csv` (843,836 B) and `TITANiAN_IM_Neoantigen.csv` (328,449 B).

> An earlier version of this notice stated that the shipped rows predated the deposit and
> "cannot be its published bytes", reasoning from the local files' 2025-05-06 modification
> time against the deposit's 2025-11-20 publication date. **That was wrong.** A zip
> preserves its members' mtimes, so extracting this deposit yields files stamped when its
> authors created them. A modification time is not provenance; a digest against a named
> upstream is.

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
