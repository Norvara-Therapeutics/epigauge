# CC BY 4.0 data

The row sets in this directory derive from sources licensed CC BY 4.0:

- IEDB / CEDAR v3 exports. Operative text from `iedb.org/terms_of_use_v3.php`:
  *"NIAID places no restrictions on the use or distribution of the data contained
  within the IEDB"*, with CC BY 4.0 stated beside the downloads.
- `hbv_hlabc_screen` — Frontiers in Immunology, open access, CC BY 4.0,
  doi:10.3389/fimmu.2025.1717543.
- `xmuimmu_noncanonical` — XMUimmu non-canonical immunopeptidome database,
  doi:10.1126/science.adk3487.

You may share and adapt these rows for any purpose, including commercially, provided you
give appropriate credit, link to the licence and indicate changes.

**Changes made:** rows were normalised into the schema in `MANIFEST.json → row_schema`,
and rows whose allele could not be resolved were dropped (80 in `hbv_hlabc_screen`).

Licence text: https://creativecommons.org/licenses/by/4.0/
