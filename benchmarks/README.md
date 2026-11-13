# epigauge benchmark data

Everything here may lawfully be redistributed. `MANIFEST.json` is the authority: per set it
carries the file path, `sha256`, byte count, row count, class balance, licence, the **basis**
for that licence, the source DOI/URL and the retrieval date, plus a per-column clearance
table for every published predictor score.

Load through the package, not by reading files directly — `load` verifies each file's
`sha256` against `MANIFEST.json` before returning a row:

```python
from epigauge.data import available, load, fetch
available()                                       # every set_id, with status and licence
rows = load("titanian_infectious")                # verifies sha256 before yielding a row
rows = load("titanian_infectious", scores=True)   # rows + cleared competitor columns
```

## What is here

| directory | licence | sets |
|---|---|---|
| `cc-by-4.0/` | CC BY 4.0 | `independent_composite`, `hbv_hlabc_screen`, `xmuimmu_noncanonical` |
| `cc0/` | CC0 1.0 (Dryad `10.5061/dryad.s7h44j1k7`) | `titanian_infectious`, `titanian_neoantigen` (+ `.scores.jsonl` sidecars) |
| `itsndb/` | **segregated — see `itsndb/NOTICE.md`** | `itsndb` |

Each directory's `NOTICE.md` carries the operative licence statement and the attribution
that licence requires. `independent_composite` is exactly the union of `hbv_hlabc_screen`
and `xmuimmu_noncanonical`; do not pool it with its own parts.

## Availability of the remaining sets

`MANIFEST.json` lists every withheld set with the operative licence clause. `fetch()`
raises `LicenceError` naming the clause; it never returns an empty result.

- **fetch-only** (not redistributable; `fetch()` pins the upstream `sha256`) — `tesla608`,
  `buckley_gbm`, `titanian_mhc1`, `titanian_mhc2`, `titanian_tcr_activation`,
  `titanian_tcr_zeroshot`, `titanian_ada`.
- **refused** — `vdjdb` (CC BY-NC-ND; ND bars derived rows), `mcpas_tcr` and
  `neoranking_muller` (no licence), `hitide_ega` and `tesla_synapse` (controlled access),
  `tchard` (non-commercial, per-row licence column), `tscape_titanian` (CC BY-NC copy of
  data shipped here under CC0).

## Predictor score columns

Of ~27 published competitor columns per arm, nine clear for redistribution. **NetMHCpan
columns are withheld** under DTU §7(v) (*"Licensee may not publish any results of benchmark
tests run on the Product to a third party without HEALTH's prior written consent"*); no
waiver exists in writing. BigMHC, MHCnuggets, MixMHCpred, PRIME and HLAthena columns are
withheld as non-commercial, unlicensed or web-service-only. Withheld columns are **absent**,
not nulled. `MANIFEST.json → predictor_columns` gives the verdict and reason for all 82.

## Packaging

This directory sits beside the `epigauge/` package. `epigauge.data` finds it via, in order:
`$EPIGAUGE_BENCHMARKS`, `epigauge/benchmarks/` inside the package, then `../benchmarks`.
Packaging config must ship it as package data for offline `pip install epigauge`.

## Citation

Cite the upstream source of every set you report on, not just this package.
`MANIFEST.json` carries the DOI for each.
