# `supplementary/`

Supplementary data for the epigauge preprint. Every artifact the paper's supplementary
index cites is here, together with the analysis code that a third party can run against
the data this package ships.

Paths in this directory and inside every artifact are relative to the package root — the
directory holding `pyproject.toml`, `leaderboards/` and `benchmarks/`. No artifact
records an absolute path. Where an artifact was computed from an input that is not
redistributed, the `script` or input field names the deposit and the `sha256` of the
exact file used, which is how a third party confirms it holds the same input.

```
supplementary/
  artifacts/              14 JSON records, 4.2 MB
    annotation/           the 7 annotation-ceiling records
  scripts/                5 runnable analysis scripts
  proteomes/              not shipped; see "The nine reference proteomes"
```

## Artifacts

### `artifacts/contamination_census.json`

The tier assignment applied mechanically to all 72 leaderboard rows of both arms, with
per-tier counts, the measured containment figure behind each tier-2 assignment, the
`sha256` of every input table, and the rows where the stated rule and the flag shipped on
disk disagree. Backs the contamination census section.

```bash
python3 -I supplementary/scripts/contamination_census.py \
    --im-test <im_test.csv> \
    --infectious-table <TITANiAN_IM_Infectious.csv> \
    --neoantigen-table <TITANiAN_IM_Neoantigen.csv>
```

Requires the three upstream CSVs, which are public and are not redistributed:
`datasets/im_test.csv` from BigMHC's Mendeley deposit `10.17632/dvmz6pkzvb`, and the two
TITANiAN evaluation tables from `test.zip` in the Dryad deposit
`10.5061/dryad.s7h44j1k7`. The script compares each file's `sha256` against the one
behind the published figures and prints a warning on a mismatch rather than correcting
it.

### `artifacts/leaderboard_uncertainty_infectious.json`, `artifacts/leaderboard_uncertainty_neoantigen.json`

Per-entry within-allele AUROC with a 95% cluster-bootstrap interval, and every pairwise
contrast with its paired standard error, its minimum detectable effect and its
Holm-adjusted *p*, for the 9 infectious and 8 neoantigen entries whose per-row scores
ship. Back the resolution-limits section and the forest-plot figure. `cluster_basis`
reads `cluster` in both: these records are on the source-protein basis, not the weaker
peptide-overlap substitute.

```bash
python3 -I supplementary/scripts/leaderboard_uncertainty.py --set titanian_infectious --B 10000
python3 -I supplementary/scripts/leaderboard_uncertainty.py --set titanian_neoantigen  --B 10000
```

Requires the nine reference proteomes. The shipped rows supply everything else.

### `artifacts/unscored_arms_scored.json`

Scores and intervals for the four shipped arms that carry no predictor columns — 441
distinct rows across `hbv_hlabc_screen`, `xmuimmu_noncanonical`,
`independent_composite` and `itsndb` — with a minimum detectable effect on each contrast.
These are the only rows in the package where a published method's figure is tier 0 under
the census. The record states the predictor version used, MHCflurry `2.3.28` with model
set `models_class1_presentation 2.3.0`, which is not the `mhcflurry_2_0` column of the
two scored arms. Backs the resolution-limits section.

```bash
python3 -I supplementary/scripts/score_unscored_arms.py --B 10000
python3 -I supplementary/scripts/score_unscored_arms.py --B 10000 --no-mhcflurry   # baselines only
```

Requires the nine reference proteomes, and MHCflurry for its column only. The two
unfitted baselines need nothing beyond the package.

### `artifacts/stratum_choice.json`

Why allele is the stratification variable, measured rather than assumed: the five
candidate variables in the shipped schema, the three the estimator's degeneracy guards
refuse and the message each guard raised, the pooled-versus-within AUROC of a lookup
predictor built on each survivor, and the rank-order shift when every shipped score
column is re-ranked under each candidate's strata. Backs the stratum-choice subsection of
the estimator section.

```bash
python3 -I supplementary/scripts/stratum_choice.py
```

Needs only the shipped rows and a standard-library interpreter.

### `artifacts/rank_reversal.json`

How far pooled and within-allele AUROC disagree about the ranking, computed on the
published leaderboard CSVs so it measures exactly what a reader sees: per-arm maximum,
mean and median absolute rank shift, the count moving five and ten places or more, and
the five largest moves. Backs the rank-reversal subsection of the estimator section and
the reversal figure. Withheld figures are excluded rather than sorted as numbers.

```bash
python3 -I supplementary/scripts/rank_reversal.py [--emit-pgfplots DIR]
```

Needs only the two shipped leaderboard CSVs and a standard-library interpreter.

### `artifacts/collapse_rule.json`

Which duplicate-collapse rule each scored arm's labels are consistent with, measured
against the annotation record on the units where the arm and the export overlap. On the
infectious arm no deterministic rule fits, which is consistent with a single-source label
or a draw; on the neoantigen arm no unit is recorded as disagreeing, so the question is
not determinable there and the record says so. Backs the collapse-rule subsection of the
annotation section. **Not reproducible from the released package:** requires the IEDB
T-cell full export, and the code is not redistributed.

## Annotation artifacts — `artifacts/annotation/`

These seven records back the annotation-ceiling section. The analysis code is **not**
redistributed, because its entry point reads a 1.35 GB IEDB T-cell full export that
cannot be redistributed and the figures are specific to one dated version of it,
identified by `source_sha256` in each record. Shipping code that cannot be run from the
package would be worse than stating the dependency, so the artifacts ship and the code
does not.

### `artifacts/annotation/iedb_unit_votes.json` — the derived per-unit table

**53,754 units, 4.0 MB, and the reason the annotation analysis is checkable without the
source export.** One record per `(peptide, mhc_allele)` unit over the human-host, MHC
class I, plain-amino-acid, length 8–11 population, carrying the audited reference count
`k`, the positive, negative and tied reference counts, the per-reference and per-row
disagreement rates `p_ref` and `p_row`, and the row count. This is the table every
downstream step of the annotation chain consumes. A third party who holds it can
recompute the ceiling, the sensitivity analysis and the collapse-rule comparison without
touching the 1.35 GB export, and can check the published figures against it directly. The
`k >= 2` subpopulation, 4,815 units, is the audited headline.

### `artifacts/annotation/iedb_disagreement_ceiling.json`

The audited filter chain from the raw export to the analysed units, with every row count
along it and the two ceiling point estimates on both the reference-voting and row-voting
bases: Ceiling-1 `0.8237`, Ceiling-2 `0.9119`, from 4,815 units of 53,754 total.

### `artifacts/annotation/iedb_ceiling_interval.json`

The cluster bootstrap over units on both ceilings and on reference-level disagreement, at
B = 10,000: Ceiling-1 `0.8237` with 95% interval `[0.8163, 0.8309]`, Ceiling-2 `0.9119`
with `[0.9068, 0.9166]`.

### `artifacts/annotation/ceiling_sensitivity.json`

Sixteen variants over the subpopulation, the duplicate-collapse rule, the tie rule and
prevalence, each with its move against the published pair, plus a gate that recomputes
the baseline and fails on drift beyond `5e-05`. The record states which assumption is
load-bearing: the duplicate-collapse rule, and it is not a sensitivity but a change of
quantity, because positive-wins and majority make the test label a deterministic function
of the unit's annotation rate and leave no annotation-noise ceiling to report.

### `artifacts/annotation/ceiling_partial_identification.json`

What the two ceilings would be at each assumed disagreement rate among the 48,939 units
measured once, which no amount of the present data observes. Backs the headroom figure.
The bound is not vacuous anywhere on the range evaluated: at the observed split rate
`0.5084` the headroom to Ceiling-1 above the in-distribution cross-validated reference
`0.6125` is still `0.2120`.

### `artifacts/annotation/ceiling_biased_draw.json`

What the ceilings become when the duplicate draw is biased towards Positive rather than
uniform, with the logit shift solved to match the measured positive-rate excess. The
direction is the finding: the published ceilings are optimistic for an arm with this bias,
Ceiling-1 falling by `0.0812` to `0.7425`.

### `artifacts/annotation/disagreement_null.json`

The hierarchical null against which the rise of disagreement with reference count is
tested, on both the audited and the vote basis, at B = 10,000. Backs the
disagreement-and-consistency figure. Of the raw `0.0515` decline in mean per-unit
consistency, a single reference-count-independent latent distribution already predicts
`0.0428`, leaving `0.0087`: the decline is the estimator's own small-sample bias, not a
property of the literature.

## Scripts

Five scripts run against data this package ships. Each inserts the package root on
`sys.path` itself, so `python3 -I` works from any working directory and no `PYTHONPATH`
is needed. Each writes into `supplementary/artifacts/`.

| script | command | needs beyond the package |
|---|---|---|
| `stratum_choice.py` | `python3 -I supplementary/scripts/stratum_choice.py` | nothing |
| `rank_reversal.py` | `python3 -I supplementary/scripts/rank_reversal.py` | nothing |
| `leaderboard_uncertainty.py` | `python3 -I supplementary/scripts/leaderboard_uncertainty.py --set titanian_infectious --B 10000` | the nine proteomes |
| `score_unscored_arms.py` | `python3 -I supplementary/scripts/score_unscored_arms.py --B 10000` | the nine proteomes, MHCflurry |
| `contamination_census.py` | `python3 -I supplementary/scripts/contamination_census.py --im-test <f> --infectious-table <f> --neoantigen-table <f>` | the three upstream CSVs |

The annotation scripts are not among them. They require the IEDB export, as stated above.

### The nine reference proteomes

`leaderboard_uncertainty.py` and `score_unscored_arms.py` cluster each row to its source
protein by exact peptide-substring match against nine UniProt reviewed (Swiss-Prot)
proteomes. The package does not redistribute them: 69 MB of third-party FASTA under their
own terms. Place all nine in `supplementary/proteomes/`, or point `--proteomes` or
`$EPIGAUGE_PROTEOMES` at the directory holding them.

Each was fetched from `https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=<QUERY>`
on 2026-11-10. The sequence count and `sha256` identify the version the published figures
were computed on; both scripts carry the same table in `PROTEOME_INDEX` and print it when
a file is missing.

| file | UniProt query | sequences | sha256 |
|---|---|---|---|
| `escherichia_reviewed.fasta` | `(taxonomy_id:561) AND (reviewed:true)` | 23996 | `4a474daee45e65ed1d221a8a0f9f320c31aaa723266ce049dbd9d1104e28fab7` |
| `human_reviewed.fasta` | `(organism_id:9606) AND (reviewed:true)` | 20429 | `cee55a0d0af00b26876126c613abdff9d2e8649eb11a73e8db6c1badcd3d8fd5` |
| `mouse_reviewed.fasta` | `(organism_id:10090) AND (reviewed:true)` | 17285 | `1d590b004e502ec8e3e1edca352afef1e71e0fded5eb2519c7a29c46feda613f` |
| `mycobacterium_reviewed.fasta` | `(taxonomy_id:1763) AND (reviewed:true)` | 10141 | `10b9d8fe3b228bc0849b1d732c9d859b66c7e93debb550cdc7a5831cbd266773` |
| `plasmodium_reviewed.fasta` | `(taxonomy_id:5820) AND (reviewed:true)` | 774 | `4448ead0918a9482e518c7f9f76598f20c7d9efb5931a9b5fbff2e0559992afd` |
| `salmonella_reviewed.fasta` | `(taxonomy_id:590) AND (reviewed:true)` | 12890 | `1dad6b2e397f67bba0189e5faa4331bba2e05e23a9cfea7d40c89827ee1ed89b` |
| `staphylococcus_reviewed.fasta` | `(taxonomy_id:1279) AND (reviewed:true)` | 13512 | `039658207b5507c535de4f7e5dbe1d7cdbdbc22cc21272b81d247042bbeed7d8` |
| `streptococcus_reviewed.fasta` | `(taxonomy_id:1301) AND (reviewed:true)` | 12805 | `cea6d6e76555bd1a76a4623c8d509c7a2859c3eaf7c8cadd8b812ae410b13508` |
| `virus_reviewed.fasta` | `(taxonomy_id:10239) AND (reviewed:true)` | 17541 | `00cfa76a9039d541ea4f1021706e66f61e21d30f77a1f215cfe0110f5fdc6239` |

**A missing proteome is a hard failure, and a partial set is refused like an empty one.**
Both scripts exit non-zero and print every missing filename with its query and checksum.
Neither falls back. The reason is quantitative: assignment against the human proteome
alone is 4.6% on the infectious arm against 85.9% for all nine, so a missing file does not
make the basis slightly worse — it makes almost every row a singleton, which is the
row-level design wearing a clustered name.

`--clusters overlap` is available on both scripts and is **not** a fallback. It unions
rows by shared k-mer with each other rather than against reference proteomes, so two
non-overlapping peptides from one protein are not unioned; its intervals are
anti-conservative relative to the proteome basis. It exists only so the substitute used
before the proteomes were available can be reproduced. Its `cluster_basis` reads
`peptide-overlap`, and it writes to a separate `_overlap` filename, so it can neither be
confused with the shipped records in a table nor overwrite them on disk.

## What a third party regenerates from the released package alone

From this package and a standard-library interpreter, with no network access and no
third-party module:

- `stratum_choice.json`, from the shipped rows.
- `rank_reversal.json`, and with `--emit-pgfplots` the reversal figure's data tables, from
  the shipped leaderboard CSVs.

With the nine reference proteomes added, and nothing else:

- `leaderboard_uncertainty_infectious.json` and `leaderboard_uncertainty_neoantigen.json`,
  and therefore the forest-plot figure, for the 9 infectious and 8 neoantigen entries
  whose per-row scores ship.
- `unscored_arms_scored.json` for the two unfitted baselines on the four small arms.
  MHCflurry adds its own column.

## What needs an input that is not redistributed

- **The contamination census** needs the two TITANiAN evaluation tables and BigMHC's
  `im_test.csv`. Both deposits are public and are fetched rather than shipped; the DOIs
  and checksums are above.
- **The nine reference proteomes** are required by the two clustered-interval scripts.
  They are public UniProt downloads; the queries and checksums are above.
- **MHCflurry** is required for the MHCflurry column of `unscored_arms_scored.json` only.
- **The annotation figures** — the disagreement-and-consistency figure and the headroom
  figure — need the IEDB T-cell full export, 1.35 GB, identified by `source_sha256` in
  each annotation record, plus NumPy and pandas. The analysis code is not redistributed.
  `artifacts/annotation/iedb_unit_votes.json` is shipped so the derived per-unit table
  behind those figures is checkable without the export.

Two rows of each leaderboard are not recomputable by anyone from the released material,
the authors included. The NetMHCpan-4.1 columns are withheld under a licence clause and
their figures are not published. The three NORVARA columns are absent from the shipped
score files, so a third party can read those figures in the leaderboard tables but cannot
recompute them, or attach an interval to them, from this package.

## Licence

The artifacts and scripts in this directory are covered by the package licence, `LICENSE`
at the package root. The upstream deposits, the UniProt proteomes and the IEDB export
carry their own terms and are not redistributed here.
