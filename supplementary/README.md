# `supplementary/`

Supplementary data for the epigauge preprint. Every artifact the paper's supplementary
index cites is here, together with the analysis code that a third party can run against
the data this package ships.

Paths in this directory and inside every artifact are relative to the package root — the
directory holding `pyproject.toml`, `leaderboards/` and `benchmarks/`. No artifact
records an absolute path. Where an artifact was computed from an input that is not
redistributed, its `script` or input field names the deposit, and every such artifact but
`annotation/ceiling_biased_draw.json` also carries the `sha256` of the exact file used,
which is how a third party confirms it holds the same input.

```
supplementary/
  artifacts/              17 JSON records, 4.2 MB
    annotation/           the 10 annotation-ceiling records, 4.0 MB
  scripts/                8 runnable analysis scripts
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

These ten records back the annotation-ceiling section, and they divide into two groups.

**Seven have no shipped code.** `iedb_unit_votes.json`, `iedb_disagreement_ceiling.json`,
`iedb_ceiling_interval.json`, `ceiling_sensitivity.json`,
`ceiling_partial_identification.json`, `ceiling_biased_draw.json` and
`disagreement_null.json` are produced by entry points that read a 1.35 GB IEDB T-cell full
export, which cannot be redistributed, and their figures are specific to one dated version
of it. Six of the seven identify it by `source_sha256`; `ceiling_biased_draw.json` records
its basis, its calibration and the solved logit shift instead. Shipping code that cannot be
run from the package would be worse than stating the dependency, so the artifacts ship and
the code does not. Each of the seven names the missing entry point in its own `script`
field.

**Three ship with their code and run on a bare interpreter.**
`ceiling_crossfit.json`, `ceiling_convergence.json` and `ceiling_deconvolved_refused.json`
are produced by `supplementary/scripts/ceiling_crossfit.py`,
`ceiling_convergence.py` and `ceiling_deconvolved.py`. All three read
`artifacts/annotation/iedb_unit_votes.json`, which ships, and need nothing else, so the
bracket on Ceiling-2 and the refusal of the deconvolution are regenerable from this
package alone.

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
along it and the two plug-in ceiling figures on both the reference-voting and row-voting
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

What the two ceilings would be at each assumed disagreement rate δ among the 48,939 units
measured once, which no amount of the present data observes. Backs the headroom figure.
Every row is an upper bound **at its own δ** and the record supports nothing stronger: the
curve decreases in δ, from Ceiling-1 `0.9822` and Ceiling-2 `0.9991` at δ = 0 to `0.8245`
and `0.9119` at the observed reference-level rate `0.5084`. No global bound is stated from
it. The headroom the paper reports is the subpopulation statement in
`ceiling_convergence.json`, on Ceiling-2.

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

### `artifacts/annotation/ceiling_crossfit.json`

The lower bound of **Ceiling-2's** bracket, and of that one only. Both published ceiling
figures are plug-ins: they rank by `p_ref` and weight by the same `p_ref`, so they inherit
its estimation error, and most audited units carry only two references. This splits each
unit's reference votes in half, scores from one half and takes the target from the other,
so no vote supplies both a score and its own target. Over 400 splits Ceiling-2 is
**0.6230** against the plug-in 0.9119. Reads `iedb_unit_votes.json` and nothing else.

```bash
python3 -I supplementary/scripts/ceiling_crossfit.py          # --splits defaults to 400
```

The bound holds for Ceiling-2 term by term: the half-B label masses are unbiased for `p`
and independent of the half-A scores, so each unordered pair's cross-fit contribution is at
most its true contribution, with equality when the two units share a rate.

**Ceiling-1 is not cross-fit here, and it has no lower bound anywhere in this package.**
Its score is the majority label, and `maj()` on a single vote in {0, ½, 1} is the identity,
so on a two-vote split the Ceiling-1 and Ceiling-2 cross-fit scorers are the same function
and 62.9% of audited units carry two references. No proof bounds Ceiling-1 below, and there
is a counterexample class: for two units with `p_u > p_v` sharing a majority label the
cross-fit contribution strictly exceeds the true one. Ceiling-1 is reported as `0.8237`, a
plug-in **upper** bound on the score of one particular predictor, which is not a bound on
model performance. `ceiling_convergence.json` records the analogous half-vote statistic as
`c1_crossfit_half_vote`, flagged `c1_crossfit_half_vote_is_a_bound: false`.

Note the definition the table actually stores: `p_ref = (n_pos_refs + n_tie_refs / 2) / k`.
It is neither `n_pos_refs / k` (which gives 0.9080) nor `n_pos_refs / (n_pos_refs +
n_neg_refs)` (0.9523). The script refuses to report unless it first reproduces the
published 0.9119 from the shipped table.

### `artifacts/annotation/ceiling_convergence.json`

Both bounds on **Ceiling-2**, and the plug-in upper bound on Ceiling-1, at seven minimum
reference counts: 2, 3, 4, 5, 6, 8 and 10. All seven strata are recorded. For Ceiling-2 the
plug-in bound is biased up and the cross-fit bound down, so the ceiling lies between them,
and both biases shrink as a unit accumulates references. The bracket narrows **overall**,
from 0.2889 at `k >= 2` to 0.1249 at `k >= 10`, and it does so non-monotonically: the width
rises across the first step (0.2889 against 0.2900 at `k >= 3`), the plug-in bound rises at
`k` 8 to 10 (0.7904 against 0.7927), and the cross-fit bound falls across the first step
(0.6230 against 0.5450). The record states each exception in its own
`c2_width_rises_at`, `c2_plugin_rises_at` and `c2_crossfit_falls_at` fields rather than
leaving a monotonicity claim to be checked by hand, and the `interpretation` string is
derived from the rows so it cannot drift from them.

**Headroom is reported on Ceiling-2 only**, because Ceiling-1 has no lower bound: above the
in-distribution cross-validated 0.6125 it is at most **0.2994**, with **0.0105** not
excluded. `c1_crossfit_half_vote` is the half-vote statistic for Ceiling-1 and is not a
bound; `c1_lower_bound` is `null` and `c1_lower_bound_note` says why.

The level at which Ceiling-2's two bounds approach each other on deeply measured units is a
property of those more-studied epitopes and does not transfer to the audited subpopulation.
The narrowing is the evidence; the level is not.

```bash
python3 -I supplementary/scripts/ceiling_convergence.py       # --splits defaults to 400
```

The record states no variance decomposition of `p_ref` and the bracket does not rest on
one. `population_shares` carries only counts: 62.9% of audited units hold two references
and 61.3% of those sit exactly at an endpoint.

### `artifacts/annotation/ceiling_deconvolved_refused.json` — a model that was rejected

A beta-binomial deconvolution of the latent disagreement rate, kept because a reader should
know it was tried. It **exits 1**, and records the refusal rather than only printing it:
the fitted parameters, the per-stratum variance comparison, both gate verdicts and the
ceilings it would have reported. The observed variance of `p_ref` falls
faster with reference count than any beta-binomial permits (at `k = 2`, 0.145 observed
against 0.116 predicted; at `k = 4`, 0.042 against 0.080), so its own goodness-of-fit gate
refuses it. Its answer of 0.8016 is held as a rejected fit, not a result, and the record's
own `note` field says so. No figure in the paper is computed from it. The bracket in
`ceiling_convergence.json` assumes no latent family and is the reportable version.

## Scripts

Eight scripts ship. Each inserts the package root on `sys.path` itself, so `python3 -I`
works from any working directory and no `PYTHONPATH` is needed. Each writes into
`supplementary/artifacts/`. Five run on the shipped data alone; three need an input this
package does not redistribute, named in the last column.

| script | command | needs beyond the package |
|---|---|---|
| `stratum_choice.py` | `python3 -I supplementary/scripts/stratum_choice.py` | nothing |
| `rank_reversal.py` | `python3 -I supplementary/scripts/rank_reversal.py` | nothing |
| `ceiling_crossfit.py` | `python3 -I supplementary/scripts/ceiling_crossfit.py` | nothing |
| `ceiling_convergence.py` | `python3 -I supplementary/scripts/ceiling_convergence.py` | nothing |
| `ceiling_deconvolved.py` | `python3 -I supplementary/scripts/ceiling_deconvolved.py` | nothing; **exits 1 by design** |
| `leaderboard_uncertainty.py` | `python3 -I supplementary/scripts/leaderboard_uncertainty.py --set titanian_infectious --B 10000` | the nine proteomes |
| `score_unscored_arms.py` | `python3 -I supplementary/scripts/score_unscored_arms.py --B 10000` | the nine proteomes, and MHCflurry for its own column |
| `contamination_census.py` | `python3 -I supplementary/scripts/contamination_census.py --im-test <f> --infectious-table <f> --neoantigen-table <f>` | the three upstream CSVs |

**The three ceiling scripts need no IEDB export.** They read the shipped per-unit vote
table, `artifacts/annotation/iedb_unit_votes.json`, and nothing else, and each refuses to
write unless it first reproduces the published plug-in figures from that table. The seven
annotation records listed above as having no shipped code are the ones whose entry points
read the export; those seven are not regenerable here.

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
- `rank_reversal.json`, and with `--emit-pgfplots DIR` the reversal figure's data tables,
  from the shipped leaderboard CSVs.
- `annotation/ceiling_crossfit.json` and `annotation/ceiling_convergence.json`, and
  therefore the bracket on Ceiling-2, from the shipped per-unit vote table.
- `annotation/ceiling_deconvolved_refused.json`, from the same table; the script exits 1
  and records the refusal.

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
- **Seven of the ten annotation records** — the per-unit vote table itself, the filter
  chain, the bootstrap interval, the sensitivity analysis, the partial-identification
  curve, the biased draw and the disagreement null — need the IEDB T-cell full export,
  1.35 GB, plus NumPy and pandas. Six identify the export by `source_sha256`. That code is
  not redistributed, so those seven ship as artifacts only, and
  `artifacts/annotation/iedb_unit_votes.json` is shipped so the derived per-unit table
  behind them is checkable without the export. The remaining three annotation records,
  `ceiling_crossfit.json`, `ceiling_convergence.json` and
  `ceiling_deconvolved_refused.json`, need nothing beyond this package.

Two rows of each leaderboard are not recomputable by anyone from the released material,
the authors included. The NetMHCpan-4.1 columns are withheld under a licence clause and
their figures are not published. The three NORVARA columns are absent from the shipped
score files, so a third party can read those figures in the leaderboard tables but cannot
recompute them, or attach an interval to them, from this package.

## Licence

The artifacts and scripts in this directory are covered by the package licence, `LICENSE`
at the package root. The upstream deposits, the UniProt proteomes and the IEDB export
carry their own terms and are not redistributed here.
