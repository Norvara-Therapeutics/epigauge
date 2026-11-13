# Submitting

**You submit by opening a pull request on this repository.** Add one file —
`submissions/<team-slug>/<YYYY-MM-DD>-<method>.yaml` — and CI validates it, scores it, and posts
the scorecard as a comment on your PR.

Predictions are **not** committed; they are large. Publish the CSV as a GitHub release asset or a Zenodo
record and link it from the YAML, pinned by sha256. Merging appends a row to `submissions/ledger.csv` — the
commit *is* the ledger entry, so the count of evaluations ever run here is public and auditable.

## The predictions CSV

```csv
row_id,score
titanian:infectious:000001,0.8214
titanian:infectious:000002,0.1077
```

## The submission file

Flat YAML: `key: value`, `  - item` lists, whole-line comments — no nesting, no inline
comments. Copy [the example](../submissions/example/2026-11-13-hydropathy-demo.yaml), which CI re-validates on every PR.

```yaml
schema_version: 1
team: example
contact: research@norvaratx.com
method: hydropathy-demo
version: 0.1.0
benchmark: titanian_infectious
predictions_url: https://github.com/norvara-therapeutics/epigauge/releases/download/example-submission/hydropathy-demo-0.1.0.csv
predictions_sha256: c19cadb2655a430dee93b3458dd906c664daf7728c662af91d838e5fd54f36ee
score_direction: higher-is-better
uses_peptide_as_input: true
training_cutoff: 2026-03-01
training_overlap: public
training_hashes_url: https://github.com/norvara-therapeutics/epigauge/releases/download/example-submission/hydropathy-demo-training-hashes.txt
training_hashes_sha256: 6826c010b26bf1de09a1449b88008771d90e58a2666dc8e858c406a334cb406a
training_sources:
  - IEDB 2026-03 class-I human T-cell assays
```

## Validation rules

Run the same code CI runs, before you open the PR:

```bash
python3 -I scripts/score_submission.py submissions/<team>/<file>.yaml --schema-only
epigauge validate my_predictions.csv --set infectious --direction higher-is-better
```

- One score per **`row_id`**, never per `(peptide, allele)` — that pair recurs across studies with different labels, so it is not a key.
- Every benchmark `row_id` present exactly once, every score finite. A partial submission is rejected, not scored on a subset.
- **A submission with no declared `score_direction` is rejected** — it is never inferred, because inferring it fits a free parameter on the evaluation labels.
- Training overlap is declared as sorted sha256 hashes of your `epigauge.canon`-canonicalised training `(peptide, HLA)` pairs, so neither side reveals data. `none` is accepted and prints `overlap: undeclared`; a high overlap is a column, not a disqualification.
- `uses_peptide_as_input: false` is a **legal** submission: scored, ranked and published like any other.
- `contact` is a GitHub handle, an https URL or an email; `predictions_url` is https on a GitHub-release or Zenodo host, and its sha256 must match the bytes CI downloads.
- The PR adds exactly one file, whose path matches the `team`, `method` and date declared inside it. A ledger row is never deleted or amended; a withdrawal is a new row.

## What the scorecard returns

The pair-weighted within-allele AUROC as the ranked figure, carrying its admissible-stratum
count and evaluable-row fraction; pooled beside it and never ranked on; the minimum detectable
difference at the top, with any gap below it printed as `=` and the measured gap instead of a
rank; the baseline roster as rows, not footnotes; and your training-overlap percentage.
Exit codes: `0` clean, `1` a check failed, `2` could not run.

Why the estimator is stratified, and why pooled and stratified ship together or neither ships:
**[METRICS.md](METRICS.md)**. Scorecard anatomy: [QUICKSTART.md](QUICKSTART.md). Overlap mechanics: [CONTAMINATION.md](CONTAMINATION.md).
