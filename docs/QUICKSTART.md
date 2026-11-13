# Quickstart

Score a predictor against the benchmark in about two minutes. Python 3.11+, no dependencies.

## Install

```bash
git clone https://github.com/norvara-therapeutics/epigauge
cd epigauge
pip install -e .
```

The benchmark rows ship with the package — nothing to download. `epigauge` is now on PATH and
works from any directory.

## See what is available

```console
$ epigauge list
titanian_infectious    2093   1488 (71.1%)  CC0 1.0   ships           TITANiAN infectious ...
titanian_neoantigen     833    178 (21.4%)  CC0 1.0   ships           TITANiAN neoantigen ...
titanian_ada            216        unknown  CC0 1.0   fetch required  TITANiAN antibody (ADA) ...
```

`ships` means on disk and sha256-verified. `fetch required` means the licence prevents
redistribution, so `epigauge fetch <set>` pulls it from the upstream.

## Score a predictor

Your predictions are a two-column CSV — one row per benchmark row, keyed by `row_id`:

```csv
row_id,score
titanian:infectious:000001,0.3333
titanian:infectious:000002,0.4444
```

Get the `row_id` values from the set itself:

```python
from epigauge import data
rows = data.load("titanian_infectious")
ids = [r["record_id"] for r in rows]          # 2,093 of them
```

Then score, declaring which direction your scores run:

```console
$ epigauge score my_predictions.csv --set titanian_infectious --direction higher-is-better

SCORECARD  my_predictions.csv
titanian_infectious  n=2093  1488 positive / 605 negative (71.1%)  78 strata
minimum detectable difference: 0.0950
estimator: pair-weighted within-allele AUROC over admissible strata only.
  Pooled is printed beside it and is never ranked on.

rank  within-allele     pooled  strata  evaluable rows     method
----  -------------  ---------  ------  -----------------  --------------------
   1         0.5264     0.5650   36/78  1787/2093 (85.4%)  BASELINE Hydropathy
   =         0.5209     0.5567   36/78  1787/2093 (85.4%)  my_predictions
      NOT RESOLVABLE against the row above: gap 0.0055 < MDE 0.0950
   3         0.4132     0.4512   36/78  1787/2093 (85.4%)  BASELINE Length
```

Read it as: your predictor landed between the hydropathy and length baselines, and the gap to
hydropathy is **smaller than this benchmark can resolve**, so the `=` is a tie rather than a rank.

## The rest of the surface

```bash
epigauge leaderboard --set titanian_infectious   # every shipped reference column
epigauge validate my_predictions.csv             # format check, no scoring
epigauge fixtures                                # 102 regression fixtures
```

Exit codes: `0` clean, `1` a check failed, `2` could not run. A `2` is never reported as a `0`.

Next: [METRICS](METRICS.md) for the estimator and its guards, [SUBMITTING](SUBMITTING.md) to put a
method on the public leaderboard.
