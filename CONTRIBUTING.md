# Contributing to epigauge

Patches welcome. CI enforces every rule below. Exceptions go in the pull request
description and are reviewed as exceptions.

## The non-negotiables

1. **The metric layer is standard library only** — `canon.py`, `metrics.py` and anything they import must run on bare CPython 3.11+; numpy code goes behind the `fit` extra.
2. **Every new check ships a known-failing fixture** — no fixture, no guard; there is no waiver list.
3. **Absent required fields RAISE** — never default a missing field to a category (`CanonError`, `DegenerateGroupingError`); callers that tolerate bad rows must catch and record the drop.
4. **No tolerant `.get` on a required key** — use `row["allele"]`; same for `dict.setdefault`, `getattr(..., default)` and `except KeyError: pass`.
5. **Any printed zero must come from a counter a fixture proves can be non-zero** — otherwise a zero is indistinguishable from a counter that never increments.
6. **One canonicalisation function per key type** — `canonicalise_allele` / `_organism` / `_peptide`; no ad-hoc `.upper()`, `.replace()`, `.strip()` or regex on a key anywhere else.

## Before you open a pull request

```sh
# Fixtures, both entry points. Both must stay green.
python3 -I epigauge/tests/test_canon_and_metrics.py
python -m pytest epigauge/tests -q

# Every registered check has a fixture that demonstrates it firing.
python3 -I scripts/check_fixture_coverage.py

# Lint.
ruff check epigauge scripts
```

## Licensing of contributions

Code is MIT. Shipped data is CC BY 4.0 or CC0 per set, with the ITSNdb slice segregated under its
own notice — do not move data between those directories. No NC or ND sources: one NC row
file makes the whole bundle non-commercial. New data requires the upstream **LICENSE
body** read and quoted in the pull request.
