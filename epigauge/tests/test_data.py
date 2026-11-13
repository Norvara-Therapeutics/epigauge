"""Fixtures for the shipped benchmark corpus and its licence boundary.

Governing rule in `test_canon_and_metrics.py`. Each test reproduces a way the data layer
could silently publish a wrong number, or someone else's data. Standalone:
`python3 -I epigauge/tests/test_data.py`.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from epigauge import data  # noqa: E402
from epigauge.canon import (  # noqa: E402
    CanonError,
    canonicalise_allele,
    canonicalise_peptide,
)
from epigauge.data import (  # noqa: E402
    DataError,
    IntegrityError,
    LicenceError,
    UnknownSetError,
    available,
    corpus_fingerprint,
    describe,
    fetch,
    load,
    manifest,
    score_columns,
)

FAILURES: list[str] = []
_STANDALONE = False


def check(cond: bool, msg: str) -> None:
    if not cond:
        FAILURES.append(msg)


def expect_raises(exc, fn, *a, **kw):
    """Return the exception instance, or None (and record a failure)."""
    try:
        fn(*a, **kw)
    except exc as e:
        return e
    except Exception as e:  # noqa: BLE001
        FAILURES.append(f"expected {exc.__name__}, got {type(e).__name__}: {e}")
        return None
    FAILURES.append(f"expected {exc.__name__}, nothing raised: {getattr(fn, '__name__', fn)}")
    return None


def _shipped() -> dict[str, dict]:
    return {k: v for k, v in manifest()["sets"].items() if v["status"] == "ship"}


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _Sandbox:
    """A throwaway copy of `benchmarks/` that a test may deliberately corrupt; the real
    tree is never written to."""

    def __init__(self) -> None:
        self.tmp = Path(tempfile.mkdtemp(prefix="epigauge-data-test-"))
        self.root = self.tmp / "benchmarks"
        shutil.copytree(data.benchmarks_dir(), self.root)
        self._prev = os.environ.get("EPIGAUGE_BENCHMARKS")

    def mutate(self, fn) -> None:
        path = self.root / "MANIFEST.json"
        m = json.loads(path.read_text(encoding="utf-8"))
        fn(m)
        path.write_text(json.dumps(m, indent=1), encoding="utf-8")

    def __enter__(self) -> "_Sandbox":
        os.environ["EPIGAUGE_BENCHMARKS"] = str(self.root)
        data._MANIFEST_CACHE.clear()
        return self

    def __exit__(self, *exc) -> None:
        if self._prev is None:
            os.environ.pop("EPIGAUGE_BENCHMARKS", None)
        else:
            os.environ["EPIGAUGE_BENCHMARKS"] = self._prev
        data._MANIFEST_CACHE.clear()
        shutil.rmtree(self.tmp, ignore_errors=True)


def test_manifest_sha256_matches_every_shipped_file() -> None:
    """A file that drifts from the manifest and loads anyway redefines the benchmark
    under a stable name; fingerprint b2c046cd6d9e72e5 is unreproducible for that reason."""
    root = data.benchmarks_dir()
    shipped = _shipped()
    check(bool(shipped), "no shipped sets in MANIFEST.json")
    for set_id, e in shipped.items():
        p = root / e["path"]
        check(p.is_file(), f"{set_id}: shipped file missing: {p}")
        if not p.is_file():
            continue
        check(_sha256(p) == e["sha256"], f"{set_id}: sha256 does not match MANIFEST.json")
        check(p.stat().st_size == e["bytes"],
              f"{set_id}: byte count {p.stat().st_size} != manifest {e['bytes']}")
        sc = e.get("scores")
        if sc:
            sp = root / sc["path"]
            check(sp.is_file(), f"{set_id}: score sidecar missing: {sp}")
            if sp.is_file():
                check(_sha256(sp) == sc["sha256"],
                      f"{set_id}: score sidecar sha256 does not match MANIFEST.json")


def test_corpus_fingerprint_recomputes_from_the_files_on_disk() -> None:
    check(corpus_fingerprint(verify=True) == manifest()["corpus_fingerprint"],
          "corpus fingerprint does not recompute from the shipped files")


def test_every_shipped_row_canonicalises_without_raising() -> None:
    """One file spelling HLA-A*01:01 and another A0101 returns 0 of 833 on the join and
    hides a contaminated benchmark, so no non-canonical row may enter."""
    for set_id, e in _shipped().items():
        rows = load(set_id)
        for r in rows:
            try:
                pep = canonicalise_peptide(r["peptide"])
                al = canonicalise_allele(r["allele"])
            except CanonError as exc:
                FAILURES.append(f"{set_id}: row {r.get('record_id')} will not "
                                f"canonicalise: {exc}")
                break
            if pep != r["peptide"]:
                FAILURES.append(f"{set_id}: {r['record_id']} peptide {r['peptide']!r} "
                                f"is not already canonical ({pep!r})")
                break
            if al != r["allele"]:
                FAILURES.append(f"{set_id}: {r['record_id']} allele {r['allele']!r} "
                                f"is not already canonical ({al!r})")
                break
            if r["label"] not in (0, 1):
                FAILURES.append(f"{set_id}: {r['record_id']} label {r['label']!r} "
                                "is not 0/1")
                break
            if r["length"] != len(pep):
                FAILURES.append(f"{set_id}: {r['record_id']} length field disagrees "
                                "with the peptide")
                break


def test_row_counts_and_class_balance_match_the_manifest() -> None:
    """A truncated file is a different benchmark."""
    for set_id, e in _shipped().items():
        rows = load(set_id)
        check(len(rows) == e["rows"],
              f"{set_id}: loaded {len(rows)} rows, manifest declares {e['rows']}")
        pos = sum(r["label"] for r in rows)
        check(pos == e["positives"],
              f"{set_id}: {pos} positives, manifest declares {e['positives']}")
        check(len(rows) - pos == e["negatives"],
              f"{set_id}: {len(rows) - pos} negatives, manifest declares {e['negatives']}")
        ids = [r["record_id"] for r in rows]
        check(len(set(ids)) == len(ids), f"{set_id}: record_id is not unique")


# Corruption that loads anyway, from four sides: wrong digest, wrong bytes, wrong row
# count, wrong sidecar. There is no --force path.

def test_corrupted_manifest_sha256_raises_rather_than_loading() -> None:
    with _Sandbox() as box:
        box.mutate(lambda m: m["sets"]["itsndb"].__setitem__("sha256", "0" * 64))
        e = expect_raises(IntegrityError, load, "itsndb")
        if e:
            check("sha256" in str(e) and "0" * 64 in str(e),
                  "IntegrityError must print both the expected and the actual digest")


def test_corrupted_data_file_raises_rather_than_loading() -> None:
    with _Sandbox() as box:
        p = box.root / manifest()["sets"]["itsndb"]["path"]
        with p.open("a", encoding="utf-8") as fh:
            fh.write('{"record_id":"itsndb:SMUGGLED","peptide":"SIINFEKL",'
                     '"allele":"HLA-A*02:01","label":1}\n')
        expect_raises(IntegrityError, load, "itsndb")


def test_wrong_row_count_in_manifest_raises() -> None:
    with _Sandbox() as box:
        # Digest kept consistent, so the count check is shown firing on its own.
        box.mutate(lambda m: m["sets"]["itsndb"].__setitem__("rows", 198))
        expect_raises(IntegrityError, load, "itsndb")


def test_corrupted_score_sidecar_raises() -> None:
    with _Sandbox() as box:
        box.mutate(lambda m: m["sets"]["titanian_neoantigen"]["scores"]
                   .__setitem__("sha256", "f" * 64))
        expect_raises(IntegrityError, load, "titanian_neoantigen", scores=True)


def test_tampered_fingerprint_is_caught_by_verification() -> None:
    with _Sandbox() as box:
        box.mutate(lambda m: m.__setitem__("corpus_fingerprint", "deadbeefdeadbeef"))
        expect_raises(IntegrityError, corpus_fingerprint, verify=True)


_MUST_REFUSE = {
    # set_id: a phrase from the operative clause that the refusal must carry
    "vdjdb": "ND",
    "mcpas_tcr": "licence",
    "neoranking_muller": "LICENSE",
    "hitide_ega": "controlled",
    "tesla_synapse": "redistribute",
    "tchard": "Non-Commercial",
    "tscape_titanian": "NC",
    "tesla608": "REDISTRIBUTION FORBIDDEN",
    "buckley_gbm": "404",
    "titanian_mhc1": "RETRIEVAL NOT RECORDED",
    "titanian_ada": "OUT OF SCOPE",
}


def test_non_shippable_sets_raise_with_the_licence_reason() -> None:
    """`load` on a set we may not redistribute must raise and name the clause."""
    for set_id, phrase in _MUST_REFUSE.items():
        e = expect_raises(LicenceError, load, set_id)
        if e is None:
            continue
        msg = str(e)
        check(phrase in msg,
              f"{set_id}: refusal must quote the licence basis; "
              f"{phrase!r} not in {msg!r}")
        check(set_id in msg, f"{set_id}: refusal must name the set")


def test_load_never_returns_empty_for_a_withheld_set() -> None:
    """An empty list is indistinguishable from a benchmark with no rows, which is how a
    licence breach becomes a published metric."""
    for set_id in _MUST_REFUSE:
        try:
            rows = load(set_id)
        except LicenceError:
            continue
        except Exception as exc:  # noqa: BLE001
            FAILURES.append(f"{set_id}: expected LicenceError, got {type(exc).__name__}")
            continue
        FAILURES.append(f"{set_id}: load() returned {len(rows)} rows for a set that "
                        "may not be redistributed")


def test_fetch_refuses_the_sets_that_cannot_be_fetched_automatically() -> None:
    # Controlled access or no licence at all: no helper, under any flag.
    for set_id in ("vdjdb", "mcpas_tcr", "neoranking_muller", "hitide_ega",
                   "tesla_synapse", "tchard", "tscape_titanian"):
        e = expect_raises(LicenceError, fetch, set_id)
        if e:
            check("No fetch helper is offered" in str(e),
                  f"{set_id}: refusal must say that no fetch helper exists")
    # A DOI landing page is not a file: refuse rather than guess a download URL.
    for set_id in ("titanian_mhc1", "titanian_mhc2", "titanian_ada"):
        e = expect_raises(LicenceError, fetch, set_id)
        if e:
            check("sha256" in str(e),
                  f"{set_id}: manual-fetch instructions must pin the sha256")


def test_fetch_refuses_a_shipped_set_rather_than_making_a_second_copy() -> None:
    e = expect_raises(LicenceError, fetch, "itsndb")
    if e:
        check("load()" in str(e), "refusal should redirect to load()")


def test_fetchable_sets_pin_an_upstream_sha256() -> None:
    """Without a pinned digest a changed upstream silently redefines the benchmark."""
    for set_id, e in manifest()["sets"].items():
        if e["status"] != "fetch":
            continue
        check(bool(e["source"].get("upstream_sha256")),
              f"{set_id}: fetchable set has no pinned upstream sha256")
        check(bool(e["licence_basis"]),
              f"{set_id}: fetchable set has no recorded licence basis")


def test_unknown_set_id_raises_and_lists_what_exists() -> None:
    e = expect_raises(UnknownSetError, load, "netmhcpan_benchmark")
    if e:
        check("itsndb" in str(e), "unknown-set error should list the known sets")


def test_no_netmhcpan_column_ships_anywhere() -> None:
    """DTU 7(v) bars publishing NetMHCpan benchmark results without written consent, and
    no waiver exists."""
    cols = manifest()["predictor_columns"]
    for name, spec in cols.items():
        if "netmhc" in name:
            check(spec["verdict"] == "exclude",
                  f"{name}: NetMHC* columns may not ship (DTU 7(v))")
    for set_id, e in _shipped().items():
        sc = e.get("scores")
        if not sc:
            continue
        for c in sc["columns_shipped"]:
            check("netmhc" not in c, f"{set_id}: ships NetMHC column {c}")
        rows = load(set_id, scores=True)
        for r in rows[:50]:
            for c in r["scores"]:
                check("netmhc" not in c,
                      f"{set_id}: row {r['record_id']} carries NetMHC column {c}")


def test_every_score_column_present_on_disk_is_cleared_in_the_manifest() -> None:
    cols = manifest()["predictor_columns"]
    for set_id, e in _shipped().items():
        if not e.get("scores"):
            continue
        rows = load(set_id, scores=True)
        seen: set[str] = set()
        for r in rows:
            seen |= set(r["scores"])
        for c in sorted(seen):
            spec = cols.get(c)
            check(spec is not None, f"{set_id}: column {c} has no manifest verdict")
            if spec:
                check(spec["verdict"] == "ship",
                      f"{set_id}: column {c} is on disk but its verdict is "
                      f"{spec['verdict']!r}")
        check(seen == set(e["scores"]["columns_shipped"]),
              f"{set_id}: columns on disk differ from columns_shipped")


def test_withheld_columns_are_absent_not_null() -> None:
    rows = load("titanian_infectious", scores=True)
    withheld = set(manifest()["sets"]["titanian_infectious"]["scores"]["columns_withheld"])
    check(bool(withheld), "titanian_infectious should withhold columns")
    for r in rows[:200]:
        leaked = withheld & set(r["scores"])
        check(not leaked,
              f"{r['record_id']}: withheld columns present as keys: {sorted(leaked)}")


def test_score_columns_explains_every_withholding() -> None:
    sc = score_columns("titanian_neoantigen")
    check(bool(sc["withheld"]), "expected withheld columns for titanian_neoantigen")
    for name, spec in sc["withheld"].items():
        check(spec["verdict"] == "exclude", f"{name}: withheld but verdict is not exclude")
        check(len(spec["reason"]) > 40, f"{name}: withholding reason is not a reason")
        check(bool(spec["licence"]), f"{name}: withheld with no licence recorded")
    check(score_columns("itsndb")["shipped"] == {},
          "itsndb ships no score sidecar")


def test_every_set_records_a_licence_and_a_basis() -> None:
    """A licence asserted without a basis is not a licence."""
    for set_id, e in manifest()["sets"].items():
        check(bool(e.get("licence")), f"{set_id}: no licence recorded")
        check(len(e.get("licence_basis", "")) > 40,
              f"{set_id}: licence basis is missing or too thin to be a basis")
        check(e["status"] in ("ship", "fetch", "excluded"),
              f"{set_id}: unknown status {e['status']!r}")


def test_itsndb_is_segregated_and_makes_no_redistribution_claim() -> None:
    """The README claimed "GPL-3.0, redistribution permitted"; the LICENSE is boilerplate
    with no copyright holder and covers the R package, not the records."""
    e = describe("itsndb")
    check(e["path"].startswith("itsndb/"),
          "ITSNdb must live in its own directory, not in the CC BY 4.0 tree")
    notice = data.benchmarks_dir() / "itsndb" / "NOTICE.md"
    check(notice.is_file(), "ITSNdb slice must carry its own NOTICE.md")
    if notice.is_file():
        text = notice.read_text(encoding="utf-8")
        check("redistribution permitted" not in text,
              "the softened ITSNdb claim must not reassert 'redistribution permitted'")
        check("no separate licence" in text,
              "the ITSNdb notice must state that the records carry no separate licence")
    check("CC BY" not in e["licence"],
          "ITSNdb must not be labelled CC BY; its position is unresolved")


def test_shipped_tier_stays_small() -> None:
    """`shipped_bytes` counts row and score files ONLY: including `MANIFEST.json` made
    the total self-referential, so editing the manifest changed what it asserted."""
    root = data.benchmarks_dir()
    m = manifest()

    counted = 0
    for entry in m["sets"].values():
        for block in (entry, entry.get("scores") or {}):
            rel = block.get("path")
            if not rel:
                continue
            f = root / rel
            if not f.is_file():
                f = root / Path(rel).name
            if f.is_file():
                counted += f.stat().st_size

    check(counted == m["shipped_bytes"],
          f"manifest shipped_bytes {m['shipped_bytes']} != row+score total {counted}")
    check("EXCLUDES this manifest" in m.get("shipped_bytes_definition", ""),
          "shipped_bytes_definition must state that the manifest is excluded, or the "
          "total silently becomes self-referential again")

    on_disk = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    check(on_disk < 10 * 1024 * 1024,
          f"shipped tier is {on_disk / 1e6:.2f} MB on disk, over the 10 MB target")


def test_available_lists_withheld_sets_rather_than_hiding_them() -> None:
    all_ids = {s.set_id for s in available()}
    check(_MUST_REFUSE.keys() <= all_ids,
          "available() must list the sets you cannot have, with their licence")
    shipped = {s.set_id for s in available(status="ship")}
    check(shipped == set(_shipped()), "available(status='ship') disagrees with manifest")
    for s in available():
        check(s.shipped == (s.status == "ship"), f"{s.set_id}: SetInfo.shipped is wrong")


def test_missing_benchmarks_directory_raises_rather_than_loading_nothing() -> None:
    prev = os.environ.get("EPIGAUGE_BENCHMARKS")
    tmp = tempfile.mkdtemp(prefix="epigauge-empty-")
    try:
        os.environ["EPIGAUGE_BENCHMARKS"] = tmp
        data._MANIFEST_CACHE.clear()
        # Package-relative fallbacks still resolve in a source checkout, so this only
        # asserts that an empty EPIGAUGE_BENCHMARKS is never silently accepted.
        d = data.benchmarks_dir()
        check(Path(d) != Path(tmp), "an empty directory was accepted as benchmarks/")
    except DataError:
        pass
    finally:
        if prev is None:
            os.environ.pop("EPIGAUGE_BENCHMARKS", None)
        else:
            os.environ["EPIGAUGE_BENCHMARKS"] = prev
        data._MANIFEST_CACHE.clear()
        shutil.rmtree(tmp, ignore_errors=True)


def test_zz_every_fixture_above_reported_no_failures() -> None:
    """Turns failures recorded by `check()` into a PYTEST failure; without it pytest
    stays green whatever was recorded. Last in the file so it runs last."""
    if _STANDALONE:
        return  # main() prints and returns the same information
    assert not FAILURES, "\n".join(["recorded failures:"] + [f"  - {f}" for f in FAILURES])


def main() -> int:
    global _STANDALONE
    _STANDALONE = True
    tests = [v for k, v in sorted(globals().items())
             if k.startswith("test_") and callable(v)]
    for t in tests:
        before = len(FAILURES)
        try:
            t()
        except Exception as e:  # noqa: BLE001
            FAILURES.append(f"{t.__name__} raised {type(e).__name__}: {e}")
        status = "ok" if len(FAILURES) == before else "FAIL"
        print(f"  [{status:>4}] {t.__name__}")
    print(f"\n{len(tests)} fixtures, {len(FAILURES)} failure(s)")
    for f in FAILURES:
        print(f"  - {f}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
