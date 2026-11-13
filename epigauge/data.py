"""Benchmark row sets: what ships, what you fetch, and what nobody may ship.

Every path verifies the pinned sha256 before yielding rows or returning one; there is no
--force. A set we may not redistribute raises `LicenceError` with the operative clause,
never `[]`. Withheld predictor columns are absent from the sidecars, not null (NetMHCpan:
DTU 7(v), no consent). `MANIFEST.json` is the authority; this module enforces it.
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "DataError",
    "LicenceError",
    "IntegrityError",
    "UnknownSetError",
    "SetInfo",
    "available",
    "describe",
    "load",
    "fetch",
    "manifest",
    "benchmarks_dir",
    "score_columns",
    "corpus_fingerprint",
]

_CHUNK = 1 << 20
_MANIFEST_NAME = "MANIFEST.json"


class DataError(Exception):
    """Base for every refusal in this module."""


class LicenceError(DataError):
    """The set may not be redistributed, or may not be fetched automatically.

    Always carries the operative clause; skipping a set means catching this explicitly.
    """


class IntegrityError(DataError):
    """Bytes on disk, or bytes just downloaded, do not match the pinned sha256."""


class UnknownSetError(DataError, KeyError):
    """No such `set_id` in the manifest."""

    def __str__(self) -> str:  # KeyError would quote the message
        return self.args[0] if self.args else ""


@dataclass(frozen=True)
class SetInfo:
    """One row of `available()`."""

    set_id: str
    status: str          # "ship" | "fetch" | "excluded"
    title: str
    licence: str
    rows: int | None
    bytes: int | None
    doi: str | None
    url: str | None

    @property
    def shipped(self) -> bool:
        return self.status == "ship"


# ---------------------------------------------------------------------------
# Locating the data
# ---------------------------------------------------------------------------

def benchmarks_dir() -> Path:
    """Where `benchmarks/` lives: `$EPIGAUGE_BENCHMARKS`, then the installed copy, then
    the sibling directory in a source checkout. Raises rather than returning a path that
    does not exist.
    """
    env = os.environ.get("EPIGAUGE_BENCHMARKS")
    here = Path(__file__).resolve().parent
    candidates = [Path(env)] if env else []
    candidates += [here / "benchmarks", here.parent / "benchmarks"]
    for c in candidates:
        if (c / _MANIFEST_NAME).is_file():
            return c
    tried = "\n  ".join(str(c) for c in candidates)
    raise DataError(
        f"benchmarks/{_MANIFEST_NAME} not found. Tried:\n  {tried}\n"
        "Set EPIGAUGE_BENCHMARKS to the directory containing MANIFEST.json, or "
        "reinstall epigauge with its package data."
    )


_MANIFEST_CACHE: dict[Path, dict] = {}


def manifest() -> dict:
    """The parsed `MANIFEST.json`. Cached; treat the result as read-only."""
    d = benchmarks_dir()
    if d not in _MANIFEST_CACHE:
        with (d / _MANIFEST_NAME).open(encoding="utf-8") as fh:
            _MANIFEST_CACHE[d] = json.load(fh)
    return _MANIFEST_CACHE[d]


def corpus_fingerprint(*, verify: bool = False) -> str:
    """The 16-hex fingerprint of the shipped corpus.

    `verify=True` recomputes it from the files on disk and raises `IntegrityError` on a
    mismatch: the one call that answers "is this the benchmark the paper reported?".
    """
    m = manifest()
    claimed = m["corpus_fingerprint"]
    if not verify:
        return claimed
    shipped = {k: v for k, v in m["sets"].items() if v["status"] == "ship"}
    parts = []
    for set_id, entry in sorted(shipped.items()):
        digest = _verified_sha256(benchmarks_dir() / entry["path"], entry["sha256"],
                                  set_id)
        line = f"{set_id}\t{digest}\t{entry['rows']}"
        if "scores" in entry:
            line += "\t" + _verified_sha256(
                benchmarks_dir() / entry["scores"]["path"],
                entry["scores"]["sha256"], f"{set_id} scores")
        parts.append(line)
    actual = hashlib.sha256("\n".join(parts).encode()).hexdigest()[:16]
    if actual != claimed:
        raise IntegrityError(
            f"corpus fingerprint mismatch: manifest claims {claimed}, the files on "
            f"disk give {actual}. The row set has changed under a stable name; this "
            "is not the benchmark the manifest describes."
        )
    return actual


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _verified_sha256(path: Path, expected: str, what: str) -> str:
    if not path.is_file():
        raise DataError(f"{what}: missing shipped file {path}")
    actual = _sha256(path)
    if actual != expected:
        raise IntegrityError(
            f"{what}: sha256 mismatch for {path}\n"
            f"  manifest expects {expected}\n"
            f"  file on disk is  {actual}\n"
            "Refusing to load. A file that does not match the manifest is a different "
            "benchmark, whatever it is called."
        )
    return actual


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _entry(set_id: str) -> dict:
    sets = manifest()["sets"]
    try:
        return sets[set_id]
    except KeyError:
        known = ", ".join(sorted(sets))
        raise UnknownSetError(f"unknown set_id {set_id!r}. Known sets: {known}") from None


def available(*, status: str | None = None) -> list[SetInfo]:
    """Every set the manifest knows, shipped or not.

    A set you may not have is listed with `status='fetch'` or `status='excluded'` and
    its licence, so the reason is discoverable before the evaluation is written.

    >>> sorted(s.set_id for s in available(status="ship"))[:2]
    ['hbv_hlabc_screen', 'independent_composite']
    """
    out = []
    for set_id, e in sorted(manifest()["sets"].items()):
        if status is not None and e["status"] != status:
            continue
        src = e.get("source") or {}
        out.append(SetInfo(
            set_id=set_id,
            status=e["status"],
            title=e["title"],
            licence=e["licence"],
            rows=e.get("rows"),
            bytes=e.get("bytes"),
            doi=src.get("doi"),
            url=src.get("url"),
        ))
    return out


def describe(set_id: str) -> dict:
    """The full manifest entry for one set: licence, basis, provenance, caveats."""
    return json.loads(json.dumps(_entry(set_id)))


def score_columns(set_id: str) -> dict:
    """Published competitor columns for a set: which shipped, which did not, and why.

    `withheld` gives each column's licence and the clause that blocks it; a baseline
    reported as "missing" is withheld.
    """
    e = _entry(set_id)
    cols = manifest()["predictor_columns"]
    sc = e.get("scores")
    if sc is None:
        return {"shipped": {}, "withheld": {}, "note":
                f"{set_id} ships no predictor score sidecar."}
    return {
        "shipped": {c: cols[c] for c in sc["columns_shipped"]},
        "withheld": {c: cols[c] for c in sc["columns_withheld"]},
        "note": sc["note"],
    }


def load(set_id: str, *, scores: bool = False) -> list[dict]:
    """Return the canonical rows of a shipped set, after verifying its sha256.

    Raises `LicenceError` for a set that is not shipped -- never an empty list -- and
    `IntegrityError` if the file does not match `MANIFEST.json`. With `scores=True`,
    each row gains a `scores` key holding the cleared competitor columns; withheld
    columns are absent, not null (`score_columns`).
    """
    e = _entry(set_id)
    if e["status"] != "ship":
        raise LicenceError(_refusal(set_id, e, verb="redistributed with epigauge"))

    root = benchmarks_dir()
    path = root / e["path"]
    _verified_sha256(path, e["sha256"], set_id)

    rows = []
    with path.open(encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if line:
                try:
                    rows.append(json.loads(line))
                except ValueError as exc:
                    raise IntegrityError(f"{set_id}: bad JSON at line {n}: {exc}") from None
    if len(rows) != e["rows"]:
        raise IntegrityError(
            f"{set_id}: manifest declares {e['rows']} rows, file holds {len(rows)}"
        )

    if scores:
        sc = e.get("scores")
        if sc is None:
            raise DataError(
                f"{set_id} ships no predictor score sidecar; call with scores=False. "
                "See score_columns() for what is withheld and why."
            )
        spath = root / sc["path"]
        _verified_sha256(spath, sc["sha256"], f"{set_id} scores")
        by_id: dict[str, dict] = {}
        with spath.open(encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    by_id[d["record_id"]] = d["scores"]
        missing = [r["record_id"] for r in rows if r["record_id"] not in by_id]
        if missing:
            raise IntegrityError(
                f"{set_id}: {len(missing)} rows have no score record "
                f"(first: {missing[0]}). The sidecar does not cover the row set."
            )
        for r in rows:
            r["scores"] = by_id[r["record_id"]]
    return rows


def fetch(set_id: str, *, dest_dir: str | Path | None = None,
          timeout: int = 120) -> Path:
    """Download a set we may not redistribute, from the upstream you may access.

    Verifies the downloaded bytes against the pinned `sha256` and deletes them on a
    mismatch. Raises `LicenceError`, naming the operative clause, when a set cannot be
    fetched automatically -- no licence permits it, access is agreement-gated, or the
    deposit is a DOI landing page rather than a direct file -- and for shipped sets,
    which `load` already verifies on disk. It never returns a path to nothing.
    """
    e = _entry(set_id)
    status = e["status"]
    if status == "ship":
        raise LicenceError(
            f"{set_id} ships with epigauge and is already sha256-verified on disk; "
            "use load(). Re-fetching it from upstream would create a second copy that "
            "can drift from the manifest."
        )
    if status == "excluded":
        raise LicenceError(_refusal(set_id, e, verb="fetched by epigauge"))

    src = e["source"]
    url = src.get("url") or ""
    expected = src.get("upstream_sha256")
    if not url.startswith(("http://", "https://")) or src.get("member"):
        raise LicenceError(
            f"{set_id} cannot be fetched automatically.\n"
            f"  licence: {e['licence']}\n"
            f"  reason:  {e['licence_basis']}\n"
            f"  source:  {url}\n"
            + (f"  file:    {src['member']}\n" if src.get("member") else "")
            + (f"  pinned sha256 of the file this benchmark was defined on:\n"
               f"           {expected}\n" if expected else "")
            + "The deposit is a DOI landing page, not a direct file. Download it "
              "yourself, confirm the sha256 above, then point epigauge at the result. "
              "epigauge will not guess a download URL for data whose provenance is "
              "not recorded."
        )
    if not expected:
        raise LicenceError(
            f"{set_id} has no pinned upstream sha256, so a fetch could not be "
            "verified and will not be attempted. Unverifiable bytes are not a "
            "benchmark."
        )

    dest_dir = Path(dest_dir) if dest_dir else _cache_dir() / set_id
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / url.rsplit("/", 1)[-1]

    if dest.is_file() and _sha256(dest) == expected:
        return dest

    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp, tmp.open("wb") as fh:
            while chunk := resp.read(_CHUNK):
                fh.write(chunk)
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise DataError(
            f"{set_id}: download failed from {url}: {exc}\n"
            f"  licence: {e['licence']}\n"
            "If the upstream has gone dark, the shipped CC BY 4.0 and CC0 tiers are "
            "self-contained and still evaluate on their own."
        ) from None

    actual = _sha256(tmp)
    if actual != expected:
        tmp.unlink(missing_ok=True)
        raise IntegrityError(
            f"{set_id}: upstream has CHANGED.\n"
            f"  url              {url}\n"
            f"  manifest expects {expected}\n"
            f"  upstream now is  {actual}\n"
            "The downloaded file has been deleted. This is not the row set this "
            "benchmark was defined on, and evaluating against it would report a "
            "different benchmark under the same name."
        )
    tmp.replace(dest)
    return dest


def _cache_dir() -> Path:
    base = os.environ.get("EPIGAUGE_CACHE")
    if base:
        return Path(base)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "epigauge"


def _refusal(set_id: str, e: dict, *, verb: str) -> str:
    src = e.get("source") or {}
    lines = [
        f"{set_id} ({e['title']}) may not be {verb}.",
        f"  licence: {e['licence']}",
        f"  reason:  {e['licence_basis']}",
    ]
    if src.get("url"):
        lines.append(f"  source:  {src['url']}")
    if e["status"] == "fetch":
        lines.append("  This set is fetchable: call fetch() instead of load().")
    else:
        lines.append(
            "  No fetch helper is offered for this set. Obtain written permission "
            "from the rights holder, or evaluate on the shipped sets -- "
            "available(status='ship') lists them."
        )
    return "\n".join(lines)
