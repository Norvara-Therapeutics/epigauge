#!/usr/bin/env python3
"""Validate and score one benchmark submission. The submission is a PULL REQUEST.

    python3 -I scripts/score_submission.py submissions/example/2026-11-13-hydropathy-demo.yaml
    python3 -I scripts/score_submission.py SUBMISSION.yaml --schema-only
    python3 -I scripts/score_submission.py SUBMISSION.yaml --predictions local.csv
    python3 -I scripts/score_submission.py SUBMISSION.yaml --comment out.md --ledger-row row.csv

THERE IS NO SERVER AND NO EMAIL ADDRESS. A submission is a pull request adding exactly
one file, `submissions/<team>/<YYYY-MM-DD>-<method>.yaml`. The predictions themselves are
NOT committed -- they are large -- so the YAML carries a public URL and the sha256 of the
bytes at it. This script downloads them, checks the digest, and scores them. Git history
is the evaluation ledger: every submission ever scored is a commit, so the count of times
this benchmark has been read is public and auditable with nothing hosted.

WHY A RESTRICTED YAML PARSER AND NOT PyYAML. The package has ZERO runtime dependencies
and keeps them (see pyproject.toml). The submission file is a FLAT subset of YAML --
whole-line comments, `key: value` scalars, and `key:` followed by `  - item` list entries
-- parsed by `parse_flat_yaml` below in about forty lines of standard library. Anything
outside that subset (nesting, flow collections `{}`/`[]`, anchors, aliases, tabs, block
scalars, inline comments, document markers) is REJECTED BY NAME rather than tolerated.
A tolerant parser on a submission file is the same bug class as a tolerant key join: it
turns a submitter's mistake into a clean-looking number on the wrong input.

WHAT IS CHECKED, AND WHY EACH ONE IS HERE.

  path matches fields   `submissions/<team>/<date>-<method>.yaml`, with `team` and
                        `method` equal to the fields inside. A file whose path disagrees
                        with its contents lands a row under the wrong name.
  no unknown keys       A misspelled key is otherwise a silently-ignored declaration.
  score direction       DECLARED, never inferred. No default exists. Reading a family of
                        percentile-rank columns with the score sign convention produced
                        an apparent 0.4650 "anti-predictive floor" that vanished when
                        flipped; inferring the direction means fitting it on the
                        evaluation labels.
  sha256 pin            The bytes scored are the bytes declared. An unpinned URL means
                        the number published and the file behind it can diverge later.
  host allowlist        Downloads only from GitHub release assets and Zenodo. CI follows
                        a submitter-supplied URL, redirects included, so the allowlist is
                        re-checked on every hop.
  one score per row_id  Delegated to `epigauge.cli.validate_submission`, which is the
                        same check `epigauge validate` runs. Not reimplemented here.

SCORING IS NOT REIMPLEMENTED. The scorecard comes from `epigauge.cli` ("score"), so a
submitter running `epigauge score` locally and CI running this script print the same
table from the same code: stratified estimator ranked, pooled beside it, baselines as
rows, and the minimum detectable difference at the top.

Exit codes follow the package: 0 clean, 1 a check FAILED, 2 could not run. 2 is never
conflated with 0 -- a download that did not happen is not a submission that passed.

Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from contextlib import redirect_stdout
from datetime import date, datetime, timezone
from pathlib import Path

# `scripts/` sits beside `epigauge/`, so the package resolves without installation.
# Required under `python3 -I`, which does not put the script's directory on sys.path --
# and that exclusion is wanted here: nothing adjacent to this file can shadow an import.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

EXIT_OK = 0
EXIT_CHECK_FAILED = 1
EXIT_COULD_NOT_RUN = 2

#: Status tokens, matching `epigauge.report`. Re-declared rather than imported so that
#: `--help` and the schema checks work even if the package import fails.
PASS, FAIL, UNRUNNABLE = "pass", "FAIL", "could-not-run"

SCHEMA_VERSION = 1

#: Predictions and training-hash files are fetched only from these hosts. A submitter
#: URL is attacker-controlled input to CI; the pin proves WHAT was scored, the allowlist
#: bounds WHERE CI will go to get it.
ASSET_HOSTS = frozenset({
    "github.com",
    "objects.githubusercontent.com",
    "release-assets.githubusercontent.com",
    "raw.githubusercontent.com",
    "zenodo.org",
    "sandbox.zenodo.org",
})

#: 64 MiB. A `row_id,score` CSV for the largest shipped set is well under a megabyte;
#: this ceiling exists so a wrong URL fails fast instead of filling the runner's disk.
MAX_BYTES = 64 * 1024 * 1024

DIRECTIONS = ("higher-is-better", "lower-is-better")
OVERLAP_DECLARATIONS = ("public", "private", "none")

_SLUG = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_HANDLE = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,37}[A-Za-z0-9])?$")
_ISO_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_EMAILISH = re.compile(r"[^\s@]+@[^\s@]+\.[^\s@]+")
_KEY_LINE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*):(?:[ ]+(.*))?$")
_ITEM_LINE = re.compile(r"^[ ]{2}-[ ]+(.*)$")

#: Scalar fields every submission must carry, and what each one is for.
REQUIRED_SCALARS = {
    "schema_version": "submission schema version; this script accepts 1",
    "team": "team slug; must equal the submissions/<team>/ directory name",
    "contact": "a GitHub handle, an https URL or an email address",
    "method": "method slug; must equal the <method> part of the filename",
    "version": "method version. Each version is a separate ledger entry",
    "benchmark": "row set name, as `epigauge list` prints it",
    "predictions_url": "public https URL of the predictions CSV (release asset or Zenodo)",
    "predictions_sha256": "sha256 of the bytes at predictions_url",
    "score_direction": f"one of {' / '.join(DIRECTIONS)}. No default exists",
    "uses_peptide_as_input": "true or false. false is a LEGAL submission",
    "training_cutoff": "YYYY-MM-DD; the cutoff date of the training data",
    "training_overlap": f"one of {' / '.join(OVERLAP_DECLARATIONS)}",
}

#: Required only when `training_overlap` is `public` or `private`.
CONDITIONAL_SCALARS = {
    "training_hashes_url": "public https URL of the sorted sha256 training-pair hashes",
    "training_hashes_sha256": "sha256 of the bytes at training_hashes_url",
}

OPTIONAL_SCALARS = {
    "notes": "one line of free text, shown on the scorecard",
}

OPTIONAL_LISTS = {
    "training_sources": "one corpus per entry, each with a version or access date",
}

ALL_KEYS = (set(REQUIRED_SCALARS) | set(CONDITIONAL_SCALARS)
            | set(OPTIONAL_SCALARS) | set(OPTIONAL_LISTS))


class SchemaError(ValueError):
    """The submission file is present but is not a submission. Exit 1."""


class CouldNotRun(RuntimeError):
    """The check could not be attempted. Always exit 2, never exit 0."""


# ---------------------------------------------------------------------------
# the restricted parser


def parse_flat_yaml(text: str, *, source: str = "<submission>") -> dict:
    """Parse the flat YAML subset. Everything outside it raises SchemaError by name.

    Accepted, and nothing else:

        # a whole-line comment
        key: scalar value
        key:
          - list item
          - list item

    Scalars may be wrapped in matching single or double quotes, which are stripped. A
    scalar with an empty value is an error, not an empty string: `key:` alone opens a
    list. There are no inline comments -- a `#` is only a comment at the start of a line,
    because a `#` inside a URL is part of the URL.
    """
    data: dict[str, object] = {}
    open_list: list[str] | None = None

    for lineno, raw in enumerate(text.splitlines(), 1):
        where = f"{source}:{lineno}"
        if "\t" in raw:
            raise SchemaError(f"{where}: tab character. Indent with spaces only.")
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        if raw.strip() in ("---", "..."):
            raise SchemaError(
                f"{where}: document marker {raw.strip()!r}. A submission file is one "
                "flat document; multi-document streams are not accepted."
            )

        item = _ITEM_LINE.match(raw)
        if item is not None:
            if open_list is None:
                raise SchemaError(
                    f"{where}: list item with no list key above it. A `  - item` line "
                    "must follow a `key:` line with an empty value."
                )
            open_list.append(_scalar(item.group(1), where))
            continue

        if raw[0] in " ":
            raise SchemaError(
                f"{where}: indented line {raw.strip()!r}. The submission schema is flat: "
                "nested mappings are not accepted, only `key: value` and `  - item`."
            )

        key_match = _KEY_LINE.match(raw.rstrip())
        if key_match is None:
            raise SchemaError(
                f"{where}: not a `key: value` line: {raw.strip()!r}. Keys are "
                "lowercase_with_underscores and a value is separated by `: `."
            )
        key, value = key_match.group(1), key_match.group(2)
        if key in data:
            raise SchemaError(
                f"{where}: duplicate key {key!r}. The later value would silently win."
            )

        if value is None or not value.strip():
            open_list = []
            data[key] = open_list
        else:
            open_list = None
            data[key] = _scalar(value, where)

    return data


def _scalar(raw: str, where: str) -> str:
    value = raw.strip()
    if not value:
        raise SchemaError(f"{where}: empty value. Remove the key or give it a value.")
    if value[0] in "\"'" or value[-1] in "\"'":
        if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
            return value[1:-1]
        raise SchemaError(f"{where}: unbalanced quote in {raw.strip()!r}.")
    if value[0] in "[{":
        raise SchemaError(
            f"{where}: flow collection {value!r}. Use `key:` plus `  - item` lines."
        )
    if value[0] in "&*!":
        raise SchemaError(
            f"{where}: YAML anchor/alias/tag {value!r} is not accepted. Write the value out."
        )
    if value[:1] in ("|", ">") and value[1:2] in ("", "-", "+"):
        raise SchemaError(
            f"{where}: block scalar {value!r} is not accepted. Values are single-line."
        )
    return value


# ---------------------------------------------------------------------------
# schema checks
#
# Every check is appended to a list and the WHOLE list is reported, so a submitter sees
# every problem in one CI run instead of one per push.


def _bool(value: str, key: str) -> bool:
    if value in ("true", "false"):
        return value == "true"
    raise SchemaError(
        f"{key}: {value!r} is not a boolean. Write exactly `true` or `false` -- "
        "`yes`/`on`/`1` are not accepted, because a tolerant reading of a declaration "
        "is how a declaration stops meaning anything."
    )


def _check_url(value: str, key: str) -> None:
    parts = urllib.parse.urlsplit(value)
    if parts.scheme != "https":
        raise SchemaError(f"{key}: must be an https:// URL, got {value!r}.")
    host = (parts.hostname or "").lower()
    if host not in ASSET_HOSTS:
        raise SchemaError(
            f"{key}: host {host!r} is not an accepted asset host. Publish the file as a "
            f"GitHub release asset or a Zenodo record. Accepted: "
            f"{', '.join(sorted(ASSET_HOSTS))}."
        )


def check_schema(data: dict, path: Path) -> list[tuple[str, str, str]]:
    """Return `[(check, status, detail), ...]`, in the shape `report.render_validation` takes."""
    checks: list[tuple[str, str, str]] = []

    def ok(name, detail):
        checks.append((name, PASS, detail))

    def bad(name, detail):
        checks.append((name, FAIL, detail))

    # --- keys present, and nothing unknown ----------------------------------
    missing = [k for k in REQUIRED_SCALARS if k not in data]
    unknown = [k for k in data if k not in ALL_KEYS]
    if missing or unknown:
        bits = []
        if missing:
            bits.append("missing: " + ", ".join(missing))
        if unknown:
            bits.append("unknown: " + ", ".join(sorted(unknown)))
        bad("every required key present, no unknown keys",
            "; ".join(bits) + ". An unknown key is a misspelled declaration that would "
            "otherwise be silently ignored.")
    else:
        ok("every required key present, no unknown keys",
           f"{len(REQUIRED_SCALARS)} required keys")

    def field(key):
        value = data.get(key)
        return value if isinstance(value, str) else None

    # --- schema version -----------------------------------------------------
    raw_version = field("schema_version")
    if raw_version != str(SCHEMA_VERSION):
        bad("schema_version is supported",
            f"{raw_version!r}; this script accepts schema_version: {SCHEMA_VERSION}.")
    else:
        ok("schema_version is supported", raw_version)

    # --- the path and the contents agree ------------------------------------
    team, method = field("team"), field("method")
    rel = path.as_posix()
    parts = path.parts
    detail = []
    if len(parts) < 3 or parts[-3] != "submissions":
        detail.append(f"{rel} is not submissions/<team>/<YYYY-MM-DD>-<method>.yaml")
    if path.suffix != ".yaml":
        detail.append(f"extension {path.suffix!r} is not .yaml")
    stem = path.stem
    file_date, file_method = (stem[:10], stem[11:]) if len(stem) > 11 else ("", "")
    if not _ISO_DATE.match(file_date) or stem[10:11] != "-":
        detail.append(f"filename {path.name!r} does not start with YYYY-MM-DD-")
    else:
        try:
            date.fromisoformat(file_date)
        except ValueError:
            detail.append(f"{file_date} is not a real date")
    if team is not None and len(parts) >= 2 and parts[-2] != team:
        detail.append(f"directory {parts[-2]!r} != team {team!r}")
    if method is not None and file_method and file_method != method:
        detail.append(f"filename method {file_method!r} != method {method!r}")
    if detail:
        bad("the path matches the declared team, method and date",
            "; ".join(detail) + ". A file whose path disagrees with its contents lands "
            "a ledger row under the wrong name.")
    else:
        ok("the path matches the declared team, method and date", rel)

    # --- slugs --------------------------------------------------------------
    for key in ("team", "method"):
        value = field(key)
        if value is None:
            continue
        if _SLUG.match(value):
            ok(f"{key} is a slug", value)
        else:
            bad(f"{key} is a slug",
                f"{value!r} is not lowercase-with-hyphens. It is used as a path and as a "
                "leaderboard label, so it has to be stable and unambiguous.")

    version = field("version")
    if version is not None:
        if version.strip() and len(version) <= 64:
            ok("version is declared", version)
        else:
            bad("version is declared", f"{version!r} is empty or over 64 characters.")

    # --- contact: a handle, an https URL or an email ------------------------
    contact = field("contact")
    if contact is not None:
        if (_HANDLE.match(contact) or contact.startswith("https://")
                or _EMAILISH.search(contact)):
            ok("contact is a GitHub handle, an https URL or an email", contact)
        else:
            bad("contact is a GitHub handle, an https URL or an email",
                f"{contact!r} is none of those.")

    # --- the rule the whole benchmark turns on ------------------------------
    direction = field("score_direction")
    if direction in DIRECTIONS:
        ok("a score direction is declared", direction)
    else:
        bad("a score direction is declared",
            f"{direction!r} is not one of {' / '.join(DIRECTIONS)}. There is no default "
            "and it is not inferred: picking the orientation that scores better fits a "
            "free parameter on the evaluation labels, and it makes a sign error "
            "indistinguishable from a result. A family of percentile-rank columns read "
            "with the score convention produced a 0.4650 'anti-predictive floor' that "
            "vanished when flipped.")

    # --- uses_peptide_as_input ----------------------------------------------
    upi = field("uses_peptide_as_input")
    if upi is not None:
        try:
            value = _bool(upi, "uses_peptide_as_input")
            ok("uses_peptide_as_input is declared",
               f"{str(value).lower()}"
               + ("" if value else "  (legal: an allele lookup with no peptide input "
                                   "reaches pooled 0.7122 at within-allele exactly 0.5000, "
                                   "and it ships as a ranked row)"))
        except SchemaError as error:
            bad("uses_peptide_as_input is declared", str(error))

    # --- dates --------------------------------------------------------------
    cutoff = field("training_cutoff")
    if cutoff is not None:
        try:
            parsed = date.fromisoformat(cutoff)
            if not _ISO_DATE.match(cutoff):
                raise ValueError(cutoff)
            if parsed > date.today():
                bad("training_cutoff is a past ISO date",
                    f"{cutoff} is in the future.")
            else:
                ok("training_cutoff is a past ISO date", cutoff)
        except ValueError:
            bad("training_cutoff is a past ISO date",
                f"{cutoff!r} is not YYYY-MM-DD.")

    # --- URLs and digests ---------------------------------------------------
    overlap = field("training_overlap")
    needs_hashes = overlap in ("public", "private")
    url_fields = [("predictions_url", "predictions_sha256")]
    if needs_hashes:
        url_fields.append(("training_hashes_url", "training_hashes_sha256"))
    elif any(k in data for k in CONDITIONAL_SCALARS):
        bad("training-hash fields match the declaration",
            f"training_overlap is {overlap!r} but training_hashes_* is present. Declare "
            "`public` or `private` to have the overlap measured.")

    if overlap in OVERLAP_DECLARATIONS:
        ok("a training-overlap declaration is present",
           overlap + ("  (reported as `overlap: undeclared`)" if overlap == "none" else ""))
    else:
        bad("a training-overlap declaration is present",
            f"{overlap!r} is not one of {' / '.join(OVERLAP_DECLARATIONS)}.")

    for url_key, sha_key in url_fields:
        url, sha = field(url_key), field(sha_key)
        if url is None:
            bad(f"{url_key} is a pinned asset URL",
                f"{url_key} is required when training_overlap is {overlap!r}.")
        else:
            try:
                _check_url(url, url_key)
                ok(f"{url_key} is a pinned asset URL", url)
            except SchemaError as error:
                bad(f"{url_key} is a pinned asset URL", str(error))
        if sha is None or not _SHA256.match(sha):
            bad(f"{sha_key} is a sha256 digest",
                f"{sha!r} is not 64 lowercase hex characters. An unpinned URL means the "
                "number published and the file behind it can diverge afterwards.")
        else:
            ok(f"{sha_key} is a sha256 digest", sha[:16] + "...")

    # --- lists --------------------------------------------------------------
    for key in OPTIONAL_LISTS:
        value = data.get(key)
        if value is None:
            continue
        if not isinstance(value, list) or not value:
            bad(f"{key} is a non-empty list",
                f"{key} is present but empty. Remove it or list the corpora.")
        else:
            ok(f"{key} is a non-empty list", f"{len(value)} entries")

    return checks


# ---------------------------------------------------------------------------
# fetching, pinned


class _AllowlistRedirects(urllib.request.HTTPRedirectHandler):
    """Re-check the host on every redirect hop, not just the first."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        _check_url(newurl, "redirect target")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def fetch_pinned(url: str, expected_sha256: str, destination: Path, *,
                 max_bytes: int = MAX_BYTES) -> int:
    """Download `url` to `destination`, verifying the digest. Returns the byte count.

    A digest mismatch raises SchemaError (exit 1: the submission is wrong). A network or
    HTTP failure raises CouldNotRun (exit 2: nothing was measured). The two are never
    collapsed -- "could not download" is not "failed the check", and neither is a pass.
    """
    _check_url(url, "predictions_url")
    opener = urllib.request.build_opener(_AllowlistRedirects)
    digest = hashlib.sha256()
    total = 0
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with opener.open(url, timeout=120) as response, destination.open("wb") as out:
            while True:
                chunk = response.read(1 << 16)
                if not chunk:
                    break
                total += len(chunk)
                if total > max_bytes:
                    raise CouldNotRun(
                        f"{url} exceeds the {max_bytes // (1024 * 1024)} MiB ceiling. A "
                        "predictions file is two columns; check the URL."
                    )
                digest.update(chunk)
                out.write(chunk)
    except urllib.error.HTTPError as error:
        raise CouldNotRun(f"{url}: HTTP {error.code} {error.reason}. Nothing was "
                          "downloaded, so nothing passed.") from None
    except (urllib.error.URLError, OSError, ValueError) as error:
        raise CouldNotRun(f"{url}: {type(error).__name__}: {error}") from None

    actual = digest.hexdigest()
    if actual != expected_sha256:
        raise SchemaError(
            f"sha256 mismatch for {url}\n"
            f"  declared: {expected_sha256}\n"
            f"  actual:   {actual}  ({total} bytes)\n"
            "The declared digest is what makes the published number reproducible. Either "
            "the file changed after the digest was taken, or the wrong file is linked."
        )
    return total


# ---------------------------------------------------------------------------
# training overlap


def read_hash_file(path: Path) -> tuple[set[str], list[str]]:
    """Sorted, deduplicated lowercase sha256 digests, one per line. Returns (set, problems)."""
    digests: set[str] = set()
    problems: list[str] = []
    previous = ""
    unsorted = False
    for lineno, line in enumerate(path.read_text(encoding="utf-8",
                                                 errors="replace").splitlines(), 1):
        token = line.strip()
        if not token:
            continue
        if not _SHA256.match(token):
            if len(problems) < 5:
                problems.append(f"line {lineno}: {token[:32]!r} is not a sha256 digest")
            continue
        if token < previous:
            unsorted = True
        previous = token
        digests.add(token)
    if unsorted:
        problems.append("the digests are not in sorted order")
    return digests, problems


def benchmark_pair_hashes(set_id: str) -> tuple[set[str], int]:
    """sha256(f'{peptide}|{allele}') over CANONICALISED benchmark keys.

    Canonicalisation is not optional here. A naive join between two spellings of the same
    reports a clean-looking zero and means nothing.
    """
    from epigauge import canon, data

    digests: set[str] = set()
    skipped = 0
    for row in data.load(set_id):
        try:
            peptide = canon.canonicalise_peptide(str(row["peptide"]))
            allele = canon.canonicalise_allele(str(row["allele"]))
        except (canon.CanonError, KeyError):
            skipped += 1
            continue
        digests.add(hashlib.sha256(f"{peptide}|{allele}".encode()).hexdigest())
    return digests, skipped


def overlap_report(declared: set[str], benchmark: set[str], skipped: int) -> str:
    """One line. A zero always arrives with the counts that prove something was compared."""
    if not benchmark:
        return ("overlap: COULD NOT RUN -- zero benchmark keys canonicalised, so nothing "
                "was compared. This is not 0%.")
    shared = declared & benchmark
    pct = 100.0 * len(shared) / len(benchmark)
    # A non-zero intersection must never print as "0.0%". A printed zero on this
    # benchmark has to come from a counter that really is zero, so the counts travel
    # with the percentage and a small non-zero share says so explicitly.
    shown = "<0.01%" if shared and pct < 0.01 else f"{pct:.2f}%"
    tail = (f"  ({skipped} benchmark row(s) not canonicalisable, excluded)"
            if skipped else "")
    return (f"overlap: {shown}  -- {len(shared)} of {len(benchmark)} benchmark "
            f"(peptide, HLA) keys appear among the {len(declared)} declared training "
            f"hashes{tail}")


# ---------------------------------------------------------------------------
# scoring, delegated to the package


def score_predictions(predictions: Path, benchmark: str, direction: str) -> tuple[int, str]:
    """Run `epigauge score` in-process and capture the scorecard. Returns (exit_code, text)."""
    try:
        from epigauge import cli
    except ImportError as error:  # pragma: no cover
        raise CouldNotRun(
            f"could not import epigauge from {ROOT}: {error}. Exiting "
            f"{EXIT_COULD_NOT_RUN} rather than reporting a figure computed without it."
        ) from None

    buffer = io.StringIO()
    with redirect_stdout(buffer):
        code = cli.main(["score", str(predictions), "--set", benchmark,
                         "--direction", direction])
    return code, buffer.getvalue()


def ledger_figures(predictions: Path, benchmark: str, direction: str) -> dict:
    """The submission's own two figures, from the same estimator the scorecard ranks on."""
    from epigauge import cli

    rows = cli._load_rows(cli._resolve_set(benchmark))
    sub = cli.read_submission(predictions, direction=direction)
    name, scores = cli._oriented(sub.label, sub.direction, cli._join(sub, rows))
    row = cli._method_row(name, cli._labels(rows), scores, cli._groups(rows),
                          sign_flipped=sub.direction == cli.LOWER,
                          direction_declared=True)
    if row.stratified is None:
        return {}
    return {
        "within_allele": f"{row.stratified.value:.4f}",
        "pooled": f"{row.pooled:.4f}",
        "strata": f"{row.stratified.n_strata_contributing}/{row.stratified.n_strata_total}",
        "evaluable_rows": f"{row.stratified.n_rows_evaluable}/{row.stratified.n_rows_total}",
    }


# ---------------------------------------------------------------------------
# rendering


def render_checks(path: str, checks: list[tuple[str, str, str]]) -> str:
    """`report.render_validation` when the package imports; a faithful fallback if not."""
    try:
        from epigauge import report
        return report.render_validation(path, checks)
    except Exception:  # pragma: no cover - only when the package is unavailable
        lines = [f"VALIDATION  {path}", ""]
        for name, status, detail in checks:
            lines.append(f"  [{status:>13}]  {name}")
            if detail:
                lines.append(f"                   {detail}")
        return "\n".join(lines) + "\n"


LEDGER_COLUMNS = (
    "scored_utc", "team", "method", "version", "benchmark", "score_direction",
    "uses_peptide_as_input", "within_allele", "pooled", "strata", "evaluable_rows",
    "training_cutoff", "training_overlap", "overlap_pct", "predictions_sha256",
)


def ledger_row(data: dict, figures: dict, overlap_pct: str) -> str:
    import csv

    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\n").writerow([
        datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        data.get("team", ""), data.get("method", ""), data.get("version", ""),
        data.get("benchmark", ""), data.get("score_direction", ""),
        data.get("uses_peptide_as_input", ""),
        figures.get("within_allele", ""), figures.get("pooled", ""),
        figures.get("strata", ""), figures.get("evaluable_rows", ""),
        data.get("training_cutoff", ""), data.get("training_overlap", ""),
        overlap_pct, data.get("predictions_sha256", ""),
    ])
    return buffer.getvalue()


def render_comment(data: dict, path: Path, body: str, verdict: str) -> str:
    """The PR comment. The scorecard goes in verbatim -- it is the artefact, not a summary."""
    return (
        f"## EpiGauge scorecard — `{path.as_posix()}`\n\n"
        f"**{verdict}**\n\n"
        f"| | |\n|---|---|\n"
        f"| team | `{data.get('team', '?')}` |\n"
        f"| method | `{data.get('method', '?')}` @ `{data.get('version', '?')}` |\n"
        f"| benchmark | `{data.get('benchmark', '?')}` |\n"
        f"| score direction | `{data.get('score_direction', '?')}` (declared) |\n"
        f"| uses peptide as input | `{data.get('uses_peptide_as_input', '?')}` |\n"
        f"| training cutoff | `{data.get('training_cutoff', '?')}` |\n"
        f"| predictions sha256 | `{data.get('predictions_sha256', '?')}` |\n\n"
        f"```\n{body.rstrip()}\n```\n\n"
        "<sub>Pooled is printed beside the stratified estimator and is never ranked on. "
        "Gaps below the minimum detectable difference print `=`, not a rank. "
        "Merging this PR appends one row to `submissions/ledger.csv`; the commit is the ledger "
        "entry.</sub>\n"
    )


# ---------------------------------------------------------------------------
# entry point


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="score_submission.py",
        description="Validate and score one submission YAML. A submission is a pull "
                    "request adding submissions/<team>/<YYYY-MM-DD>-<method>.yaml; the "
                    "predictions are downloaded from the URL it declares and checked "
                    "against the sha256 it declares. There is no server to maintain.",
        epilog="exit codes: 0 clean, 1 a check failed, 2 could not run. 2 is never "
               "conflated with 0.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("submission", metavar="SUBMISSION.yaml", nargs="?",
                        help="path to the submission file, relative to the repo root")
    parser.add_argument("--ledger-header", action="store_true",
                        help="print the submissions/ledger.csv header row and exit. Used to "
                             "create the ledger on first merge, so the columns are "
                             "defined in exactly one place")
    parser.add_argument("--schema-only", action="store_true",
                        help="validate the YAML and stop: no download, no scoring")
    parser.add_argument("--predictions", metavar="FILE", default=None,
                        help="score this local CSV instead of downloading. The declared "
                             "sha256 is still checked against it")
    parser.add_argument("--training-hashes", metavar="FILE", default=None,
                        help="use this local hash file instead of downloading it")
    parser.add_argument("--comment", metavar="FILE", default=None,
                        help="write the scorecard as PR-comment markdown")
    parser.add_argument("--ledger-row", metavar="FILE", default=None,
                        help="write the one submissions/ledger.csv row for this submission")
    parser.add_argument("--work-dir", metavar="DIR", default=None,
                        help="directory for downloaded files (default: a temp dir)")
    parser.add_argument("--max-bytes", type=int, default=MAX_BYTES, metavar="N",
                        help=f"download ceiling in bytes (default {MAX_BYTES})")
    return parser


def run(args) -> int:
    if args.ledger_header:
        print(",".join(LEDGER_COLUMNS))
        return EXIT_OK
    if not args.submission:
        raise CouldNotRun("no submission file given. Pass SUBMISSION.yaml, or "
                          "--ledger-header to print the ledger columns.")
    path = Path(args.submission)
    if not path.is_file():
        raise CouldNotRun(f"no such submission file: {path}. Nothing was read, so "
                          "nothing passed.")
    try:
        rel = path.resolve().relative_to(ROOT)
    except ValueError:
        rel = path

    data = parse_flat_yaml(path.read_text(encoding="utf-8"), source=path.as_posix())
    checks = check_schema(data, rel)
    sys.stdout.write(render_checks(rel.as_posix(), checks))

    if any(status == FAIL for _, status, _ in checks):
        print("refusing to score a submission that failed the schema check.")
        return EXIT_CHECK_FAILED
    if args.schema_only:
        print("schema check only; nothing was downloaded and no AUROC was computed.")
        return EXIT_OK

    import tempfile
    temp = None
    if args.work_dir:
        work = Path(args.work_dir)
        work.mkdir(parents=True, exist_ok=True)
    else:
        temp = tempfile.TemporaryDirectory(prefix="epigauge-submission-")
        work = Path(temp.name)

    try:
        # --- predictions, pinned -------------------------------------------
        label = f"{data['method']}-{data['version']}"
        predictions = work / f"{label}.csv"
        if args.predictions:
            source = Path(args.predictions)
            if not source.is_file():
                raise CouldNotRun(f"no such predictions file: {source}")
            predictions.write_bytes(source.read_bytes())
            actual = hashlib.sha256(predictions.read_bytes()).hexdigest()
            if actual != data["predictions_sha256"]:
                raise SchemaError(
                    f"sha256 mismatch for {source}\n"
                    f"  declared: {data['predictions_sha256']}\n"
                    f"  actual:   {actual}"
                )
            print(f"\npredictions: {source} (local), sha256 verified")
        else:
            size = fetch_pinned(data["predictions_url"], data["predictions_sha256"],
                                predictions, max_bytes=args.max_bytes)
            print(f"\npredictions: {data['predictions_url']}\n"
                  f"             {size} bytes, sha256 verified")

        # --- training overlap ----------------------------------------------
        overlap_line = ("overlap: undeclared  (training_overlap: none -- accepted, and "
                        "shown as `undeclared` on the leaderboard)")
        overlap_pct = "undeclared"
        if data.get("training_overlap") in ("public", "private"):
            hashes = work / "training-hashes.txt"
            if args.training_hashes:
                hashes.write_bytes(Path(args.training_hashes).read_bytes())
            else:
                fetch_pinned(data["training_hashes_url"], data["training_hashes_sha256"],
                             hashes, max_bytes=args.max_bytes)
            declared, problems = read_hash_file(hashes)
            if not declared:
                raise SchemaError(
                    "the declared training-hash file contains zero sha256 digests. A "
                    "declaration that compares nothing reports 0% overlap and means "
                    "nothing; declare `training_overlap: none` instead."
                )
            try:
                bench, skipped = benchmark_pair_hashes(data["benchmark"])
                overlap_line = overlap_report(declared, bench, skipped)
                if bench:
                    overlap_pct = f"{100.0 * len(declared & bench) / len(bench):.2f}"
            except Exception as error:
                overlap_line = (f"overlap: COULD NOT RUN -- {type(error).__name__}: "
                                f"{error}. This is not 0%.")
            for problem in problems:
                print(f"  note: training hash file: {problem}")

        # --- the scorecard, from epigauge.cli ------------------------------
        code, body = score_predictions(predictions, data["benchmark"],
                                       data["score_direction"])
        # The scorecard titles itself with the file it read, which here is a scratch
        # path in a temp directory. Show the submission instead: a PR comment carrying a
        # runner's temp path tells the reader nothing about what was scored.
        body = body.replace(str(predictions), f"{data['team']}/{label}")
        body = body.rstrip() + "\n\n" + overlap_line
        if data.get("notes"):
            body += f"\nsubmitter note: {data['notes']}"
        sys.stdout.write("\n" + body + "\n")

        verdict = {
            EXIT_OK: "Scored. Merging appends this row to `submissions/ledger.csv`.",
            EXIT_CHECK_FAILED: "A check FAILED. Not scored.",
        }.get(code, "COULD NOT RUN. This is not a pass.")

        if args.comment:
            Path(args.comment).write_text(render_comment(data, rel, body, verdict),
                                          encoding="utf-8")
            print(f"\nwrote PR comment markdown to {args.comment}")

        if args.ledger_row and code == EXIT_OK:
            figures = ledger_figures(predictions, data["benchmark"],
                                     data["score_direction"])
            Path(args.ledger_row).write_text(ledger_row(data, figures, overlap_pct),
                                             encoding="utf-8")
            print(f"wrote the submissions/ledger.csv row to {args.ledger_row}")

        return code
    finally:
        if temp is not None:
            temp.cleanup()


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except SchemaError as error:
        print(f"score_submission: CHECK FAILED: {error}", file=sys.stderr)
        return EXIT_CHECK_FAILED
    except CouldNotRun as error:
        print(f"score_submission: COULD NOT RUN: {error}", file=sys.stderr)
        return EXIT_COULD_NOT_RUN


if __name__ == "__main__":
    raise SystemExit(main())
