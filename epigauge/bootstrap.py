"""Clustered resampling, intervals, and the power gate. Standard library only.

Labels cluster by source protein: within-protein concordance is 0.815 against 0.635
expected under independence, over 507 pairs in 398 proteins, so the row is not the
independent unit. The design effect is MEASURED per comparison by `paired`, never scaled
by a stored factor; across the 833 rows it runs 1.00-1.08x (BigMHC_IM 1.00, MHCflurry
1.01, NetMHCpan-4.1 1.03, PRIME-2.0 1.04, TITANiAN-IM 1.05) and sub-1.0 is reported as is.
"""

from __future__ import annotations

import hashlib
import inspect
import math
import random
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import NormalDist

__all__ = [
    "MixedBasisError",
    "PreRegistrationUnderpowered",
    "UnmeasuredBasis",
    "Unmeasurable",
    "UNMEASURABLE_STUDY_CONSTANT",
    "SeedRecipe",
    "ClusterMap",
    "DrawMatrix",
    "Estimate",
    "Delta",
    "MDEReport",
    "assign_clusters",
    "study_level_basis",
    "draws",
    "interval",
    "paired",
    "mde",
    "format_table",
    "z_sum",
]


# ---------------------------------------------------------------------------
# Refusals. Each one is a raise, never a warning and never a default value.
# ---------------------------------------------------------------------------


class MixedBasisError(ValueError):
    """Two cells of one table do not share a resampling basis. Printed tables are refused."""


class PreRegistrationUnderpowered(ValueError):
    """A declared success threshold sits below the minimum detectable effect."""


class UnmeasuredBasis(ValueError):
    """A quantity was asked for on a basis that was never resampled."""


# ---------------------------------------------------------------------------
# The "unmeasurable" sentinel
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Unmeasurable:
    """A measurement the data cannot support, carrying WHY.

    Returned in a required field rather than omitted, so the field reads `unmeasurable`
    instead of "not relevant". Canonical instance: `study_id` is constant across all 2,093
    infectious rows, so a study bootstrap there has exactly one cluster.
    """

    reason: str

    def __str__(self) -> str:
        return "unmeasurable"

    def __bool__(self) -> bool:
        return False

    def explain(self) -> str:
        return f"unmeasurable: {self.reason}"


UNMEASURABLE_STUDY_CONSTANT = Unmeasurable(
    "study_id is constant across every row; a study-level bootstrap would have exactly "
    "one cluster, which is not a wide interval but no interval at all"
)


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SeedRecipe:
    """Order-independent seeds. `seed_b = first 64 bits of sha256(version|set|cmp|b)`.

    `comparison_id` enters the per-replicate seed so comparisons draw independently of
    each other and of evaluation order, but NOT `digest`: requiring cells to share a
    comparison id would make every multi-arm table a mixed-basis error.
    """

    version: str
    set_id: str
    comparison_id: str

    def seed(self, b: int) -> int:
        h = hashlib.sha256(
            f"{self.version}|{self.set_id}|{self.comparison_id}|{b}".encode()
        ).digest()
        return int.from_bytes(h[:8], "big")

    def digest(self, n: int = 8) -> str:
        """Short digest over `version|set_id` -- the comparability-defining prefix."""
        return hashlib.sha256(f"{self.version}|{self.set_id}".encode()).hexdigest()[:n]

    def derive(self, suffix: str) -> SeedRecipe:
        """A sibling recipe for a companion resample (the row-level arm of a design effect)."""
        return SeedRecipe(self.version, self.set_id, f"{self.comparison_id}|{suffix}")


# ---------------------------------------------------------------------------
# Cluster assignment
# ---------------------------------------------------------------------------


@dataclass
class ClusterMap:
    """Row -> cluster, with the assignment rate that makes it readable.

    `assigned_fraction` is carried, not derivable later: the same 2,093 rows give 4.6%
    against the human proteome and 85.9% against all nine, and the two maps are otherwise
    indistinguishable from their cluster lists.
    """

    cluster_of_row: list[int]
    members: list[list[int]]
    keys: list[str]
    n_rows: int
    n_assigned: int
    n_proteins: int
    n_singleton_clusters: int
    unit: str = "cluster"
    notes: dict = field(default_factory=dict)

    @property
    def n_clusters(self) -> int:
        return len(self.members)

    @property
    def assigned_fraction(self) -> float:
        return self.n_assigned / self.n_rows if self.n_rows else 0.0

    @property
    def provisional(self) -> bool:
        """Below 80% assigned the clustered label is aspirational; say so in the table."""
        return self.assigned_fraction < 0.80

    def __str__(self) -> str:
        flag = "  [PROVISIONAL: assigned fraction below 0.80]" if self.provisional else ""
        return (
            f"{self.n_clusters} clusters over {self.n_rows} rows  "
            f"assigned {self.n_assigned}/{self.n_rows} = {self.assigned_fraction:.3f}  "
            f"proteins {self.n_proteins}  singletons {self.n_singleton_clusters}{flag}"
        )

    @classmethod
    def from_keys(cls, keys, *, n_assigned: int | None = None) -> ClusterMap:
        """Build from an explicit per-row cluster key. Used by tests and by saved maps."""
        keys = list(keys)
        order: dict[object, int] = {}
        cluster_of_row: list[int] = []
        members: list[list[int]] = []
        for i, k in enumerate(keys):
            if k is None or k == "":
                raise ValueError(f"row {i} has an empty cluster key; use a SINGLETON key")
            if k not in order:
                order[k] = len(members)
                members.append([])
            j = order[k]
            cluster_of_row.append(j)
            members[j].append(i)
        singletons = sum(1 for m in members if len(m) == 1)
        n_assigned = len(keys) if n_assigned is None else n_assigned
        return cls(
            cluster_of_row=cluster_of_row,
            members=members,
            keys=[str(k) for k in order],
            n_rows=len(keys),
            n_assigned=n_assigned,
            n_proteins=len(members) - singletons,
            n_singleton_clusters=singletons,
        )

    @classmethod
    def rows_as_clusters(cls, n_rows: int) -> ClusterMap:
        """The row-level design as the degenerate clustering: the design-effect denominator,
        and the only legitimate one-row-per-cluster map."""
        m = cls(
            cluster_of_row=list(range(n_rows)),
            members=[[i] for i in range(n_rows)],
            keys=[f"ROW:{i}" for i in range(n_rows)],
            n_rows=n_rows,
            n_assigned=n_rows,
            n_proteins=0,
            n_singleton_clusters=n_rows,
            unit="row",
        )
        return m


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, x: str) -> None:
        self.parent.setdefault(x, x)

    def find(self, x: str) -> str:
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            # Deterministic merge target: the lexicographically smaller accession, so the
            # cluster key does not depend on FASTA file order.
            lo, hi = (ra, rb) if ra < rb else (rb, ra)
            self.parent[hi] = lo


def _parse_fasta(path):
    """Yield (accession, sequence). `>sp|P01234|NAME_HUMAN desc` -> `P01234`."""
    acc = None
    chunks: list[str] = []
    with open(path, "r", encoding="utf-8", errors="strict") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(chunks)
                header = line[1:].strip()
                parts = header.split("|")
                acc = parts[1] if len(parts) >= 3 else (header.split() or ["?"])[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if acc is not None:
        yield acc, "".join(chunks)


def assign_clusters(
    peptides,
    proteome_paths,
    *,
    k_window: int = 9,
    short_window: int = 8,
) -> ClusterMap:
    """Source-protein clusters by exact peptide-substring match against FASTA files.

    A peptide of length >= `k_window` is assigned to every protein containing any of its
    `k_window`-mer windows; one of exactly `short_window` residues is matched whole, since
    the 99 8-mers in the infectious set are part of the gap between 4.6% assignment (human
    proteome alone) and 85.9% (all nine). Shorter raises: a silent skip is indistinguishable
    from an unmatched peptide in `assigned_fraction`, the number this function exists to
    report.

    Multi-protein hits are merged by union-find, never resolved to one winner -- an
    arbitrary winner lets a replicate draw both proteins and count the shared row twice.
    Merging also makes the result independent of FASTA file order.

    Unassigned rows become singletons keyed by PEPTIDE, so duplicates of one unmatched
    peptide stay in one unit. Pooling them instead would make a single 295-row cluster, the
    largest in the design at 14.1% of the weight, out of which proteomes we happen to ship;
    dropping them would make `assigned_fraction` unobservable.
    """
    peptides = [p.upper() for p in peptides]
    paths = [Path(p) for p in proteome_paths]
    if not paths:
        raise UnmeasuredBasis(
            "no proteome supplied; assignment would be 0.000 and every row a singleton, "
            "which is the row-level design under a clustered name"
        )
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"proteome not found: {p}")

    want9: dict[str, set[int]] = defaultdict(set)
    want8: dict[str, set[int]] = defaultdict(set)
    for i, pep in enumerate(peptides):
        if len(pep) >= k_window:
            for j in range(len(pep) - k_window + 1):
                want9[pep[j : j + k_window]].add(i)
        elif len(pep) == short_window:
            want8[pep].add(i)
        else:
            raise ValueError(
                f"peptide {pep!r} (row {i}) has length {len(pep)}, below the "
                f"{short_window}-mer window; refusing to skip it silently because a "
                "skipped row is indistinguishable from an unmatched one in the "
                "assigned fraction"
            )

    hits: list[set[str]] = [set() for _ in peptides]
    seen_proteins: set[str] = set()
    for path in paths:
        for acc, seq in _parse_fasta(path):
            n = len(seq)
            if want9 and n >= k_window:
                for j in range(n - k_window + 1):
                    rows = want9.get(seq[j : j + k_window])
                    if rows:
                        seen_proteins.add(acc)
                        for r in rows:
                            hits[r].add(acc)
            if want8 and n >= short_window:
                for j in range(n - short_window + 1):
                    rows = want8.get(seq[j : j + short_window])
                    if rows:
                        seen_proteins.add(acc)
                        for r in rows:
                            hits[r].add(acc)

    uf = _UnionFind()
    for acc in sorted(seen_proteins):
        uf.add(acc)
    for h in hits:
        if len(h) > 1:
            ordered = sorted(h)
            first = ordered[0]
            for other in ordered[1:]:
                uf.union(first, other)

    keys: list[str] = []
    n_assigned = 0
    for i, h in enumerate(hits):
        if h:
            n_assigned += 1
            keys.append(f"PROTEIN:{uf.find(sorted(h)[0])}")
        else:
            keys.append(f"SINGLETON:{peptides[i]}")

    cmap = ClusterMap.from_keys(keys, n_assigned=n_assigned)
    cmap.n_proteins = len({k for k in cmap.keys if k.startswith("PROTEIN:")})
    cmap.n_singleton_clusters = len({k for k in cmap.keys if k.startswith("SINGLETON:")})
    cmap.notes = {
        "k_window": k_window,
        "short_window": short_window,
        "proteomes": [p.name for p in paths],
        "n_proteins_hit_before_merge": len(seen_proteins),
    }
    return cmap


def study_level_basis(study_ids, *, label: str = "study") -> str | Unmeasurable:
    """The study-clustering field, which is REQUIRED and may be a sentinel.

    Returns `UNMEASURABLE_STUDY_CONSTANT` on the infectious set, where `study_id` is
    constant across all 2,093 rows. Reported rather than omitted so a reader can see the
    question was asked: study- and leakage-group-disjointness are incompatible here
    (15,032 of 16,021 pathogen rows sit in one connected component).
    """
    if study_ids is None:
        return Unmeasurable("no study_id column was supplied to this estimate")
    ids = list(study_ids)
    if not ids:
        return Unmeasurable("no rows")
    distinct = set(ids)
    if len(distinct) == 1:
        return UNMEASURABLE_STUDY_CONSTANT
    sizes = sorted((sum(1 for x in ids if x == d) for d in distinct), reverse=True)
    if sizes[0] / len(ids) > 0.99:
        return Unmeasurable(
            f"one study holds {sizes[0]} of {len(ids)} rows; the study bootstrap is a "
            "one-cluster bootstrap in all but name"
        )
    return f"{label}/{len(distinct)}-clusters"


# ---------------------------------------------------------------------------
# Draws
# ---------------------------------------------------------------------------


@dataclass
class DrawMatrix:
    """B x K matrix of CLUSTER indices. Rows are expanded at use, not stored."""

    matrix: list[list[int]]
    n_clusters: int
    recipe: SeedRecipe
    unit: str = "cluster"
    provenance: str = "fresh"

    @property
    def B(self) -> int:
        return len(self.matrix)

    def __str__(self) -> str:
        return (
            f"{self.B} x {self.n_clusters} {self.unit} draws  "
            f"{self.provenance}/{self.recipe.digest()}"
        )


def draws(
    n_clusters: int,
    B: int = 10000,
    *,
    seed_recipe: SeedRecipe,
    unit: str = "cluster",
    provenance: str = "fresh",
) -> DrawMatrix:
    """Deterministic cluster-index matrix: B replicates, each K clusters with replacement.

    Stage one of the two-stage draw; `_expand` is stage two. Carrying no labels or scores,
    the matrix can be shared verbatim by both arms of a paired comparison.

    Replicate b is reseeded from `sha256(version|set_id|comparison_id|b)` rather than
    advancing one RNG, so a table recomputed in a different order, or with one cell
    omitted, reproduces byte for byte. The alternative fails invisibly -- every arm still
    looks internally consistent.
    """
    if n_clusters < 1:
        raise UnmeasuredBasis(f"cannot draw from {n_clusters} clusters")
    if B < 1:
        raise ValueError(f"B must be at least 1, got {B}")
    if unit not in ("cluster", "row"):
        raise ValueError(f"unit must be 'cluster' or 'row', got {unit!r}")
    if provenance not in ("fresh", "saved"):
        raise ValueError(f"provenance must be 'fresh' or 'saved', got {provenance!r}")

    matrix: list[list[int]] = []
    for b in range(B):
        rnd = random.Random(seed_recipe.seed(b))
        matrix.append([rnd.randrange(n_clusters) for _ in range(n_clusters)])
    return DrawMatrix(
        matrix=matrix,
        n_clusters=n_clusters,
        recipe=seed_recipe,
        unit=unit,
        provenance=provenance,
    )


# ---------------------------------------------------------------------------
# Estimates
# ---------------------------------------------------------------------------


def _percentile(sorted_vals: list[float], p: float) -> float:
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    pos = p * (n - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return sorted_vals[int(pos)]
    frac = pos - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def _quantile_mc_se(sorted_vals: list[float], p: float) -> float:
    """Monte-Carlo SE of a bootstrap percentile, from the order statistics themselves.

    The p-quantile's rank among B replicates is Binomial(B, p), s.d. sqrt(p(1-p)B) in rank
    units; the local spacing of the order statistics converts that to value units. Reported
    because a bound quoted to four decimals off B = 10,000 replicates is not good to four.
    """
    n = len(sorted_vals)
    if n < 3:
        return float("nan")
    sd_rank = math.sqrt(p * (1.0 - p) * n)
    centre = p * (n - 1)
    lo_i = max(0, min(n - 1, int(round(centre - sd_rank))))
    hi_i = max(0, min(n - 1, int(round(centre + sd_rank))))
    return (sorted_vals[hi_i] - sorted_vals[lo_i]) / 2.0


def _as_float(result) -> float:
    return float(result.value) if hasattr(result, "value") else float(result)


def _metric_arity(fn) -> int:
    try:
        params = inspect.signature(fn).parameters.values()
    except (TypeError, ValueError):
        return 3
    n = sum(
        1
        for p in params
        if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD) and p.default is p.empty
    )
    return 2 if n <= 2 else 3


def _call_metric(fn, arity, labels, scores, groups) -> float:
    if arity == 2:
        return _as_float(fn(labels, scores))
    return _as_float(fn(labels, scores, groups))


def _expand(cmap: ClusterMap, drawn: list[int]) -> list[int]:
    """Stage two: every row of every drawn cluster, so the replicate row count is ragged.

    That raggedness IS the cluster-size variance the design exists to capture; do not
    renormalise to a fixed n, which restores the row-level variance clustering avoids.
    """
    out: list[int] = []
    for c in drawn:
        out.extend(cmap.members[c])
    return out


def _basis(unit: str, B: int, provenance: str, digest: str, assigned) -> str:
    af = str(assigned) if isinstance(assigned, Unmeasurable) else f"{float(assigned):.3f}"
    return f"{unit}/{B}/{provenance}/{digest}/{af}"


@dataclass
class Estimate:
    """A point estimate with its interval AND the basis the interval was computed on.

    Every field is required; `study_level` has no default so the sentinel cannot be
    forgotten.
    """

    value: float
    lo: float
    hi: float
    mc_se: tuple
    basis: str
    study_level: object
    n_rows: int
    n_clusters: int
    assigned_fraction: object
    n_replicates_used: int
    n_replicates_failed: int
    provisional: bool

    def __post_init__(self) -> None:
        if self.study_level is None:
            raise UnmeasuredBasis(
                "study_level is required; pass a basis string or an Unmeasurable sentinel "
                "with a reason. A missing key reads as 'not relevant', which is the one "
                "thing it never means."
            )

    @property
    def width(self) -> float:
        return self.hi - self.lo

    def __str__(self) -> str:
        flag = "  [PROVISIONAL]" if self.provisional else ""
        return (
            f"{self.value:.6f} [{self.lo:.6f}, {self.hi:.6f}]  "
            f"mc_se ({self.mc_se[0]:.6f}, {self.mc_se[1]:.6f})  "
            f"basis {self.basis}{flag}"
        )


def interval(
    metric_fn,
    labels,
    scores,
    groups,
    clusters: ClusterMap,
    draw_matrix: DrawMatrix,
    *,
    alpha: float = 0.05,
    study_ids=None,
    max_failed_fraction: float = 0.5,
) -> Estimate:
    """Percentile interval over the two-stage cluster bootstrap.

    `metric_fn` is called as `fn(labels, scores, groups)`, or `fn(labels, scores)` when it
    takes two arguments, so `stratified_auroc` and `auroc` both drop in.

    Replicates the metric refuses are counted on `n_replicates_failed`, never patched to
    0.5. Above `max_failed_fraction` this raises: an interval over whichever replicates
    happened to be admissible is an interval over a different population.
    """
    labels, scores, groups = list(labels), list(scores), list(groups)
    n = len(labels)
    if not (n == len(scores) == len(groups)):
        raise ValueError("labels, scores and groups must be the same length")
    if clusters.n_rows != n:
        raise ValueError(
            f"cluster map covers {clusters.n_rows} rows, data has {n}"
        )
    if draw_matrix.n_clusters != clusters.n_clusters:
        raise UnmeasuredBasis(
            f"draw matrix was built for {draw_matrix.n_clusters} clusters, the map has "
            f"{clusters.n_clusters}; these draws do not describe this design"
        )

    arity = _metric_arity(metric_fn)
    point = _call_metric(metric_fn, arity, labels, scores, groups)

    reps: list[float] = []
    failed = 0
    for drawn in draw_matrix.matrix:
        idx = _expand(clusters, drawn)
        try:
            reps.append(
                _call_metric(
                    metric_fn,
                    arity,
                    [labels[i] for i in idx],
                    [scores[i] for i in idx],
                    [groups[i] for i in idx],
                )
            )
        except Exception:  # noqa: BLE001 -- the metric's own refusal, recorded not patched
            failed += 1

    B = draw_matrix.B
    if failed > max_failed_fraction * B:
        raise UnmeasuredBasis(
            f"{failed} of {B} replicates were refused by the metric; the surviving "
            "replicates describe a different population from the one being estimated"
        )
    if not reps:
        raise UnmeasuredBasis("every replicate was refused by the metric")

    reps.sort()
    lo_p, hi_p = alpha / 2.0, 1.0 - alpha / 2.0
    return Estimate(
        value=point,
        lo=_percentile(reps, lo_p),
        hi=_percentile(reps, hi_p),
        mc_se=(_quantile_mc_se(reps, lo_p), _quantile_mc_se(reps, hi_p)),
        basis=_basis(
            draw_matrix.unit,
            B,
            draw_matrix.provenance,
            draw_matrix.recipe.digest(),
            clusters.assigned_fraction,
        ),
        study_level=study_level_basis(study_ids),
        n_rows=n,
        n_clusters=clusters.n_clusters,
        assigned_fraction=clusters.assigned_fraction,
        n_replicates_used=len(reps),
        n_replicates_failed=failed,
        provisional=clusters.provisional,
    )


# ---------------------------------------------------------------------------
# Paired comparison
# ---------------------------------------------------------------------------


@dataclass
class Delta:
    """A paired difference, its design effect, and the effective n that follows."""

    d: float
    lo: float
    hi: float
    p: float
    se: float
    design_effect: float
    n_eff: float
    mc_se: tuple
    basis: str
    row_basis: str
    study_level: object
    cluster_width: float
    row_width: float
    n_rows: int
    n_clusters: int
    B: int
    n_replicates_failed: int
    draws: object = None

    def __post_init__(self) -> None:
        if self.study_level is None:
            raise UnmeasuredBasis("study_level is required on a Delta; see Estimate")

    def __str__(self) -> str:
        return (
            f"{self.d:+.6f} [{self.lo:+.6f}, {self.hi:+.6f}]  p={self.p:.4f}  "
            f"se={self.se:.6f}  DE={self.design_effect:.3f}  n_eff={self.n_eff:.0f}  "
            f"basis {self.basis}"
        )


def paired(
    metric_fn,
    labels,
    a_scores,
    b_scores,
    groups,
    clusters: ClusterMap,
    draw_matrix: DrawMatrix,
    *,
    alpha: float = 0.05,
    row_draws: DrawMatrix | None = None,
    study_ids=None,
    max_failed_fraction: float = 0.5,
) -> Delta:
    """Shared-index paired bootstrap: d = metric(A) - metric(B) on ONE draw matrix.

    Both arms score the SAME resampled rows each replicate, so shared row noise cancels.
    Comparing two independent marginal CIs is weaker and gave "+0.0246 over PRIME-2.0,
    p = 0.078" where paired and clustered it is +0.0187 [-0.0111, +0.0466], p = 0.215.

    `p = (#{delta <= 0} + 1) / (B + 1)`, one-sided. The +1s give the honest floor 1/(B+1);
    without them a comparison no replicate crossed prints p = 0.0000, which B = 10,000
    cannot support -- this project retracted a p < 0.0005 for exactly that.

    `design_effect = (cluster width / row width)^2`, measured on this comparison with the
    row arm resampled from a sibling recipe, never read from a table: the widening is
    model-dependent (1.00 for one published column, 1.05 for another on the same 833
    rows). `n_eff = n / DE`, and DE below 1 is reported as measured, not clamped.
    """
    labels = list(labels)
    a_scores, b_scores, groups = list(a_scores), list(b_scores), list(groups)
    n = len(labels)
    if not (n == len(a_scores) == len(b_scores) == len(groups)):
        raise ValueError("labels, both score arms and groups must be the same length")
    if clusters.n_rows != n:
        raise ValueError(f"cluster map covers {clusters.n_rows} rows, data has {n}")
    if draw_matrix.unit != "cluster":
        raise UnmeasuredBasis(
            f"paired() expects a cluster draw matrix, got unit={draw_matrix.unit!r}"
        )
    if draw_matrix.n_clusters != clusters.n_clusters:
        raise UnmeasuredBasis(
            f"draw matrix was built for {draw_matrix.n_clusters} clusters, the map has "
            f"{clusters.n_clusters}"
        )

    arity = _metric_arity(metric_fn)
    d_point = _call_metric(metric_fn, arity, labels, a_scores, groups) - _call_metric(
        metric_fn, arity, labels, b_scores, groups
    )

    def _deltas(cmap: ClusterMap, dm: DrawMatrix) -> tuple[list[float], int]:
        out: list[float] = []
        bad = 0
        for drawn in dm.matrix:
            idx = _expand(cmap, drawn)
            y = [labels[i] for i in idx]
            g = [groups[i] for i in idx]
            try:
                out.append(
                    _call_metric(metric_fn, arity, y, [a_scores[i] for i in idx], g)
                    - _call_metric(metric_fn, arity, y, [b_scores[i] for i in idx], g)
                )
            except Exception:  # noqa: BLE001
                bad += 1
        return out, bad

    cluster_reps, failed = _deltas(clusters, draw_matrix)
    if not cluster_reps:
        raise UnmeasuredBasis("every clustered replicate was refused by the metric")
    # Without this guard, with 477 of 600 replicates refused (79.5%), `paired` returned
    # p = 0.113 where the survivor-consistent value is 0.548: p divides by B, not by the
    # survivors, so every refused replicate silently votes "delta > 0".
    if failed > max_failed_fraction * draw_matrix.B:
        raise UnmeasuredBasis(
            f"{failed} of {draw_matrix.B} clustered replicates were refused by the "
            f"metric; the surviving {len(cluster_reps)} are not a sample of the design "
            "and no delta, interval or p-value is reportable from them"
        )

    if row_draws is None:
        row_draws = draws(
            n,
            B=draw_matrix.B,
            seed_recipe=draw_matrix.recipe.derive("row"),
            unit="row",
            provenance=draw_matrix.provenance,
        )
    row_map = ClusterMap.rows_as_clusters(n)
    if row_draws.n_clusters != n:
        raise UnmeasuredBasis(
            f"row draw matrix covers {row_draws.n_clusters} units, data has {n} rows"
        )
    row_reps, _ = _deltas(row_map, row_draws)
    if not row_reps:
        raise UnmeasuredBasis("every row-level replicate was refused by the metric")

    cs = sorted(cluster_reps)
    rs = sorted(row_reps)
    lo_p, hi_p = alpha / 2.0, 1.0 - alpha / 2.0
    lo, hi = _percentile(cs, lo_p), _percentile(cs, hi_p)
    cluster_width = hi - lo
    row_width = _percentile(rs, hi_p) - _percentile(rs, lo_p)

    if row_width <= 0.0:
        raise UnmeasuredBasis(
            "the row-level interval has zero width; the design effect is undefined and "
            "cannot be defaulted to 1.0"
        )
    de = (cluster_width / row_width) ** 2

    mean = sum(cs) / len(cs)
    se = math.sqrt(sum((x - mean) ** 2 for x in cs) / (len(cs) - 1)) if len(cs) > 1 else 0.0
    B = draw_matrix.B
    p = (sum(1 for x in cluster_reps if x <= 0.0) + 1) / (B + 1)

    return Delta(
        d=d_point,
        lo=lo,
        hi=hi,
        p=p,
        se=se,
        design_effect=de,
        n_eff=n / de if de > 0 else float("inf"),
        mc_se=(_quantile_mc_se(cs, lo_p), _quantile_mc_se(cs, hi_p)),
        basis=_basis(
            "cluster", B, draw_matrix.provenance, draw_matrix.recipe.digest(),
            clusters.assigned_fraction,
        ),
        row_basis=_basis(
            "row", row_draws.B, row_draws.provenance, row_draws.recipe.digest(), 1.0
        ),
        study_level=study_level_basis(study_ids),
        cluster_width=cluster_width,
        row_width=row_width,
        n_rows=n,
        n_clusters=clusters.n_clusters,
        B=B,
        n_replicates_failed=failed,
        draws=draw_matrix,
    )


# ---------------------------------------------------------------------------
# The power gate
# ---------------------------------------------------------------------------


def z_sum(alpha: float = 0.025, power: float = 0.80) -> float:
    """`z_{1-alpha} + z_{1-power}`; 2.801585 at alpha=0.025, power=0.80.

    Recomputed rather than copied from `power_assessment.json`, so the gate cannot drift
    from the filed assessment undetected.
    """
    nd = NormalDist()
    return nd.inv_cdf(1.0 - alpha) + nd.inv_cdf(power)


@dataclass
class MDEReport:
    mde: float
    se_delta: float
    z_sum: float
    alpha: float
    power: float
    n_eff: float
    design_effect: float
    declared_threshold: object
    basis: str
    study_level: object
    verdict: str

    def __str__(self) -> str:
        decl = (
            "none declared"
            if self.declared_threshold is None
            else f"declared {self.declared_threshold:+.4f}"
        )
        return (
            f"MDE {self.mde:.6f} at alpha={self.alpha} power={self.power} "
            f"(z_sum {self.z_sum:.6f}, se {self.se_delta:.6f})  "
            f"n_eff {self.n_eff:.0f}  DE {self.design_effect:.3f}  "
            f"{decl}  {self.verdict}  basis {self.basis}"
        )


def mde(
    delta: Delta,
    *,
    alpha: float = 0.025,
    power: float = 0.80,
    declared_threshold: float | None = None,
) -> MDEReport:
    """Minimum detectable paired delta `= (z_{1-alpha} + z_{1-power}) * se_delta`.

    Raises `UnmeasuredBasis` without cluster draws (a row-level SE understates the
    requirement by the design effect), and `PreRegistrationUnderpowered` when
    `abs(declared_threshold) < mde` -- never a warning. Worked instance: the +0.0187
    margin over PRIME-2.0 sits against a row-level SE of ~0.0144, so even the row-level
    MDE is 2.801585 * 0.0144 = 0.0403 and the clustered one is larger.
    """
    if delta.draws is None:
        raise UnmeasuredBasis(
            "this Delta carries no draw matrix; the MDE would be computed from an SE "
            "whose resampling unit is unknown"
        )
    if getattr(delta.draws, "unit", None) != "cluster":
        raise UnmeasuredBasis(
            f"cluster draws do not exist for this comparison (unit="
            f"{getattr(delta.draws, 'unit', None)!r}); a row-level MDE understates the "
            "requirement by the design effect"
        )
    if not (delta.se > 0.0) or not math.isfinite(delta.se):
        raise UnmeasuredBasis(f"delta SE is {delta.se!r}; the MDE is undefined")

    zs = z_sum(alpha, power)
    value = zs * delta.se

    verdict = "no threshold declared"
    if declared_threshold is not None:
        if abs(declared_threshold) < value:
            raise PreRegistrationUnderpowered(
                f"declared threshold {declared_threshold:+.4f} is below the minimum "
                f"detectable paired delta {value:.4f} "
                f"(z_sum {zs:.6f} x se {delta.se:.6f}, basis {delta.basis}, "
                f"n_eff {delta.n_eff:.0f} after design effect {delta.design_effect:.3f}). "
                "A success criterion below the MDE is not pre-registerable: the design "
                "cannot distinguish it from zero, so 'success' and 'no effect' are the "
                "same outcome."
            )
        verdict = f"declared {declared_threshold:+.4f} clears the MDE"

    return MDEReport(
        mde=value,
        se_delta=delta.se,
        z_sum=zs,
        alpha=alpha,
        power=power,
        n_eff=delta.n_eff,
        design_effect=delta.design_effect,
        declared_threshold=declared_threshold,
        basis=delta.basis,
        study_level=delta.study_level,
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# The mixed-basis guard
# ---------------------------------------------------------------------------


def format_table(cells, *, title: str | None = None) -> str:
    """Format a table of Estimates, or refuse it.

    `cells` is a mapping or a sequence of (name, Estimate) pairs. Every cell must carry
    the same basis string -- the whole string, since B and the assigned fraction move
    widths too -- and the first disagreement raises `MixedBasisError`. Raised rather than
    footnoted: at 1.00x for one published column and 1.05x for another on the same 833
    rows, a row-level interval beside a clustered one does not look wrong, but reorders
    the ranking by the size of the margins being ranked. The study-level field prints for
    every row, sentinel included.
    """
    items = list(cells.items()) if hasattr(cells, "items") else [tuple(x) for x in cells]
    if not items:
        raise ValueError("no cells to format; an empty table cannot be checked")
    for name, est in items:
        if not isinstance(est, Estimate):
            raise TypeError(f"cell {name!r} is {type(est).__name__}, expected Estimate")

    ref_name, ref = items[0]
    for name, est in items[1:]:
        if est.basis != ref.basis:
            raise MixedBasisError(
                f"cell {name!r} has basis {est.basis!r} but cell {ref_name!r} has "
                f"{ref.basis!r}. These intervals are not comparable: the widening under "
                "clustering is model-dependent, so the mixed table silently reorders the "
                "ranking. Recompute every cell on one basis, or print two tables."
            )

    width = max(len(str(n)) for n, _ in items)
    lines: list[str] = []
    if title:
        lines.append(title)
    lines.append(f"basis: {ref.basis}")
    for name, est in items:
        study = (
            est.study_level.explain()
            if isinstance(est.study_level, Unmeasurable)
            else str(est.study_level)
        )
        lines.append(
            f"  {str(name):<{width}}  {est.value:.4f} [{est.lo:.4f}, {est.hi:.4f}]  "
            f"+/-mc {max(est.mc_se):.4f}  study: {study}"
        )
    return "\n".join(lines)
