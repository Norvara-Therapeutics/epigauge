"""Metrics. Standard library only, so this runs on the corpus-building machine.

`stratified_auroc` is the primary comparability estimator. `pooled_rank_auroc` is the
project's historical `within_allele_auroc`, retained to reproduce published figures and
diagnostic-only: it rank-normalises within group then pools, so cross-group pairs still
contribute and it admits a mixture effect -- a pooled 0.5835 came out above both its
constituents (pathogen 0.5811, self 0.5730), and singletons are injected at 0.5.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

__all__ = [
    "DegenerateGroupingError",
    "auroc",
    "stratified_auroc",
    "pooled_rank_auroc",
    "mixture_decomposition",
    "StratifiedResult",
    "MixtureResult",
]


class DegenerateGroupingError(ValueError):
    """The grouping cannot support a stratified estimate. Raised, never defaulted."""


def _midranks(scores: list[float]) -> list[float]:
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        shared = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def auroc(labels, scores) -> float:
    """Tie-corrected Mann-Whitney. A constant score gives exactly 0.5, never 0 or 1."""
    labels, scores = list(labels), list(scores)
    if len(labels) != len(scores):
        raise ValueError(f"length mismatch: {len(labels)} labels, {len(scores)} scores")
    pos = sum(1 for y in labels if y == 1)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        raise DegenerateGroupingError(f"AUROC undefined: {pos} positive, {neg} negative")
    ranks = _midranks(scores)
    rsum = sum(r for r, y in zip(ranks, labels) if y == 1)
    return (rsum - pos * (pos + 1) / 2.0) / (pos * neg)


@dataclass
class StratifiedResult:
    value: float
    n_strata_contributing: int
    n_strata_total: int
    n_rows_evaluable: int
    n_rows_total: int
    herfindahl: float
    top3_pair_weight: float
    provisional: bool
    per_stratum: dict = field(default_factory=dict)

    @property
    def evaluable_fraction(self) -> float:
        return self.n_rows_evaluable / self.n_rows_total if self.n_rows_total else 0.0

    def __str__(self) -> str:
        flag = "  [PROVISIONAL: top-3 strata hold >50% of pair weight]" if self.provisional else ""
        return (
            f"{self.value:.6f}  "
            f"strata {self.n_strata_contributing}/{self.n_strata_total}  "
            f"rows {self.n_rows_evaluable}/{self.n_rows_total}  "
            f"H={self.herfindahl:.4f}{flag}"
        )


def stratified_auroc(labels, scores, groups, *, min_per_class: int = 1) -> StratifiedResult:
    """Pair-weighted within-group AUROC. The primary comparability estimator.

    A_wa = sum_a w_a * C_a,  w_a = n+_a * n-_a / sum_b n+_b * n-_b

    Pair-weighted rather than unweighted macro: three alleles hold 49% of the evaluable
    rows on the infectious set and 42 of 78 contribute nothing, so a macro would give a
    2-row allele the weight of a 300-row one. Strata with fewer than `min_per_class` of
    either class get weight ZERO; they are never mapped to 0.5.
    """
    labels, scores, groups = list(labels), list(scores), list(groups)
    n = len(labels)
    if not (n == len(scores) == len(groups)):
        raise ValueError("labels, scores and groups must be the same length")
    if n == 0:
        raise DegenerateGroupingError("no rows")

    by_group: dict[object, list[int]] = defaultdict(list)
    for i, g in enumerate(groups):
        if g is None or g == "":
            raise DegenerateGroupingError(
                f"row {i} has an empty group key. A uniformly-missing grouping silently "
                "returns the pooled statistic -- refusing rather than reporting it."
            )
        by_group[g].append(i)

    # Raises rather than returning a number: a uniformly-missing allele column makes
    # within-allele equal pooled exactly (0.584138 vs the correct 0.492317), and
    # all-singleton grouping returns exactly 0.5 -- both look publishable.
    if len(by_group) == 1:
        raise DegenerateGroupingError(
            f"all {n} rows share one group key ({next(iter(by_group))!r}); the stratified "
            "estimate would equal pooled AUROC exactly"
        )
    if all(len(ix) == 1 for ix in by_group.values()):
        raise DegenerateGroupingError(
            f"all {len(by_group)} groups are singletons; the stratified estimate would be "
            "exactly 0.5 regardless of the scores"
        )

    per: dict[object, dict] = {}
    pair_counts: dict[object, int] = {}
    for g, ix in by_group.items():
        pos = sum(1 for i in ix if labels[i] == 1)
        neg = len(ix) - pos
        admissible = pos >= min_per_class and neg >= min_per_class
        entry = {"n": len(ix), "n_pos": pos, "n_neg": neg, "admissible": admissible,
                 "auroc": None, "pair_weight": 0.0}
        if admissible:
            entry["auroc"] = auroc([labels[i] for i in ix], [scores[i] for i in ix])
            pair_counts[g] = pos * neg
        per[g] = entry

    total_pairs = sum(pair_counts.values())
    if total_pairs == 0:
        raise DegenerateGroupingError(
            f"no group has at least {min_per_class} of each class across "
            f"{len(by_group)} groups and {n} rows"
        )

    value = 0.0
    for g, pairs in pair_counts.items():
        w = pairs / total_pairs
        per[g]["pair_weight"] = w
        value += w * per[g]["auroc"]

    weights = sorted((e["pair_weight"] for e in per.values()), reverse=True)
    herf = sum(w * w for w in weights)
    top3 = sum(weights[:3])
    rows_eval = sum(e["n"] for e in per.values() if e["admissible"])

    return StratifiedResult(
        value=value,
        n_strata_contributing=len(pair_counts),
        n_strata_total=len(by_group),
        n_rows_evaluable=rows_eval,
        n_rows_total=n,
        herfindahl=herf,
        top3_pair_weight=top3,
        provisional=top3 > 0.5,
        per_stratum=per,
    )


def pooled_rank_auroc(labels, scores, groups) -> float:
    """DIAGNOSTIC ONLY. Rank-normalise within group, then pool.

    The historical `within_allele_auroc`, kept to reproduce published figures. Never report
    it as a within-allele number: cross-group pairs contribute.
    """
    labels, scores, groups = list(labels), list(scores), list(groups)
    by_group: dict[object, list[int]] = defaultdict(list)
    for i, g in enumerate(groups):
        by_group[g].append(i)
    normed = [0.0] * len(labels)
    for ix in by_group.values():
        if len(ix) == 1:
            normed[ix[0]] = 0.5
            continue
        r = _midranks([scores[i] for i in ix])
        for k, i in enumerate(ix):
            normed[i] = (r[k] - 0.5) / len(ix)
    return auroc(labels, normed)


@dataclass
class MixtureResult:
    pooled: float
    diagonal_mass: float
    diagonal_only: float
    off_diagonal_mean: float
    max_constituent: float
    excess: float
    artifact: bool
    blocks: dict = field(default_factory=dict)
    cross: dict = field(default_factory=dict)

    def __str__(self) -> str:
        tag = "  <-- MIXTURE ARTIFACT" if self.artifact else ""
        return (
            f"pooled {self.pooled:.6f}  max constituent {self.max_constituent:.6f}  "
            f"excess {self.excess:+.6f}  diag mass {self.diagonal_mass:.4f}{tag}"
        )


def mixture_decomposition(labels, scores, blocks) -> MixtureResult:
    """Exact algebraic decomposition of pooled AUROC over a partition.

    A_pool = sum_{j,k} w_jk C_jk ,  w_jk = n+_j n-_k / (N+ N-)

    `excess = A_pool - max_j C_jj`, flagged against the max rather than a weighted average.
    A positive excess is a composition effect, not a peptide-level gain: measured at
    pathogen 0.5811, self 0.5730, pooled 0.5835, excess +0.0024.
    """
    labels, scores, blocks = list(labels), list(scores), list(blocks)
    pos_by: dict[object, list[int]] = defaultdict(list)
    neg_by: dict[object, list[int]] = defaultdict(list)
    for i, (y, b) in enumerate(zip(labels, blocks)):
        (pos_by if y == 1 else neg_by)[b].append(i)

    n_pos = sum(len(v) for v in pos_by.values())
    n_neg = sum(len(v) for v in neg_by.values())
    if n_pos == 0 or n_neg == 0:
        raise DegenerateGroupingError(f"decomposition undefined: {n_pos} pos, {n_neg} neg")

    cross: dict[tuple, dict] = {}
    pooled = 0.0
    diag_mass = 0.0
    for j, pix in pos_by.items():
        for k, nix in neg_by.items():
            w = (len(pix) * len(nix)) / (n_pos * n_neg)
            wins = ties = 0
            for a in pix:
                sa = scores[a]
                for b in nix:
                    if sa > scores[b]:
                        wins += 1
                    elif sa == scores[b]:
                        ties += 1
            c = (wins + 0.5 * ties) / (len(pix) * len(nix))
            cross[(j, k)] = {"weight": w, "c": c}
            pooled += w * c
            if j == k:
                diag_mass += w

    diag = {j: cross[(j, j)]["c"] for j in pos_by if (j, j) in cross}
    if not diag:
        raise DegenerateGroupingError("no block carries both classes; diagonal is empty")
    diag_only = (
        sum(cross[(j, j)]["weight"] * c for j, c in diag.items()) / diag_mass
        if diag_mass > 0 else float("nan")
    )
    off = [v["c"] for kk, v in cross.items() if kk[0] != kk[1]]
    max_c = max(diag.values())
    excess = pooled - max_c

    return MixtureResult(
        pooled=pooled,
        diagonal_mass=diag_mass,
        diagonal_only=diag_only,
        off_diagonal_mean=(sum(off) / len(off)) if off else float("nan"),
        max_constituent=max_c,
        excess=excess,
        artifact=excess > 0,
        blocks={j: {"auroc": c, "weight": cross[(j, j)]["weight"]} for j, c in diag.items()},
        cross={f"{a}|{b}": v for (a, b), v in cross.items()},
    )
