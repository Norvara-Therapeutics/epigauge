"""Trivial baselines, shipped as leaderboard rows rather than footnotes.

Measured figures on the 2,093-row infectious benchmark: see `MEASURED_INFECTIOUS`.
Standard library only. OLS is solved by hand through the normal equations (3-4
predictors), and a singular design raises rather than being ridged into a number.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from epigauge.canon import canonicalise_allele, canonicalise_peptide

# One implementation of mid-ranks per package: four copies of allele normalisation each
# produced a wrong published number (canon.py).
from epigauge.metrics import _midranks as midranks

__all__ = [
    "BaselineError",
    "EvaluationLeakError",
    "SingularDesignError",
    "KD",
    "PWM_LENGTHS",
    "MEASURED_INFECTIOUS",
    "Baseline",
    "AlleleOnly",
    "Hydropathy",
    "HydropathyLengthAllele",
    "PWM",
    "Length",
    "ShuffledSentinel",
    "SentinelDistribution",
    "sentinel_distribution",
    "within_allele_rank",
    "residualise",
    "mean_kd",
    "mean_kd_non_anchor",
    "ols",
]


class BaselineError(ValueError):
    """A baseline was asked for something it cannot honestly produce."""


class EvaluationLeakError(BaselineError):
    """A fitting call was handed a row from the declared evaluation set."""


class SingularDesignError(BaselineError):
    """The OLS design is rank-deficient. Raised rather than regularised silently."""


# Kyte-Doolittle hydropathy. Kyte J, Doolittle RF, J Mol Biol 157:105-132 (1982).
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Class-I lengths covering 100% of both benchmarks; PWM refuses anything outside.
PWM_LENGTHS = (8, 9, 10, 11)

# Shrinkage for the allele target encoding, fixed at the value behind every measured
# figure here. Not a knob to tune against the evaluation set.
ALLELE_SHRINKAGE = 10.0

# Measured on the 2,093-row infectious benchmark, not computed here: a reimplementation
# that lands far from these is visibly wrong. Cite with the corpus fingerprint.
MEASURED_INFECTIOUS = {
    "HydropathyLengthAllele": {"within": 0.5514, "pooled": 0.7084},
    "Hydropathy": {"within": 0.5264, "pooled": None},
    # PWM: no archived artifact for its within-allele figure, so none is quoted.
    "PWM": {"within": None, "pooled": None},
    "AlleleOnly": {"within": 0.5000, "pooled": 0.7122},
    "Length": {"within": 0.4132, "pooled": None},
    "ShuffledSentinel[order=before]": {"within": 0.4751, "pooled": 0.4448},
    "ShuffledSentinel[order=after]": {"within": 0.5120, "pooled": 0.6083},
}


# ---------------------------------------------------------------------------
# Row access. One place for a malformed row to fail.
# ---------------------------------------------------------------------------

def _field(row, key):
    try:
        return row[key]
    except (KeyError, TypeError) as e:
        raise BaselineError(f"row is missing required field {key!r}: {row!r}") from e


def _peptide(row) -> str:
    return canonicalise_peptide(_field(row, "peptide"))


def _allele(row) -> str:
    return canonicalise_allele(_field(row, "mhc_allele"))


def _label(row) -> int:
    y = _field(row, "label")
    try:
        f = float(y)
    except (TypeError, ValueError) as e:
        raise BaselineError(f"label {y!r} is not numeric") from e
    if f not in (0.0, 1.0):
        raise BaselineError(f"label must be 0 or 1, got {y!r}")
    return int(f)


def _record_id(row):
    rid = _field(row, "record_id")
    if rid is None or rid == "":
        raise BaselineError(f"row has an empty record_id: {row!r}")
    return rid


def mean_kd(peptide: str) -> float:
    """Mean Kyte-Doolittle hydropathy over the whole peptide."""
    p = canonicalise_peptide(peptide)
    return sum(KD[a] for a in p) / len(p)


def mean_kd_non_anchor(peptide: str) -> float:
    """Mean KD over the TCR-facing residues: everything but P2 and the C-terminus.

    P2 and P-omega are the class-I anchors, read by the MHC rather than the receptor.
    Pairing this with whole-peptide KD measures 0.5514 against 0.5264 for KD alone.
    """
    p = canonicalise_peptide(peptide)
    anchors = {1, len(p) - 1}
    rest = [KD[a] for i, a in enumerate(p) if i not in anchors]
    if not rest:
        raise BaselineError(f"peptide {p!r} has no non-anchor positions")
    return sum(rest) / len(rest)


# ---------------------------------------------------------------------------
# OLS by hand. Normal equations plus Gaussian elimination with partial pivoting.
# ---------------------------------------------------------------------------

def ols(design, target, *, tol: float = 1e-10) -> list[float]:
    """Least squares by the normal equations. Raises on a rank-deficient design.

    A singular or near-singular design raises `SingularDesignError` rather than being
    dropped or ridged: a constant predictor (all peptides one length, say) is a fact
    about the data, and residualising past it reports a removal that never happened.
    """
    X = [list(map(float, r)) for r in design]
    y = [float(v) for v in target]
    n = len(X)
    if n != len(y):
        raise BaselineError(f"design has {n} rows, target has {len(y)}")
    if n == 0:
        raise SingularDesignError("empty design")
    p = len(X[0])
    if any(len(r) != p for r in X):
        raise BaselineError("design rows are not all the same width")
    if n < p:
        raise SingularDesignError(f"{n} rows cannot determine {p} coefficients")

    for j in range(p):
        col = [r[j] for r in X]
        if max(col) == min(col) and j != 0:
            raise SingularDesignError(
                f"predictor {j} is constant at {col[0]!r}; it is collinear with the "
                "intercept, so it carries nothing to remove. Drop it deliberately "
                "rather than having it dropped for you."
            )

    A = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    rhs = [sum(X[i][a] * y[i] for i in range(n)) for a in range(p)]
    scale = max((abs(v) for row in A for v in row), default=0.0)
    if scale == 0.0:
        raise SingularDesignError("design is identically zero")

    aug = [A[i] + [rhs[i]] for i in range(p)]
    for c in range(p):
        pivot = max(range(c, p), key=lambda r: abs(aug[r][c]))
        if abs(aug[pivot][c]) <= tol * scale:
            raise SingularDesignError(
                f"design is rank-deficient at column {c} "
                f"(pivot {aug[pivot][c]:.3e} against scale {scale:.3e})"
            )
        aug[c], aug[pivot] = aug[pivot], aug[c]
        for r in range(p):
            if r == c:
                continue
            f = aug[r][c] / aug[c][c]
            if f == 0.0:
                continue
            for k in range(c, p + 1):
                aug[r][k] -= f * aug[c][k]
    return [aug[i][p] / aug[i][i] for i in range(p)]


# ---------------------------------------------------------------------------
# The baseline protocol.
# ---------------------------------------------------------------------------

class Baseline:
    """`fit(rows) -> self`, `score(rows) -> list[float]`, and a leak guard on `fit`.

    `eval_ids` is required: the record_ids `fit` may not learn from. `eval_ids=()`
    disables the guard and must be typed by hand. `_fit_encoders` (label-derived
    features, e.g. the allele target encoding) is split from `_fit_weights` so
    `ShuffledSentinel` can permute at that boundary -- worth 0.4448 against 0.6083
    pooled, so it has to be nameable.
    """

    name = "baseline"
    uses_labels = True
    uses_peptide = True

    def __init__(self, *, eval_ids):
        if eval_ids is None:
            raise BaselineError(
                "eval_ids is required. Pass the evaluation record_ids, or pass "
                "eval_ids=() to declare in writing that no evaluation set is held out."
            )
        self.eval_ids = frozenset(eval_ids)
        self.is_fitted = False

    # -- guard ------------------------------------------------------------
    def check_fitting_pool(self, rows) -> list:
        """Raise if the fitting pool intersects the declared evaluation set."""
        rows = list(rows)
        if not rows:
            raise BaselineError(f"{self.name}: empty fitting pool")
        hits = []
        for r in rows:
            rid = _record_id(r)
            if rid in self.eval_ids:
                hits.append(rid)
        if hits:
            shown = sorted(map(str, set(hits)))[:5]
            raise EvaluationLeakError(
                f"{self.name}: fit() was handed {len(hits)} row(s) whose record_id is in "
                f"the declared evaluation set ({', '.join(shown)}"
                f"{', ...' if len(set(hits)) > 5 else ''}). The fitting pool and the "
                "evaluation set must be disjoint by record_id, not by peptide and not "
                "by allele."
            )
        return rows

    # -- fit --------------------------------------------------------------
    def fit(self, rows):
        rows = self.check_fitting_pool(rows)
        self._fit_encoders(rows)
        self._fit_weights(rows)
        self.is_fitted = True
        return self

    def _fit_encoders(self, rows) -> None:
        """Fit label-derived features. Override where a target encoding exists."""

    def _fit_weights(self, rows) -> None:
        """Fit everything downstream of the encoders."""

    # -- score ------------------------------------------------------------
    def score(self, rows) -> list[float]:
        raise NotImplementedError

    def __repr__(self) -> str:
        return f"<{type(self).__name__} name={self.name!r} fitted={self.is_fitted}>"


class _AlleleEncoder:
    """Shrunk per-allele positive rate: rate(a) = (pos_a + k * global) / (n_a + k).

    Unseen alleles fall back to the global rate, not 0.5, and are counted in
    `n_fallback` rather than passing silently.
    """

    def __init__(self, k: float = ALLELE_SHRINKAGE):
        self.k = float(k)
        self.rate: dict[str, float] = {}
        self.n: dict[str, int] = {}
        self.global_rate: float | None = None
        self.n_fallback = 0

    def fit(self, rows) -> "_AlleleEncoder":
        pos: dict[str, int] = defaultdict(int)
        tot: dict[str, int] = defaultdict(int)
        gp = gt = 0
        for r in rows:
            a = _allele(r)
            y = _label(r)
            pos[a] += y
            tot[a] += 1
            gp += y
            gt += 1
        if gt == 0:
            raise BaselineError("allele encoder: no rows to fit")
        self.global_rate = gp / gt
        self.rate = {a: (pos[a] + self.k * self.global_rate) / (tot[a] + self.k)
                     for a in tot}
        self.n = dict(tot)
        self.n_fallback = 0
        return self

    def lookup(self, allele: str) -> float:
        if self.global_rate is None:
            raise BaselineError("allele encoder used before it was fitted")
        if allele in self.rate:
            return self.rate[allele]
        self.n_fallback += 1
        return self.global_rate


class AlleleOnly(Baseline):
    """Per-allele positive rate from the fitting pool. No peptide input at all.

    Within-allele AUROC is exactly 0.5000 by construction -- the score is constant
    inside a stratum -- against ~0.70 pooled on the infectious benchmark. A table that
    has never seen a peptide reaches 0.70, so the pooled column must not be ranked on.
    """

    name = "AlleleOnly"
    uses_peptide = False

    def __init__(self, *, eval_ids, k: float = ALLELE_SHRINKAGE):
        super().__init__(eval_ids=eval_ids)
        self.encoder = _AlleleEncoder(k)

    def _fit_encoders(self, rows) -> None:
        self.encoder.fit(rows)

    def score(self, rows) -> list[float]:
        if not self.is_fitted:
            raise BaselineError(f"{self.name}: score() before fit()")
        return [self.encoder.lookup(_allele(r)) for r in rows]


class Hydropathy(Baseline):
    """Mean Kyte-Doolittle hydropathy. Unfitted: no labels are read, ever.

    0.5264 within-allele, which outranks most published tools in the same table. `fit`
    only runs the leak guard; it cannot change a score.
    """

    name = "Hydropathy"
    uses_labels = False

    def fit(self, rows):
        # Label-free, but still held to the disjointness contract.
        self.check_fitting_pool(rows)
        self.is_fitted = True
        return self

    def score(self, rows) -> list[float]:
        return [mean_kd(_peptide(r)) for r in rows]


class Length(Baseline):
    """Peptide length alone. Measured 0.4132 within-allele -- strongly anti-predictive.

    Kept in the table because 0.41 is as far from chance as 0.59.
    """

    name = "Length"
    uses_labels = False

    def fit(self, rows):
        self.check_fitting_pool(rows)
        self.is_fitted = True
        return self

    def score(self, rows) -> list[float]:
        return [float(len(_peptide(r))) for r in rows]


class PWM(Baseline):
    """Per-length, per-position add-one log-odds over lengths 8-11.

    20 residues * (8+9+10+11) positions = 760 parameters. Not conditioned on the allele,
    so it floors peptide-level signal rather than mixing the two. An unfitted length
    raises instead of scoring a plausible-looking 0.0.
    """

    name = "PWM"

    def __init__(self, *, eval_ids, lengths=PWM_LENGTHS, min_rows: int = 500,
                 strict: bool = True):
        super().__init__(eval_ids=eval_ids)
        self.lengths = tuple(lengths)
        self.min_rows = int(min_rows)
        self.strict = bool(strict)
        self.matrices: dict[int, list[dict[str, float]]] = {}
        self.skipped: dict[int, str] = {}
        self.n_unscorable = 0

    def _fit_weights(self, rows) -> None:
        pos: dict[int, list[Counter]] = {}
        neg: dict[int, list[Counter]] = {}
        total: Counter = Counter()
        positives: Counter = Counter()
        for r in rows:
            p = _peptide(r)
            y = _label(r)
            L = len(p)
            if L not in self.lengths:
                continue
            if L not in pos:
                pos[L] = [Counter() for _ in range(L)]
                neg[L] = [Counter() for _ in range(L)]
            tgt = pos if y == 1 else neg
            for i, a in enumerate(p):
                tgt[L][i][a] += 1
            total[L] += 1
            positives[L] += y

        self.matrices, self.skipped = {}, {}
        for L in self.lengths:
            n = total[L]
            if n < self.min_rows:
                self.skipped[L] = f"{n} rows below min_rows={self.min_rows}"
                continue
            if not 0 < positives[L] < n:
                self.skipped[L] = f"single-class at length {L} ({positives[L]}/{n} positive)"
                continue
            mat = []
            for i in range(L):
                pt = sum(pos[L][i].values()) + len(AMINO_ACIDS)
                nt = sum(neg[L][i].values()) + len(AMINO_ACIDS)
                mat.append({a: math.log(((pos[L][i][a] + 1.0) / pt)
                                        / ((neg[L][i][a] + 1.0) / nt))
                            for a in AMINO_ACIDS})
            self.matrices[L] = mat
        if not self.matrices:
            raise BaselineError(
                f"{self.name}: no length in {self.lengths} cleared min_rows="
                f"{self.min_rows} with both classes present: {self.skipped}"
            )

    @property
    def n_parameters(self) -> int:
        """20 * sum of fitted lengths. 760 when all of 8-11 are fitted."""
        return sum(len(AMINO_ACIDS) * L for L in self.matrices)

    def score(self, rows) -> list[float]:
        if not self.is_fitted:
            raise BaselineError(f"{self.name}: score() before fit()")
        out = []
        for r in rows:
            p = _peptide(r)
            mat = self.matrices.get(len(p))
            if mat is None:
                self.n_unscorable += 1
                if self.strict:
                    raise BaselineError(
                        f"{self.name}: no matrix for length {len(p)} (peptide {p!r}); "
                        f"fitted lengths are {sorted(self.matrices)}, skipped "
                        f"{self.skipped}. Scoring it 0.0 would put a constant in a "
                        "ranked column; exclude the row deliberately instead."
                    )
                out.append(float("nan"))
                continue
            out.append(sum(mat[i][a] for i, a in enumerate(p)))
        return out


class HydropathyLengthAllele(Baseline):
    """Four inputs: mean KD, non-anchor mean KD, length, allele rate. Measured 0.5514
    within-allele / 0.7084 pooled.

        score(row) = rate(allele) + b0 + b1*meanKD + b2*meanKD_nonanchor + b3*length
        b fitted by OLS of (label - rate(allele)) on [1, meanKD, meanKD_nonanchor, length]

    The allele rate enters as an offset with a fixed unit coefficient: it estimates
    P(y=1 | allele), so refitting its sign discards the one thing it knows. 4 predictors,
    5 coefficients, plus one rate per allele in the fitting pool.
    """

    name = "HydropathyLengthAllele"

    def __init__(self, *, eval_ids, k: float = ALLELE_SHRINKAGE):
        super().__init__(eval_ids=eval_ids)
        self.encoder = _AlleleEncoder(k)
        self.coefficients: list[float] | None = None

    PREDICTORS = ("intercept", "mean_kd", "mean_kd_non_anchor", "length")

    def _features(self, row) -> list[float]:
        p = _peptide(row)
        return [1.0, mean_kd(p), mean_kd_non_anchor(p), float(len(p))]

    def _fit_encoders(self, rows) -> None:
        self.encoder.fit(rows)

    def _fit_weights(self, rows) -> None:
        X = [self._features(r) for r in rows]
        y = [_label(r) - self.encoder.lookup(_allele(r)) for r in rows]
        self.coefficients = ols(X, y)

    def score(self, rows) -> list[float]:
        if not self.is_fitted or self.coefficients is None:
            raise BaselineError(f"{self.name}: score() before fit()")
        b = self.coefficients
        out = []
        for r in rows:
            x = self._features(r)
            out.append(self.encoder.lookup(_allele(r))
                       + sum(bi * xi for bi, xi in zip(b, x)))
        return out


# ---------------------------------------------------------------------------
# The sentinel.
# ---------------------------------------------------------------------------

class ShuffledSentinel(Baseline):
    """Any baseline, refitted on permuted labels. The null every column is read against.

        order="before"  default, the only reportable one. Every label-derived step,
            target encoding included, sees permuted labels. Measured pooled 0.4448,
            within-allele 0.4751.
        order="after"   the defect, reproducible on demand. The target encoding sees the
            true labels and only the weight fit sees the permutation. Measured pooled
            0.6083, within-allele 0.5120 -- looks like a model, is not a null.

    One draw is not a result: report `sentinel_distribution` over several seeds.
    """

    name = "ShuffledSentinel"

    def __init__(self, factory, *, eval_ids, seed: int = 0, order: str = "before",
                 **inner_kwargs):
        super().__init__(eval_ids=eval_ids)
        if order not in ("before", "after"):
            raise BaselineError(f"order must be 'before' or 'after', got {order!r}")
        self.order = order
        self.seed = int(seed)
        self.factory = factory
        self.inner = factory(eval_ids=eval_ids, **inner_kwargs)
        self.name = f"ShuffledSentinel[{self.inner.name}, order={order}, seed={seed}]"
        if not self.inner.uses_labels:
            raise BaselineError(
                f"{self.inner.name} reads no labels, so permuting them cannot change "
                "its scores. A sentinel over a label-free baseline is not a null; it is "
                "the baseline with a misleading label."
            )

    def permute(self, rows) -> list[dict]:
        """Return copies of `rows` with the labels permuted. Inputs are never mutated."""
        rows = list(rows)
        labels = [_label(r) for r in rows]
        rng = random.Random(self.seed)
        rng.shuffle(labels)
        out = []
        for r, y in zip(rows, labels):
            c = dict(r)
            c["label"] = y
            out.append(c)
        return out

    def fit(self, rows):
        rows = self.check_fitting_pool(rows)
        shuffled = self.permute(rows)
        if self.order == "before":
            self.inner._fit_encoders(shuffled)
        else:
            # The defect, reproduced exactly: the encoding is built from the true labels.
            self.inner._fit_encoders(rows)
        self.inner._fit_weights(shuffled)
        self.inner.is_fitted = True
        self.is_fitted = True
        return self

    def score(self, rows) -> list[float]:
        if not self.is_fitted:
            raise BaselineError(f"{self.name}: score() before fit()")
        return self.inner.score(rows)


@dataclass
class SentinelDistribution:
    """A sentinel's figures over seeds. Never report one element of this."""

    order: str
    inner: str
    seeds: list[int]
    values: list[float] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.values)

    @property
    def mean(self) -> float:
        return sum(self.values) / len(self.values)

    @property
    def sd(self) -> float:
        if len(self.values) < 2:
            return float("nan")
        m = self.mean
        return math.sqrt(sum((v - m) ** 2 for v in self.values) / (len(self.values) - 1))

    def quantile(self, q: float) -> float:
        s = sorted(self.values)
        if len(s) == 1:
            return s[0]
        pos = q * (len(s) - 1)
        lo = int(math.floor(pos))
        hi = min(lo + 1, len(s) - 1)
        return s[lo] + (pos - lo) * (s[hi] - s[lo])

    def __str__(self) -> str:
        return (f"{self.inner} order={self.order}  n={self.n}  "
                f"mean {self.mean:.4f}  sd {self.sd:.4f}  "
                f"[{min(self.values):.4f}, {max(self.values):.4f}]  "
                f"p05 {self.quantile(0.05):.4f}  p95 {self.quantile(0.95):.4f}")


def sentinel_distribution(factory, fit_rows, eval_rows, evaluate, *, eval_ids,
                          seeds=range(20), order: str = "before",
                          **inner_kwargs) -> SentinelDistribution:
    """Refit `factory` on permuted labels once per seed and summarise the figures.

    `evaluate(scores) -> float` is the caller's, so this module never picks between the
    stratified and pooled estimators. Fewer than two seeds raises: one draw at 0.6083
    looks like a finding.
    """
    seeds = list(seeds)
    if len(seeds) < 2:
        raise BaselineError(
            f"a sentinel distribution needs at least 2 seeds, got {len(seeds)}. "
            "One permutation is a draw, not a null."
        )
    vals = []
    name = None
    for s in seeds:
        sent = ShuffledSentinel(factory, eval_ids=eval_ids, seed=s, order=order,
                                **inner_kwargs)
        name = sent.inner.name
        sent.fit(fit_rows)
        vals.append(float(evaluate(sent.score(eval_rows))))
    return SentinelDistribution(order=order, inner=name or "?", seeds=seeds, values=vals)


# ---------------------------------------------------------------------------
# Residualisation.
# ---------------------------------------------------------------------------

def within_allele_rank(scores, rows) -> list[float]:
    """Mid-rank within allele, mapped to (0, 1) as `(rank - 0.5) / n`.

    A singleton allele lands at exactly 0.5: right for a rank transform, wrong for an
    AUROC, which is why `stratified_auroc` excludes singleton strata.
    """
    scores, rows = list(scores), list(rows)
    if len(scores) != len(rows):
        raise BaselineError(f"{len(scores)} scores against {len(rows)} rows")
    by: dict[str, list[int]] = defaultdict(list)
    for i, r in enumerate(rows):
        by[_allele(r)].append(i)
    out = [0.0] * len(scores)
    for ix in by.values():
        if len(ix) == 1:
            out[ix[0]] = 0.5
            continue
        r = midranks([scores[i] for i in ix])
        for k, i in enumerate(ix):
            out[i] = (r[k] - 0.5) / len(ix)
    return out


_RESIDUAL_BASIS = {
    "kd": lambda p: mean_kd(p),
    "kd2": lambda p: mean_kd(p) ** 2,
    "length": lambda p: float(len(p)),
}


def residualise(scores, rows, *, predictors=("kd", "kd2", "length")) -> list[float]:
    """Within-allele-rank the score, then OLS out [mean KD, KD^2, length].

    Measured costs on the infectious benchmark: PRIME-2.0 0.5348 -> 0.5183,
    TITANiAN-IM (leaked) 0.7040 -> 0.7113. An artifact shrinks, a leak does not move,
    and what is left is signal.

    CAVEAT, which must travel with any residualised number: residualising mean KD
    against itself gives 0.4072 on infectious (n=2,093) but 0.5748 on neoantigen
    (n=833), where a pure function of the basis must come out at chance. No residualised
    neoantigen figure from this function is trustworthy; report the self-check beside
    any number you publish.

    A rank-deficient design raises `SingularDesignError` rather than silently dropping a
    predictor.
    """
    scores, rows = list(scores), list(rows)
    bad = [p for p in predictors if p not in _RESIDUAL_BASIS]
    if bad:
        raise BaselineError(f"unknown residual predictors {bad}; "
                            f"available: {sorted(_RESIDUAL_BASIS)}")
    if not predictors:
        raise BaselineError("residualise with no predictors removes nothing")
    ranked = within_allele_rank(scores, rows)
    peps = [_peptide(r) for r in rows]
    X = [[1.0] + [_RESIDUAL_BASIS[p](pep) for p in predictors] for pep in peps]
    beta = ols(X, ranked)
    return [ranked[i] - sum(b * x for b, x in zip(beta, X[i])) for i in range(len(rows))]
