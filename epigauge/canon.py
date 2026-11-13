"""Canonicalisation. One function per key type, and nothing else may normalise a key.

Unparseable input raises; there is no coercion and no fallback value. Two spellings of
one allele made two copies of the same 833-row table join at
0 of 833 matches, and a missed allele join scored exactly 0.5000. Organism folding is
one trailing parenthetical and nothing else -- never substring: "human" occurs in 16
organism strings here and 13 of them are viruses.
"""

from __future__ import annotations

import re

__all__ = [
    "CanonError",
    "canonicalise_allele",
    "canonicalise_organism",
    "canonicalise_peptide",
    "AA",
]

AA = "ACDEFGHIKLMNPQRSTVWY"
_AA_SET = frozenset(AA)


class CanonError(ValueError):
    """A key could not be canonicalised. Never caught and defaulted inside this package."""


# HLA-A*02:01, A*02:01, A0201, HLA-A02:01, DRB1*01:01, H-2-Kb (mouse, passed through).
# Group and protein fields are 2-3 digits each; a trailing expression suffix (N/L/S/C/A/Q)
# is dropped because it does not change the binding groove.
_HLA = re.compile(
    r"^(?:HLA[-_])?"
    r"(?P<locus>[A-Z]{1,4}\d?)"
    r"\*?"
    r"(?P<group>\d{2,3})"
    r":?"
    r"(?P<protein>\d{2,3})"
    r"(?::\d{2,3})*"
    r"[NLSCAQ]?$"
)
_MOUSE = re.compile(r"^(?:H2|H-2)[-_]?(?P<allele>[A-Za-z0-9]+)$")


def canonicalise_allele(raw: str, *, fields: int = 2) -> str:
    """Two-field allele key: `'A0101'`, `'hla-a0101'`, `'HLA-A*01:01'` -> `'HLA-A*01:01'`.

    Also `'DRB1*01:01'` -> `'HLA-DRB1*01:01'`, `'H-2-Kb'` -> `'H2-Kb'`. Raises CanonError
    on anything unparseable; a caller tolerating bad rows must catch it and record the drop.
    """
    if fields != 2:
        raise CanonError(f"only two-field normalisation is implemented, got fields={fields}")
    if not isinstance(raw, str):
        raise CanonError(f"allele must be str, got {type(raw).__name__}")
    s = raw.strip().upper().replace(" ", "")
    if not s:
        raise CanonError("empty allele string")

    m = _MOUSE.match(s)
    if m:
        a = m.group("allele")
        return f"H2-{a[0].upper()}{a[1:].lower()}"

    m = _HLA.match(s)
    if not m:
        raise CanonError(f"unparseable allele {raw!r}")
    return f"HLA-{m.group('locus')}*{m.group('group')}:{m.group('protein')}"


def canonicalise_organism(raw: str) -> str:
    """Fold exactly one trailing parenthetical and lowercase. Nothing else.

    `'Homo sapiens (human)'` -> `'homo sapiens'`; `'Influenza A virus (A/PR/8/34(H1N1))'`
    -> `'influenza a virus'`; `'Human betaherpesvirus 5'` stays distinct, which substring
    matching on "human" would not.
    """
    if not isinstance(raw, str):
        raise CanonError(f"organism must be str, got {type(raw).__name__}")
    s = raw.strip()
    if not s:
        raise CanonError("empty organism string")
    # One trailing parenthetical, balanced, possibly containing nested parens.
    if s.endswith(")"):
        depth = 0
        for i in range(len(s) - 1, -1, -1):
            if s[i] == ")":
                depth += 1
            elif s[i] == "(":
                depth -= 1
                if depth == 0:
                    s = s[:i].strip()
                    break
    out = " ".join(s.lower().split())
    if not out:
        raise CanonError(f"organism {raw!r} reduced to empty after folding")
    return out


def canonicalise_peptide(raw: str, *, min_len: int = 8, max_len: int = 15) -> str:
    """Uppercase, strip whitespace, assert canonical residues and a length window.

    Non-canonical residues raise rather than being dropped or masked: the metrics here are
    defined over the 20 canonical letters, so an X or a modification must be excluded
    deliberately by the caller.
    """
    if not isinstance(raw, str):
        raise CanonError(f"peptide must be str, got {type(raw).__name__}")
    s = "".join(raw.split()).upper()
    if not s:
        raise CanonError("empty peptide string")
    bad = sorted(set(s) - _AA_SET)
    if bad:
        raise CanonError(f"non-canonical residue(s) {bad} in peptide {raw!r}")
    if not (min_len <= len(s) <= max_len):
        raise CanonError(f"peptide {s!r} length {len(s)} outside [{min_len}, {max_len}]")
    return s
