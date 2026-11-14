#!/usr/bin/env python3
"""Per-entry intervals and all-pairs contrasts for the arms whose raw scores ship.

The referee objection this answers: the leaderboards carry rank, within-allele and
pooled AUROC and no uncertainty at all, in a paper whose thesis is that this field ranks
methods on margins it cannot resolve. And the resolution limit itself came from four
hand-picked contrasts at B=2,000 in a script that calls itself a feasibility read.

WHAT THIS CAN AND CANNOT COVER. Only the score columns the package ships can be
recomputed, so only those get an interval here. The remaining leaderboard entries are
reference figures transcribed from their sources, whose per-row scores are either
withheld under licence or were never published; no interval is computable for them by
anyone, including the authors.

THE REFERENCE PROTEOMES ARE REQUIRED AND ARE NOT REDISTRIBUTED. The clustering basis is
source protein, assigned by exact peptide-substring match against nine UniProt reviewed
proteomes. The package does not redistribute them: they are 69 MB of third-party FASTA
with their own terms. Fetch the nine files listed in PROTEOME_INDEX below into
`supplementary/proteomes/` and the default basis works. Absent even one of them this
script EXITS NON-ZERO and names the missing files, their UniProt queries and their
checksums. It does not fall back.

TWO CLUSTERING BASES, AND THE BASIS STRING ALWAYS SAYS WHICH.

`--clusters proteome` (the DEFAULT) is the basis the package specifies:
`epigauge.bootstrap.assign_clusters` assigns each row to its source protein by exact
peptide-substring match against the nine reference proteomes, merging multi-protein hits
by union-find. Its basis string is whatever `assign_clusters`/`ClusterMap` produce
naturally, i.e. `cluster/...`.

`--clusters overlap` is the WEAKER substitute that was used when the proteomes were
believed to be unavailable. It unions rows by shared k-mer *with each other*: peptides
from one source protein that overlap land in one cluster, but two non-overlapping
peptides from the same protein do not. It therefore captures part of the dependence the
proteome map captures and none that it misses, so its intervals are conservative relative
to row resampling and anti-conservative relative to the proteome map. Its basis string
says `peptide-overlap`, never `source-protein`, and its artifact is written to a
`_overlap` filename so it can never overwrite the shipped proteome-basis record.

ONE THING THE BASIS STRING DOES NOT SAY. `cluster_basis` in the artifact, and the
`basis`/`notes` on the ClusterMap itself, read `peptide-overlap` under `--clusters
overlap` and `cluster` under `--clusters proteome`. The PER-CELL strings on
`Estimate.basis` and `Delta.basis` do not: `paired` hardcodes the unit token to `cluster`
and `interval` takes it from the DrawMatrix, whose unit must be `cluster` for `paired` to
accept it at all. So a per-cell basis reads `cluster/B/provenance/digest/assigned` on
BOTH bases, and the two are told apart only by the trailing assigned fraction (0.336 vs
0.859 on the infectious arm) and the seed digest. A table that prints per-cell bases must
therefore also carry `cluster_basis`, or those two cells look interchangeable.

    python3 -I supplementary/scripts/leaderboard_uncertainty.py \
        [--set titanian_infectious] [--B 10000] [--proteomes DIR] \
        [--clusters proteome|overlap]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
from pathlib import Path

# The released package root: `supplementary/scripts/x.py` -> two levels up. Inserted
# explicitly rather than relying on PYTHONPATH, so `python3 -I` works from any cwd.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from epigauge import bootstrap as bs
    from epigauge import data, metrics
except ImportError as exc:  # pragma: no cover - environment, not a data finding
    sys.stderr.write(
        f"cannot import epigauge from {ROOT}: {exc}\n"
        "Run this from an unpacked epigauge source distribution, or install the "
        "package first.\n"
    )
    raise SystemExit(2)

SCRIPT_REL = "supplementary/scripts/leaderboard_uncertainty.py"
OUT_DIR = ROOT / "supplementary" / "artifacts"

# ------------------------------------------------------------------ the nine proteomes
#
# UniProt reviewed (Swiss-Prot) proteomes, downloaded 2026-11-10. The query is the exact
# UniProt REST search behind each file; `sequences` and `sha256` identify the version,
# because the clustering and therefore every interval below is specific to it.
PROTEOME_INDEX = {
    "escherichia_reviewed.fasta": {
        "query": "(taxonomy_id:561) AND (reviewed:true)", "sequences": 23996,
        "sha256": "4a474daee45e65ed1d221a8a0f9f320c31aaa723266ce049dbd9d1104e28fab7"},
    "human_reviewed.fasta": {
        "query": "(organism_id:9606) AND (reviewed:true)", "sequences": 20429,
        "sha256": "cee55a0d0af00b26876126c613abdff9d2e8649eb11a73e8db6c1badcd3d8fd5"},
    "mouse_reviewed.fasta": {
        "query": "(organism_id:10090) AND (reviewed:true)", "sequences": 17285,
        "sha256": "1d590b004e502ec8e3e1edca352afef1e71e0fded5eb2519c7a29c46feda613f"},
    "mycobacterium_reviewed.fasta": {
        "query": "(taxonomy_id:1763) AND (reviewed:true)", "sequences": 10141,
        "sha256": "10b9d8fe3b228bc0849b1d732c9d859b66c7e93debb550cdc7a5831cbd266773"},
    "plasmodium_reviewed.fasta": {
        "query": "(taxonomy_id:5820) AND (reviewed:true)", "sequences": 774,
        "sha256": "4448ead0918a9482e518c7f9f76598f20c7d9efb5931a9b5fbff2e0559992afd"},
    "salmonella_reviewed.fasta": {
        "query": "(taxonomy_id:590) AND (reviewed:true)", "sequences": 12890,
        "sha256": "1dad6b2e397f67bba0189e5faa4331bba2e05e23a9cfea7d40c89827ee1ed89b"},
    "staphylococcus_reviewed.fasta": {
        "query": "(taxonomy_id:1279) AND (reviewed:true)", "sequences": 13512,
        "sha256": "039658207b5507c535de4f7e5dbe1d7cdbdbc22cc21272b81d247042bbeed7d8"},
    "streptococcus_reviewed.fasta": {
        "query": "(taxonomy_id:1301) AND (reviewed:true)", "sequences": 12805,
        "sha256": "cea6d6e76555bd1a76a4623c8d509c7a2859c3eaf7c8cadd8b812ae410b13508"},
    "virus_reviewed.fasta": {
        "query": "(taxonomy_id:10239) AND (reviewed:true)", "sequences": 17541,
        "sha256": "00cfa76a9039d541ea4f1021706e66f61e21d30f77a1f215cfe0110f5fdc6239"},
}


def proteome_dir(cli_value: str | None) -> Path:
    """`--proteomes`, then `$EPIGAUGE_PROTEOMES`, then `supplementary/proteomes/`."""
    if cli_value:
        return Path(cli_value).expanduser().resolve()
    env = os.environ.get("EPIGAUGE_PROTEOMES")
    if env:
        return Path(env).expanduser().resolve()
    return ROOT / "supplementary" / "proteomes"


class MissingProteomes(FileNotFoundError):
    """Raised when the source-protein basis cannot be built. Never caught here."""


def proteome_paths(directory: Path) -> list[Path]:
    """The nine reference proteomes, or a hard failure naming every one that is absent.

    A PARTIAL set is the dangerous case and is refused like an empty one. Assignment
    against the human proteome alone is 4.6% on the infectious arm against 85.9% for all
    nine, so a missing file does not make the basis slightly worse: it makes almost every
    row a singleton, which is the row-level design wearing a clustered name.
    """
    missing = [n for n in PROTEOME_INDEX if not (directory / n).is_file()]
    if missing:
        lines = [
            f"cannot build the source-protein clustering basis: {len(missing)} of "
            f"{len(PROTEOME_INDEX)} reference proteomes are missing from {directory}",
            "",
            "The package does not redistribute them. They are UniProt reviewed "
            "(Swiss-Prot) proteomes, 69 MB in total, fetched from the UniProt REST API:",
            "",
            "    https://rest.uniprot.org/uniprotkb/stream"
            "?format=fasta&query=<QUERY>",
            "",
            "Missing files, with the exact query and the checksum of the version these "
            "figures were computed on:",
        ]
        for n in missing:
            e = PROTEOME_INDEX[n]
            lines += [
                f"  {n}",
                f"      query     {e['query']}",
                f"      sequences {e['sequences']}",
                f"      sha256    {e['sha256']}",
            ]
        lines += [
            "",
            f"Put all {len(PROTEOME_INDEX)} files in {directory}, or point "
            "--proteomes / $EPIGAUGE_PROTEOMES at the directory holding them.",
            "",
            "`--clusters overlap` is NOT a fallback. It is a different and weaker basis "
            "whose intervals are anti-conservative relative to this one; it is offered "
            "only so the substitute used before the proteomes were available can be "
            "reproduced, and it writes to a separate `_overlap` filename.",
        ]
        raise MissingProteomes("\n".join(lines))
    return [directory / n for n in sorted(PROTEOME_INDEX)]


def proteome_clusters(peptides: list[str], directory: Path) -> bs.ClusterMap:
    """Source-protein clusters against all nine reference proteomes.

    Nothing is relabelled: `unit` stays `cluster` and `notes` stay whatever
    `assign_clusters` recorded, so the basis string is the package's own.
    """
    paths = proteome_paths(directory)
    cm = bs.assign_clusters(peptides, paths)
    cm.notes["n_proteomes"] = len(paths)
    return cm


def basis_label(cm: bs.ClusterMap) -> str:
    """The summary basis label: `peptide-overlap` for the substitute, `cluster` otherwise."""
    return cm.notes.get("basis") or cm.unit


def overlap_clusters(peptides: list[str]) -> bs.ClusterMap:
    """Union-find over peptides sharing a k-mer window (9, or 8 for 8-mers).

    Mirrors the window rule in `assign_clusters` so the two maps differ only in what
    they union against.
    """
    parent = list(range(len(peptides)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    index: dict[str, int] = {}
    for i, pep in enumerate(peptides):
        k = 8 if len(pep) == 8 else 9
        for s in range(0, max(1, len(pep) - k + 1)):
            kmer = pep[s:s + k]
            if len(kmer) < k:
                continue
            if kmer in index:
                union(i, index[kmer])
            else:
                index[kmer] = i

    keys = [f"ov:{find(i)}" for i in range(len(peptides))]
    n_non_singleton_rows = sum(
        1 for i in range(len(peptides))
        if sum(1 for j in range(len(peptides)) if find(j) == find(i)) > 1
    ) if len(peptides) <= 4000 else None

    cm = bs.ClusterMap.from_keys(
        keys,
        n_assigned=n_non_singleton_rows if n_non_singleton_rows is not None else len(keys),
    )
    cm.unit = "peptide-overlap"   # the basis string must not read "cluster"
    cm.notes["basis"] = "peptide-overlap"
    cm.notes["caveat"] = (
        "unions rows by shared k-mer with each other, not against reference proteomes; "
        "two non-overlapping peptides from one protein are not unioned"
    )
    return cm


def columns_of(set_id: str) -> list[str]:
    rows = data.load(set_id, scores=True)
    return sorted({c for r in rows for c in r.get("scores", {})})


def wa(labels, scores, groups) -> float:
    return metrics.stratified_auroc(labels, scores, groups).value


def artifact_name(set_id: str, mode: str) -> str:
    """`leaderboard_uncertainty_infectious.json`, as the released layout spells it.

    The `titanian_` prefix is dropped because the shipped artifact names carry the arm,
    not the set id. An overlap run gets its own `_overlap` stem so it cannot overwrite
    the proteome-basis record that ships.
    """
    arm = set_id[len("titanian_"):] if set_id.startswith("titanian_") else set_id
    suffix = "" if mode == "proteome" else "_overlap"
    return f"leaderboard_uncertainty_{arm}{suffix}.json"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", dest="set_id", default="titanian_infectious")
    ap.add_argument("--B", type=int, default=10_000)
    ap.add_argument("--clusters", choices=("proteome", "overlap"), default="proteome",
                    help="resampling basis: source-protein against the nine proteomes "
                         "(default), or the weaker peptide-overlap substitute")
    ap.add_argument("--proteomes", default=None,
                    help="directory holding the nine reference FASTA files "
                         "(default: supplementary/proteomes/, or $EPIGAUGE_PROTEOMES)")
    args = ap.parse_args(argv)

    rows = data.load(args.set_id, scores=True)
    cols = columns_of(args.set_id)
    if not cols:
        print(f"could not run: {args.set_id} ships no score columns", file=sys.stderr)
        return 2

    pdir = proteome_dir(args.proteomes)
    peptides = [r["peptide"] for r in rows]
    try:
        clusters = (proteome_clusters(peptides, pdir) if args.clusters == "proteome"
                    else overlap_clusters(peptides))
    except MissingProteomes as exc:
        print(exc, file=sys.stderr)
        return 3

    print(f"{args.set_id}: {len(rows)} rows, {len(cols)} shipped score columns")
    print(f"clusters: {clusters}")
    print(f"  basis: {basis_label(clusters)}")
    if args.clusters == "proteome":
        print(f"  proteomes: {len(clusters.notes['proteomes'])} "
              f"({', '.join(clusters.notes['proteomes'])})")
        print(f"  proteome directory: {pdir}")
        print(f"  proteins hit before merge: {clusters.notes['n_proteins_hit_before_merge']}")
    print()

    recipe = bs.SeedRecipe(version="1", set_id=args.set_id, comparison_id="leaderboard-ci")
    dm = bs.draws(clusters.n_clusters, B=args.B, seed_recipe=recipe)

    labels = [int(r["label"]) for r in rows]
    groups = [r["allele"] for r in rows]

    # --- per-entry intervals
    entries = {}
    print(f"{'entry':44s} {'within-allele':>14s}  95% CI")
    for c in cols:
        sc = [r.get("scores", {}).get(c) for r in rows]
        if any(v is None for v in sc):
            print(f"{c[:44]:44s}   incomplete column, skipped")
            continue
        est = bs.interval(wa, labels, [float(v) for v in sc], groups, clusters, dm)
        entries[c] = {
            "value": round(est.value, 4), "lo": round(est.lo, 4), "hi": round(est.hi, 4),
            "mc_se": [round(v, 5) for v in est.mc_se], "basis": est.basis,
            "n_clusters": est.n_clusters, "assigned_fraction": round(est.assigned_fraction, 4),
            "replicates_failed": est.n_replicates_failed, "provisional": est.provisional,
        }
        print(f"{c[:44]:44s} {est.value:14.4f}  [{est.lo:.4f}, {est.hi:.4f}]"
              f"{'  [PROVISIONAL]' if est.provisional else ''}")

    # --- all pairs
    names = [c for c in cols if c in entries]
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    print(f"\nall {len(pairs)} pairwise contrasts over the {len(names)} recomputable entries")

    contrasts = []
    for a, b in pairs:
        # Orient each pair so the higher point estimate is `a`. `paired` reports a
        # one-sided p for "a beats b", so an unoriented pair returns p ~ 1 by
        # construction and any family-wise correction over a mix of orientations is
        # meaningless.
        if entries[b]["value"] > entries[a]["value"]:
            a, b = b, a
        sa = [float(r["scores"][a]) for r in rows]
        sb = [float(r["scores"][b]) for r in rows]
        rec = bs.SeedRecipe(version="1", set_id=args.set_id, comparison_id=f"{a}|{b}")
        d = bs.paired(wa, labels, sa, sb, groups, clusters,
                      bs.draws(clusters.n_clusters, B=args.B, seed_recipe=rec))
        m = bs.mde(d)
        contrasts.append({
            "a": a, "b": b, "delta": round(d.d, 4), "lo": round(d.lo, 4),
            "hi": round(d.hi, 4), "p": round(d.p, 4), "se": round(d.se, 4),
            "design_effect": round(d.design_effect, 3) if d.design_effect else None,
            "n_eff": round(d.n_eff, 1) if d.n_eff else None,
            "mde": round(m.mde, 4), "resolvable": bool(abs(d.d) >= m.mde),
            "basis": d.basis,
        })

    # Holm-Bonferroni over the family, because 36 contrasts on one arm is a family.
    ordered = sorted(range(len(contrasts)), key=lambda i: contrasts[i]["p"])
    n = len(ordered)
    running = 0.0
    for rank, i in enumerate(ordered):
        adj = min(1.0, max(running, contrasts[i]["p"] * (n - rank)))
        running = adj
        contrasts[i]["p_holm"] = round(adj, 4)

    ses = sorted(c["se"] for c in contrasts)
    mdes = sorted(c["mde"] for c in contrasts)
    resolvable = [c for c in contrasts if c["resolvable"]]
    sig_holm = [c for c in contrasts if c["p_holm"] < 0.05]

    print(f"  paired SE      min {ses[0]:.4f}  median {ses[len(ses)//2]:.4f}  max {ses[-1]:.4f}")
    print(f"  MDE            min {mdes[0]:.4f}  median {mdes[len(mdes)//2]:.4f}  max {mdes[-1]:.4f}")
    print(f"  |delta| >= MDE                       {len(resolvable)} of {len(contrasts)}")
    print(f"  Holm-adjusted p < 0.05               {len(sig_holm)} of {len(contrasts)}")
    if resolvable:
        print("  resolvable contrasts:")
        for c in sorted(resolvable, key=lambda c: -abs(c["delta"])):
            print(f"    {c['a'][:30]:32s} vs {c['b'][:30]:32s} "
                  f"d={c['delta']:+.4f} MDE={c['mde']:.4f} p_holm={c['p_holm']:.4f}")

    rec = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": SCRIPT_REL,
        "set_id": args.set_id,
        "B": args.B,
        "cluster_mode": args.clusters,
        "cluster_basis": basis_label(clusters),
        "cluster_caveat": clusters.notes.get("caveat"),
        "cluster_notes": {k: v for k, v in clusters.notes.items() if k != "caveat"},
        "cluster_provisional": clusters.provisional,
        "n_clusters": clusters.n_clusters,
        "n_rows": clusters.n_rows,
        "assigned_fraction": round(clusters.assigned_fraction, 4),
        "n_recomputable_entries": len(entries),
        "entries": entries,
        "n_contrasts": len(contrasts),
        "se_min": ses[0], "se_median": ses[len(ses) // 2], "se_max": ses[-1],
        "mde_min": mdes[0], "mde_median": mdes[len(mdes) // 2], "mde_max": mdes[-1],
        "n_resolvable": len(resolvable),
        "n_significant_holm": len(sig_holm),
        "contrasts": contrasts,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / artifact_name(args.set_id, args.clusters)
    out.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"\nARTIFACT WRITTEN: {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
