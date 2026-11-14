#!/usr/bin/env python3
"""Score the four shipped arms that carry no predictor columns.

WHY THESE FOUR MATTER. Every competitor figure on the two scored arms is reached by the
contamination census: the predictor columns arrived inside the evaluation tables, and on
the neoantigen arm nine entries are scored for methods whose own deposit contains every
row. These four arms -- 441 distinct rows, zero overlap with `im_test` and absent from the
evaluation deposit -- are the only rows in the package where a published method's figure
would be tier 0. Without them the package has two arms where nearly every figure is
affected and none where any is not.

WHAT IS SCORED, AND UNDER WHICH LICENCE. MHCflurry's licence is Apache-2.0 and the project
licence review records output redistribution as unrestricted, which is why it is the
designated presentation baseline. The unfitted baselines are functions of the row's own
fields and carry no licence question at all.

VERSION PRECISION. The figures that ship were produced with **MHCflurry 2.3.28** and model
set **models_class1_presentation 2.3.0**. That is NOT "MHCflurry 2.0", and the two scored
arms' shipped `mhcflurry_2_0_*` columns came from a different release. A figure produced
here must never be placed in the same column as those. The version actually installed is
recorded in the artifact, so a re-run under a different release is visible rather than
silently substituted.

TWO EXTERNAL DEPENDENCIES, NEITHER REDISTRIBUTED.

  1. The nine UniProt reviewed reference proteomes, listed in PROTEOME_INDEX below. They
     carry the default clustering basis. Absent even one of them this script EXITS
     NON-ZERO and names the missing files, their queries and their checksums. It does not
     fall back to a weaker basis.
  2. MHCflurry, for its column only. `pip install mhcflurry && mhcflurry-downloads fetch
     models_class1_presentation`. Without it, pass `--no-mhcflurry` and the two unfitted
     baselines still run; the MHCflurry entry is then absent rather than imputed.

CLUSTERING BASIS. `--clusters proteome` (the DEFAULT) uses the basis the package
specifies: `epigauge.bootstrap.assign_clusters` against the nine reference proteomes,
whose basis string is the package's own `cluster/...`. `--clusters overlap` keeps the
weaker substitute used when those proteomes were believed unavailable; its basis string
reads `peptide-overlap` and its artifact goes to a separate `_overlap` filename, so the
two can never be confused in a table or overwrite one another on disk.

    python3 -I supplementary/scripts/score_unscored_arms.py \
        [--B 10000] [--no-mhcflurry] [--proteomes DIR] \
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

SCRIPT_REL = "supplementary/scripts/score_unscored_arms.py"
OUT_DIR = ROOT / "supplementary" / "artifacts"

ARMS = ("hbv_hlabc_screen", "xmuimmu_noncanonical", "independent_composite", "itsndb")

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


def load_release():
    try:
        from epigauge import baselines, bootstrap, data, metrics
    except ImportError as exc:  # pragma: no cover - environment, not a data finding
        sys.stderr.write(
            f"cannot import epigauge from {ROOT}: {exc}\n"
            "Run this from an unpacked epigauge source distribution, or install the "
            "package first.\n"
        )
        raise SystemExit(2)
    return baselines, bootstrap, data, metrics


def proteome_dir(cli_value):
    """`--proteomes`, then `$EPIGAUGE_PROTEOMES`, then `supplementary/proteomes/`."""
    if cli_value:
        return Path(cli_value).expanduser().resolve()
    env = os.environ.get("EPIGAUGE_PROTEOMES")
    if env:
        return Path(env).expanduser().resolve()
    return ROOT / "supplementary" / "proteomes"


class MissingProteomes(FileNotFoundError):
    """Raised when the source-protein basis cannot be built. Never caught here."""


def proteome_paths(directory):
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


def proteome_clusters(bootstrap, peptides, directory):
    """Source-protein clusters against all nine reference proteomes.

    Nothing is relabelled: `unit` stays `cluster` and `notes` stay whatever
    `assign_clusters` recorded, so the basis string is the package's own.
    """
    paths = proteome_paths(directory)
    cm = bootstrap.assign_clusters(peptides, paths)
    cm.notes["n_proteomes"] = len(paths)
    return cm


def basis_label(cm):
    """`peptide-overlap` for the substitute, `cluster` for the package's own basis."""
    return cm.notes.get("basis") or cm.unit


def overlap_clusters(bootstrap, peptides):
    """Union rows sharing a k-mer window. Same construction as leaderboard_uncertainty."""
    parent = list(range(len(peptides)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    index = {}
    for i, pep in enumerate(peptides):
        k = 8 if len(pep) == 8 else 9
        for s in range(0, max(1, len(pep) - k + 1)):
            kmer = pep[s:s + k]
            if len(kmer) < k:
                continue
            if kmer in index:
                ra, rb = find(i), find(index[kmer])
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
            else:
                index[kmer] = i

    keys = [f"ov:{find(i)}" for i in range(len(peptides))]
    sizes = {}
    for k in keys:
        sizes[k] = sizes.get(k, 0) + 1
    cm = bootstrap.ClusterMap.from_keys(
        keys, n_assigned=sum(1 for k in keys if sizes[k] > 1))
    cm.unit = "peptide-overlap"
    cm.notes["basis"] = "peptide-overlap"
    return cm


def mhcflurry_scores(rows):
    """Presentation score per row, plus the coverage the predictor actually supported.

    Returns (scores, n_supported, n_total, version, model_version). A row whose allele the
    predictor does not support gets None and is excluded from that arm's figure rather than
    imputed -- an unsupported allele is missing data, not a prediction of 0.5.
    """
    import mhcflurry
    from mhcflurry import Class1PresentationPredictor

    pred = Class1PresentationPredictor.load()
    supported = set(pred.supported_alleles)

    # `predict` treats an allele LIST as one individual's set and takes the strongest
    # prediction across it, capped at 6. So group by allele and call once per allele;
    # each call then has a single-allele set and the score is unambiguous.
    by_allele = {}
    for i, r in enumerate(rows):
        if r["allele"] in supported:
            by_allele.setdefault(r["allele"], []).append(i)

    out = [None] * len(rows)
    n = 0
    for allele, idxs in by_allele.items():
        peps = [rows[i]["peptide"] for i in idxs]
        df = pred.predict(peptides=peps, alleles=[allele], verbose=0)
        col = "presentation_score" if "presentation_score" in df.columns else df.columns[-1]
        vals = df[col].tolist()
        if len(vals) != len(idxs):
            raise RuntimeError(
                f"{allele}: predictor returned {len(vals)} scores for {len(idxs)} peptides")
        for i, v in zip(idxs, vals):
            out[i] = float(v)
        n += len(idxs)
    return out, n, len(rows), mhcflurry.__version__, "models_class1_presentation 2.3.0"


MHCFLURRY_HINT = (
    "install it with `pip install mhcflurry` followed by `mhcflurry-downloads fetch "
    "models_class1_presentation`, or pass --no-mhcflurry to run the unfitted baselines "
    "alone"
)


def evaluate(metrics, bootstrap, rows, scores, cm, dm, name):
    """within-allele + pooled on the rows this scorer actually covered."""
    keep = [i for i, v in enumerate(scores) if v is not None]
    if len(keep) < 20:
        return {"entry": name, "status": "too few scored rows", "n_scored": len(keep)}
    lab = [int(rows[i]["label"]) for i in keep]
    sc = [float(scores[i]) for i in keep]
    grp = [rows[i]["allele"] for i in keep]
    if len(set(lab)) < 2:
        return {"entry": name, "status": "single-class after coverage", "n_scored": len(keep)}

    rec = {"entry": name, "n_scored": len(keep), "n_total": len(rows),
           "coverage": round(len(keep) / len(rows), 4)}
    try:
        st = metrics.stratified_auroc(lab, sc, grp)
        rec["within_allele"] = round(st.value, 4)
        rec["strata"] = f"{st.n_strata_contributing}/{st.n_strata_total}"
        rec["rows_evaluable"] = f"{st.n_rows_evaluable}/{st.n_rows_total}"
        rec["provisional"] = st.provisional
    except Exception as e:
        rec["within_allele"] = None
        rec["within_allele_refused"] = f"{type(e).__name__}: {e}"
    rec["pooled"] = round(metrics.auroc(lab, sc), 4)

    # An interval only where the cluster map lines up with the kept rows.
    if len(keep) == len(rows):
        try:
            est = bootstrap.interval(
                lambda L, S, G: metrics.stratified_auroc(L, S, G).value,
                lab, sc, grp, cm, dm)
            rec["ci95"] = [round(est.lo, 4), round(est.hi, 4)]
            rec["basis"] = est.basis
            rec["interval_provisional"] = est.provisional
        except Exception as e:
            rec["ci95"] = None
            rec["interval_refused"] = f"{type(e).__name__}: {e}"
    else:
        rec["ci95"] = None
        rec["interval_refused"] = (
            "partial allele coverage; the cluster map is built over all rows and "
            "resampling a subset of it would not be the stated design")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--B", type=int, default=10_000)
    ap.add_argument("--no-mhcflurry", action="store_true")
    ap.add_argument("--clusters", choices=("proteome", "overlap"), default="proteome",
                    help="resampling basis: source-protein against the nine proteomes "
                         "(default), or the weaker peptide-overlap substitute")
    ap.add_argument("--proteomes", default=None,
                    help="directory holding the nine reference FASTA files "
                         "(default: supplementary/proteomes/, or $EPIGAUGE_PROTEOMES)")
    args = ap.parse_args(argv)

    baselines, bootstrap, data, metrics = load_release()
    pdir = proteome_dir(args.proteomes)

    # Fail before any arm is scored, so a missing proteome costs nothing and cannot be
    # mistaken for a partial result.
    if args.clusters == "proteome":
        try:
            proteome_paths(pdir)
        except MissingProteomes as exc:
            print(exc, file=sys.stderr)
            return 3

    results = {}
    for arm in ARMS:
        rows = data.load(arm)
        peptides = [r["peptide"] for r in rows]
        cm = (proteome_clusters(bootstrap, peptides, pdir) if args.clusters == "proteome"
              else overlap_clusters(bootstrap, peptides))
        dm = bootstrap.draws(
            cm.n_clusters, B=args.B,
            seed_recipe=bootstrap.SeedRecipe(version="1", set_id=arm,
                                             comparison_id="unscored-arm"))
        entries = []

        # Unfitted baselines: functions of the row's own fields, no fitting pool needed.
        for cls, label in ((baselines.Length, "BASELINE — peptide length"),
                           (baselines.Hydropathy, "BASELINE — mean hydropathy")):
            b = cls(eval_ids=())
            b.fit(rows)
            entries.append(evaluate(metrics, bootstrap, rows, b.score(rows), cm, dm, label))

        if not args.no_mhcflurry:
            try:
                sc, nsup, ntot, ver, mver = mhcflurry_scores(rows)
                e = evaluate(metrics, bootstrap, rows, sc, cm, dm,
                             f"MHCflurry {ver} (presentation score)")
                e["alleles_supported_rows"] = f"{nsup}/{ntot}"
                e["model_version"] = mver
                entries.append(e)
            except Exception as exc:
                entries.append({
                    "entry": "MHCflurry",
                    "status": f"failed: {type(exc).__name__}: {exc}",
                    "hint": MHCFLURRY_HINT,
                })

        # Paired contrasts between every pair of entries that produced a figure. On an
        # arm this small the question is not the ranking but whether anything separates.
        scored = [(e["entry"], e) for e in entries if e.get("within_allele") is not None]
        contrasts = []
        if len(scored) > 1:
            series = {}
            for cls, label in ((baselines.Length, "BASELINE — peptide length"),
                               (baselines.Hydropathy, "BASELINE — mean hydropathy")):
                b = cls(eval_ids=())
                b.fit(rows)
                series[label] = b.score(rows)
            if not args.no_mhcflurry:
                try:
                    sc, _n, _t, ver, _m = mhcflurry_scores(rows)
                    if all(v is not None for v in sc):
                        series[f"MHCflurry {ver} (presentation score)"] = sc
                except Exception:
                    pass
            lab = [int(r["label"]) for r in rows]
            grp = [r["allele"] for r in rows]
            names = [n for n, _ in scored if n in series]
            def wa_fn(labels, scores, groups):
                return metrics.stratified_auroc(labels, scores, groups).value
            for i, a in enumerate(names):
                for bname in names[i + 1:]:
                    hi, lo = (a, bname)
                    if dict(scored)[bname]["within_allele"] > dict(scored)[a]["within_allele"]:
                        hi, lo = bname, a
                    rec_seed = bootstrap.SeedRecipe(version="1", set_id=arm,
                                                    comparison_id=f"{hi}|{lo}")
                    try:
                        d = bootstrap.paired(
                            wa_fn, lab, [float(v) for v in series[hi]],
                            [float(v) for v in series[lo]], grp, cm,
                            bootstrap.draws(cm.n_clusters, B=args.B, seed_recipe=rec_seed))
                        m = bootstrap.mde(d)
                        contrasts.append({
                            "a": hi, "b": lo, "delta": round(d.d, 4),
                            "lo": round(d.lo, 4), "hi": round(d.hi, 4),
                            "p": round(d.p, 4), "se": round(d.se, 4),
                            "mde": round(m.mde, 4), "resolvable": bool(abs(d.d) >= m.mde),
                            "basis": d.basis,
                        })
                    except Exception as e:
                        contrasts.append({"a": hi, "b": lo,
                                          "refused": f"{type(e).__name__}: {e}"})

        results[arm] = {
            "n_rows": len(rows),
            "n_positive": sum(1 for r in rows if int(r["label"]) == 1),
            "n_alleles": len({r["allele"] for r in rows}),
            "clusters": cm.n_clusters,
            "assigned_fraction": round(cm.assigned_fraction, 4),
            "cluster_basis": basis_label(cm),
            "cluster_provisional": cm.provisional,
            "entries": entries,
            "contrasts": contrasts,
        }

        print(f"\n=== {arm}: {len(rows)} rows, {results[arm]['n_alleles']} alleles, "
              f"{cm.n_clusters} clusters (assigned {cm.assigned_fraction:.3f}, "
              f"basis {basis_label(cm)}"
              f"{', PROVISIONAL' if cm.provisional else ''}) ===")
        for e in entries:
            if e.get("within_allele") is None:
                print(f"  {e['entry'][:46]:48s} {e.get('status') or e.get('within_allele_refused','')}")
                continue
            ci = f"[{e['ci95'][0]:.4f}, {e['ci95'][1]:.4f}]" if e.get("ci95") else "(no interval)"
            print(f"  {e['entry'][:46]:48s} wa={e['within_allele']:.4f} "
                  f"pooled={e['pooled']:.4f} {ci} strata={e.get('strata','')} "
                  f"cov={e.get('coverage')}")
        for c in contrasts:
            if "refused" in c:
                print(f"    {c['a'][:26]:28s} vs {c['b'][:26]:28s} REFUSED {c['refused'][:40]}")
            else:
                print(f"    {c['a'][:26]:28s} vs {c['b'][:26]:28s} "
                      f"d={c['delta']:+.4f} MDE={c['mde']:.4f} p={c['p']:.4f} "
                      f"{'RESOLVABLE' if c['resolvable'] else '-'}")

    rec = {
        "produced_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "script": SCRIPT_REL,
        "B": args.B,
        "cluster_mode": args.clusters,
        "cluster_basis": ("peptide-overlap" if args.clusters == "overlap" else "cluster"),
        "proteomes": ([p.name for p in proteome_paths(pdir)]
                      if args.clusters == "proteome" else None),
        "note": ("these arms carry no shipped predictor columns; every figure here is "
                 "computed by this script and is tier 0 under the contamination census"),
        "arms": results,
    }
    suffix = "" if args.clusters == "proteome" else "_overlap"
    out = OUT_DIR / f"unscored_arms_scored{suffix}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"\nARTIFACT WRITTEN: {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
