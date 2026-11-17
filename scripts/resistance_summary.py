#!/usr/bin/env python3
"""
TyphiWatch - resistance and typing summary.

Reads AMRFinderPlus output (results/resistance/*.amrfinder.tsv), typing results
(MLST, SeqSero2, Mykrobe/GenoTyphi) and sample metadata, and writes:

  results/resistance/gene_matrix.tsv        sample x gene presence/absence
  results/resistance/resistance_summary.tsv per-sample AMR calls and MDR/XDR status
  results/resistance/sample_summary.tsv     AMR + typing + metadata, one row per sample

IMPORTANT: these are genotype-based predictions from known markers. They are not
laboratory susceptibility results. See docs/interpretation.md.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import pandas as pd


# ----------------------------------------------------------------------------
# AMRFinderPlus parsing
# ----------------------------------------------------------------------------
def first_col(columns, candidates):
    for c in candidates:
        if c in columns:
            return c
    return None


def read_amrfinder(path):
    """Return a list of dicts: symbol, cls, subclass for AMR-type elements."""
    path = Path(path)
    if path.stat().st_size == 0:
        return []
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        return []
    cols = list(df.columns)
    c_sym = first_col(cols, ["Element symbol", "Gene symbol"])
    c_cls = first_col(cols, ["Class"])
    c_sub = first_col(cols, ["Subclass"])
    c_typ = first_col(cols, ["Type", "Element type"])
    if c_sym is None:
        raise ValueError(f"{path}: cannot find gene symbol column in {cols}")

    hits = []
    for _, r in df.iterrows():
        if c_typ and r[c_typ].upper() != "AMR":
            continue  # drop STRESS / VIRULENCE / etc. reported by --plus
        hits.append({
            "symbol": r[c_sym],
            "cls": r[c_cls].upper() if c_cls else "",
            "subclass": r[c_sub].upper() if c_sub else "",
        })
    return hits


def classify(hits):
    """Apply simple, transparent marker rules (documented in docs/interpretation.md)."""
    sym = [h["symbol"] for h in hits]
    cls = {h["cls"] for h in hits}
    sub = {h["subclass"] for h in hits}

    chloramphenicol = any(s.startswith("cat") or s.startswith("cml") or s.startswith("floR") for s in sym) \
        or any("PHENICOL" in c for c in cls)
    ampicillin = any(s.startswith("bla") for s in sym)
    trimethoprim = any(s.startswith("dfr") for s in sym) or any("TRIMETHOPRIM" in c for c in cls)
    sulfonamide = any(s.startswith("sul") for s in sym) or any("SULFONAMIDE" in c for c in cls)
    sxt = trimethoprim and sulfonamide

    qrdr = sorted(s for s in sym if re.match(r"^(gyrA|gyrB|parC|parE)_", s))
    qnr = [s for s in sym if s.startswith("qnr")]
    fq_ns = bool(qrdr or qnr)

    ceph3 = any(s.startswith(("blaCTX-M", "blaCMY", "blaDHA")) for s in sym) \
        or any("CEPHALOSPORIN" in x for x in sub)
    azithro = any(s.startswith("acrB_R717") or s in ("mph(A)", "ermB") for s in sym)

    mdr = chloramphenicol and ampicillin and sxt
    xdr = mdr and fq_ns and ceph3

    if xdr:
        profile = "XDR"
    elif mdr:
        profile = "MDR"
    elif fq_ns:
        profile = "FQ-NS"
    elif sym:
        profile = "Other AMR markers"
    else:
        profile = "No known AMR markers"

    return {
        "n_amr_markers": len(sym),
        "amr_markers": ";".join(sorted(set(sym))),
        "amr_classes": ";".join(sorted(c for c in cls if c)),
        "chloramphenicol": chloramphenicol,
        "ampicillin": ampicillin,
        "trimethoprim_sulfamethoxazole": sxt,
        "fluoroquinolone_ns": fq_ns,
        "qrdr_mutations": ";".join(qrdr),
        "qrdr_mutation_count": len(qrdr),
        "ceftriaxone_3gc": ceph3,
        "azithromycin": azithro,
        "mdr": mdr,
        "xdr": xdr,
        "amr_profile": profile,
    }


# ----------------------------------------------------------------------------
# Typing parsers
# ----------------------------------------------------------------------------
def parse_mlst(path):
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                return {"mlst_scheme": parts[1], "st": parts[2]}
    return {"mlst_scheme": "NA", "st": "NA"}


def parse_seqsero2(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return {"serotype": "NA"}
    row = rows[0]
    key = next((k for k in row if "predicted serotype" in k.lower()), None)
    return {"serotype": row.get(key, "NA") if key else "NA"}


def find_values(obj, key):
    """Recursively yield values stored under `key` in nested JSON."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key:
                yield v
            yield from find_values(v, key)
    elif isinstance(obj, list):
        for item in obj:
            yield from find_values(item, key)


def parse_mykrobe(path):
    """Best-effort GenoTyphi genotype extraction from Mykrobe typhi-panel JSON.

    The JSON layout differs between Mykrobe versions; check one real output file
    and adjust the candidate keys below if the genotype column comes out as NA.
    """
    try:
        data = json.loads(Path(path).read_text())
    except (json.JSONDecodeError, OSError):
        return {"genotype": "NA"}
    for key in ("genotype", "lineage"):
        for value in find_values(data, key):
            if isinstance(value, list) and value:
                value = value[0]
            if isinstance(value, str) and value:
                return {"genotype": value}
    return {"genotype": "NA"}


def collect_typing(typing_dir):
    typing_dir = Path(typing_dir)
    out = {}

    def add(sample, d):
        out.setdefault(sample, {}).update(d)

    for p in sorted((typing_dir / "mlst").glob("*.mlst.tsv")):
        add(p.name.removesuffix(".mlst.tsv"), parse_mlst(p))
    for p in sorted((typing_dir / "seqsero2").glob("*.seqsero2.tsv")):
        add(p.name.removesuffix(".seqsero2.tsv"), parse_seqsero2(p))
    for p in sorted((typing_dir / "mykrobe").glob("*.mykrobe.json")):
        add(p.name.removesuffix(".mykrobe.json"), parse_mykrobe(p))
    return out


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amr-dir", default="results/resistance")
    ap.add_argument("--typing-dir", default="results/typing")
    ap.add_argument("--metadata", default="data/metadata/metadata.tsv")
    ap.add_argument("--outdir", default="results/resistance")
    args = ap.parse_args()

    amr_files = sorted(Path(args.amr_dir).glob("*.amrfinder.tsv"))
    if not amr_files:
        raise SystemExit(f"No *.amrfinder.tsv files found in {args.amr_dir}")

    records, presence = [], {}
    for f in amr_files:
        sample = f.name.removesuffix(".amrfinder.tsv")
        hits = read_amrfinder(f)
        row = {"sample": sample}
        row.update(classify(hits))
        records.append(row)
        presence[sample] = {h["symbol"]: 1 for h in hits}

    res = pd.DataFrame(records)
    matrix = pd.DataFrame(presence).T.fillna(0).astype(int)
    matrix.index.name = "sample"
    matrix = matrix.reindex(res["sample"]).fillna(0).astype(int)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(outdir / "gene_matrix.tsv", sep="\t")
    res.to_csv(outdir / "resistance_summary.tsv", sep="\t", index=False)

    summary = res.copy()
    typing = collect_typing(args.typing_dir)
    if typing:
        summary = summary.merge(
            pd.DataFrame.from_dict(typing, orient="index").rename_axis("sample").reset_index(),
            on="sample", how="left",
        )
    meta_path = Path(args.metadata)
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t", dtype=str)
        keep = [c for c in ["sample", "country", "location", "year", "collection_date", "study_accession"] if c in meta.columns]
        summary = summary.merge(meta[keep], on="sample", how="left")
    summary.to_csv(outdir / "sample_summary.tsv", sep="\t", index=False)

    print(f"Samples: {len(res)}")
    print(res["amr_profile"].value_counts().to_string())
    print(f"Wrote {outdir / 'resistance_summary.tsv'}, {outdir / 'gene_matrix.tsv'}, {outdir / 'sample_summary.tsv'}")


if __name__ == "__main__":
    main()
