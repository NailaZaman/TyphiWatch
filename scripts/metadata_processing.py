#!/usr/bin/env python3
"""
TyphiWatch - metadata and input preparation (standard library only).

Subcommands
-----------
search     Find candidate S. Typhi Illumina paired-end runs in ENA
           (e.g. by country) and write them as an accession list.
resolve    Turn data/accession_list.tsv (run / experiment / study accessions)
           into run-level metadata.tsv and the pipeline samplesheet.csv.
reference  Download the reference genome (default: S. Typhi CT18, AL513382.1)
           from NCBI as GenBank.

Examples
--------
python scripts/metadata_processing.py search --country Pakistan --max 60
python scripts/metadata_processing.py resolve --max-per-accession 10
python scripts/metadata_processing.py reference
"""
import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TYPHI_TAXID = "90370"  # Salmonella enterica subsp. enterica serovar Typhi

FIELDS = [
    "run_accession", "sample_accession", "study_accession", "scientific_name",
    "country", "collection_date", "instrument_model", "library_layout",
    "base_count", "fastq_ftp",
]

METADATA_COLUMNS = [
    "sample", "run_accession", "sample_accession", "study_accession",
    "scientific_name", "country", "location", "collection_date", "year",
    "instrument_model", "library_layout", "base_count", "fastq_1", "fastq_2",
]


def http_get(url, retries=3, pause=3.0, timeout=180):
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "TyphiWatch/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            time.sleep(pause * attempt)
    raise RuntimeError(f"Request failed after {retries} attempts:\n  {url}\n  {last}")


def read_accessions(path):
    accessions = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        first = re.split(r"[\t,]", line)[0].strip()
        if first.lower() in ("accession", "run_accession"):
            continue
        accessions.append(first)
    return list(dict.fromkeys(accessions))  # de-duplicate, keep order


def parse_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def fetch_runs(accession):
    """Return ENA read_run records for a run / experiment / sample / study accession."""
    params = urllib.parse.urlencode(
        {"accession": accession, "result": "read_run", "fields": ",".join(FIELDS), "format": "tsv"}
    )
    return parse_tsv(http_get(f"{ENA_API}/filereport?{params}"))


def clean_row(rec):
    country_raw = rec.get("country", "") or ""
    country, _, location = country_raw.partition(":")
    year_match = re.search(r"(18\d{2}|19\d{2}|20\d{2})", rec.get("collection_date", "") or "")
    urls = [u for u in (rec.get("fastq_ftp", "") or "").split(";") if u]
    r1 = next((u for u in urls if u.endswith("_1.fastq.gz")), "")
    r2 = next((u for u in urls if u.endswith("_2.fastq.gz")), "")
    run = rec.get("run_accession", "")
    return {
        "sample": re.sub(r"[^A-Za-z0-9_]", "_", run),
        "run_accession": run,
        "sample_accession": rec.get("sample_accession", ""),
        "study_accession": rec.get("study_accession", ""),
        "scientific_name": rec.get("scientific_name", ""),
        "country": country.strip() or "NA",
        "location": location.strip(),
        "collection_date": rec.get("collection_date", "") or "NA",
        "year": year_match.group(1) if year_match else "NA",
        "instrument_model": rec.get("instrument_model", ""),
        "library_layout": rec.get("library_layout", ""),
        "base_count": rec.get("base_count", ""),
        "fastq_1": f"https://{r1}" if r1 else "",
        "fastq_2": f"https://{r2}" if r2 else "",
    }


def subsample(records, n):
    """Keep at most n records, evenly spaced through the list (more diverse than the first n)."""
    if not n or len(records) <= n:
        return records
    step = len(records) / n
    return [records[int(i * step)] for i in range(n)]


def cmd_resolve(args):
    accessions = read_accessions(args.accessions)
    if not accessions:
        sys.exit(f"No accessions found in {args.accessions}")

    rows, seen = [], set()
    for acc in accessions:
        print(f"[resolve] {acc}", file=sys.stderr)
        try:
            records = fetch_runs(acc)
        except RuntimeError as exc:
            print(f"  WARNING: {exc}", file=sys.stderr)
            continue
        if not records:
            print(f"  WARNING: no runs returned for {acc}", file=sys.stderr)
        if args.max_per_accession and len(records) > args.max_per_accession:
            print(f"  {len(records)} runs found; keeping {args.max_per_accession} evenly spaced", file=sys.stderr)
            records = subsample(records, args.max_per_accession)
        for rec in records:
            row = clean_row(rec)
            if row["run_accession"] in seen:
                continue
            if row["library_layout"] != "PAIRED" or not (row["fastq_1"] and row["fastq_2"]):
                print(f"  skipped {row['run_accession']}: not paired-end with _1/_2 FASTQ", file=sys.stderr)
                continue
            seen.add(row["run_accession"])
            rows.append(row)
        time.sleep(0.3)  # be polite to the ENA API

    if not rows:
        sys.exit("No usable paired-end runs found.")

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / "metadata.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=METADATA_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with open(outdir / "samplesheet.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fastq_1", "fastq_2"])
        for r in rows:
            writer.writerow([r["sample"], r["fastq_1"], r["fastq_2"]])

    print(f"Wrote {len(rows)} runs -> {outdir / 'metadata.tsv'} and {outdir / 'samplesheet.csv'}")


def cmd_search(args):
    query = (
        f"tax_tree({TYPHI_TAXID}) AND library_layout=\"PAIRED\" "
        f"AND instrument_platform=\"ILLUMINA\""
    )
    if args.country:
        query += f" AND country=\"{args.country}*\""
    params = urllib.parse.urlencode(
        {"result": "read_run", "query": query, "fields": ",".join(FIELDS),
         "limit": args.max, "format": "tsv"}
    )
    records = parse_tsv(http_get(f"{ENA_API}/search?{params}"))
    if not records:
        sys.exit("No runs found. Try removing --country or check the ENA query syntax.")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("# Candidate S. Typhi runs from ENA search. Edit freely, then run: resolve\n")
        fh.write("accession\n")
        for rec in records:
            fh.write(rec["run_accession"] + "\n")
    print(f"Wrote {len(records)} candidate run accessions to {out}")
    for rec in records[:10]:
        print("  ", rec["run_accession"], rec.get("country", ""), rec.get("collection_date", ""))


def cmd_reference(args):
    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": args.accession, "rettype": "gbwithparts", "retmode": "text"}
    )
    text = http_get(f"{EFETCH}?{params}")
    if not text.startswith("LOCUS"):
        sys.exit(f"Unexpected response for {args.accession}:\n{text[:200]}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"Wrote reference {args.accession} -> {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("resolve", help="accession list -> metadata.tsv + samplesheet.csv")
    p.add_argument("--accessions", default="data/accession_list.tsv")
    p.add_argument("--out-dir", default="data/metadata")
    p.add_argument("--max-per-accession", type=int, default=None,
                   help="for study-level accessions, keep at most N runs (evenly spaced)")
    p.set_defaults(func=cmd_resolve)

    p = sub.add_parser("search", help="search ENA for candidate S. Typhi runs")
    p.add_argument("--country", default=None, help="e.g. Pakistan (prefix match)")
    p.add_argument("--max", type=int, default=50)
    p.add_argument("--out", default="data/accession_list.candidates.tsv")
    p.set_defaults(func=cmd_search)

    p = sub.add_parser("reference", help="download reference genome (GenBank)")
    p.add_argument("--accession", default="AL513382.1", help="default: S. Typhi CT18 chromosome")
    p.add_argument("--out", default="data/reference/CT18_AL513382.1.gbk")
    p.set_defaults(func=cmd_reference)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
