#!/usr/bin/env python3
"""
parse_genbank_phage_virome_fetch_cached.py
==========================================

WHAT THIS SCRIPT DOES (Beginner-friendly overview)
--------------------------------------------------
This script helps you *download* GenBank records from NCBI **by accession or by search query**,
**stores them locally** (your "bank"), and then **summarizes** them with a focus on
**bacteriophages and viromes**. It can also work on local files you already have.

It prints a small table with, for each file:
- phage_name (taken from the filename)
- genome_length (total length in base pairs)
- num_genes
- num_trnas
- num_proteins (CDS features)
- n_hypotheticals (how many CDS products contain the phrase "hypothetical protein")

It can also save outputs:
- CSV (spreadsheet)
- JSON (machine-readable)
- TXT list of all "hypothetical protein" product lines

BEGINNER NOTES
--------------
- **NCBI Entrez API** lets us download records for free. Please provide your email with --entrez-email.
- For faster, more reliable requests, you can also provide an --api-key (optional but recommended).
- The **cache** means: if a file for a given accession already exists in your save folder,
  we will **reuse it and skip re-downloading** (unless you use --force). This saves time and quota.

EXAMPLES
--------
# 1) Fetch by accessions (download -> cache -> summarize)
python parse_genbank_phage_virome_fetch_cached.py \
  --entrez-email you@school.edu \
  --acc NC_001416 MN508617 \
  --save-dir ./bank \
  --csv summary.csv

# 2) Search NCBI, fetch top 5, save + summarize
python parse_genbank_phage_virome_fetch_cached.py \
  --entrez-email you@school.edu \
  --search "bacteriophage complete genome" \
  --max-hits 5 \
  --save-dir ./bank \
  --json summary.json

# 3) Use local files + turn off the phage filter
python parse_genbank_phage_virome_fetch_cached.py \
  --no-filter local1.gbk local2.gbk

# 4) Pick local files interactively from a folder
python parse_genbank_phage_virome_fetch_cached.py --pick ./bank

# 5) Use a GUI file picker (Linux may need: sudo apt install python3-tk)
python parse_genbank_phage_virome_fetch_cached.py --gui

SETUP
-----
pip install biopython
# (Optional on Linux for GUI picker) sudo apt install python3-tk
"""

# --------------- 1) IMPORTS (built-in + Biopython) ---------------
from __future__ import annotations

import argparse   # Read command-line options like --acc, --search, etc.
import csv        # Write CSV files
import json       # Write JSON files
import os         # Work with files and folders
import re         # Search text using regular expressions
import sys        # For error messages and exiting
import time       # For polite delays between NCBI requests
from typing import Dict, List, Iterable

# Biopython pieces we use:
from Bio import SeqIO          # To parse GenBank files
from Bio import Entrez         # To talk to NCBI (download/search)


# --------------- 2) SIMPLE PHAGE/VIROME FILTER HELPERS ---------------
# We want to *prefer* bacteriophages/viromes. These keywords help us detect that.
PHAGE_WORDS = {"phage", "bacteriophage"}
VIROME_WORDS = {"virome", "metagenome"}
TAXON_HINTS = {"viruses", "caudoviricetes"}  # common phage taxonomy words

def _safe_lower(s: str | None) -> str:
    """Return a lowercase string; if s is None, return empty string."""
    return (s or "").lower()

def record_looks_like_phage_or_virome(record) -> bool:
    """
    Heuristic (educated guess) to decide if a GenBank record is a bacteriophage/virome.
    We check several places for hints: description, organism, taxonomy, keywords, source.
    This is imperfect but works well for a lot of student use cases.
    """
    # 1) Description text
    desc = _safe_lower(getattr(record, "description", ""))
    if any(w in desc for w in PHAGE_WORDS | VIROME_WORDS):
        return True

    # 2) Record annotations: organism name, taxonomy, keywords
    ann = getattr(record, "annotations", {}) or {}

    organism = _safe_lower(ann.get("organism", ""))
    if any(w in organism for w in PHAGE_WORDS):
        return True

    taxonomy = [t.lower() for t in ann.get("taxonomy", []) if isinstance(t, str)]
    if any(h in taxonomy for h in TAXON_HINTS):
        return True

    keywords = [k.lower() for k in ann.get("keywords", []) if isinstance(k, str)]
    if any(k in PHAGE_WORDS or k in VIROME_WORDS for k in keywords):
        return True

    # 3) "source" feature may contain organism/host/isolation_source strings
    for feat in getattr(record, "features", []):
        if getattr(feat, "type", "").lower() == "source":
            quals = getattr(feat, "qualifiers", {}) or {}
            src_orgs = [_safe_lower(x) for x in quals.get("organism", [])]
            hosts = [_safe_lower(x) for x in quals.get("host", [])]
            iso = [_safe_lower(x) for x in quals.get("isolation_source", [])]
            txt = " ".join(src_orgs + hosts + iso)
            if any(w in txt for w in PHAGE_WORDS | VIROME_WORDS):
                return True

    # If we saw no hints, return False
    return False


# --------------- 3) CORE SUMMARIZATION OF ONE FILE ---------------
def summarize_file(filepath: str, filter_to_phage: bool = True) -> Dict:
    """
    Read one GenBank file (may have multiple records) and compute counts.
    If filter_to_phage is True, only include records that look like phage/virome.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")

    # phage_name = file name without extension
    phage_name = os.path.splitext(os.path.basename(filepath))[0]

    # Initialize counters we will add to as we scan records
    genome_length = 0
    gene_count = 0
    trna_count = 0
    protein_count = 0
    hypothetical_products: List[str] = []

    # Pre-compile a regex to detect "hypothetical protein" (case-insensitive)
    hypo_re = re.compile(r"\bhypothetical\s+protein\b", re.IGNORECASE)

    # Read the GenBank file into a list of records
    records = list(SeqIO.parse(filepath, "genbank"))
    if not records:
        raise RuntimeError(f"No GenBank records found in {filepath!r}")

    included = 0  # how many records we actually counted
    for record in records:
        # If filtering is ON and this record doesn't look like phage/virome, skip it
        if filter_to_phage and not record_looks_like_phage_or_virome(record):
            continue

        included += 1
        genome_length += len(record.seq)

        # Walk through each annotated feature and count types of interest
        for feature in record.features:
            ftype = (feature.type or "").lower()  # normalize to lowercase

            if ftype == "gene":
                gene_count += 1

            elif ftype == "trna":
                trna_count += 1

            elif ftype == "cds":
                protein_count += 1

                # Look for "hypothetical protein" in product qualifiers
                for prod in feature.qualifiers.get("product", []):
                    if prod and hypo_re.search(prod):
                        # Collect useful IDs if available for context
                        locus = ",".join(feature.qualifiers.get("locus_tag", [])) or ""
                        prot_id = ",".join(feature.qualifiers.get("protein_id", [])) or ""
                        tag = locus or prot_id
                        hypothetical_products.append(f"{tag} :: {prod}" if tag else prod)

    if filter_to_phage and included == 0:
        raise RuntimeError(f"{filepath!r}: no records looked like phage/virome (use --no-filter to include all).")

    # Return a dictionary of results
    return {
        "phage_name": phage_name,
        "genome_length": genome_length,
        "num_genes": gene_count,
        "num_trnas": trna_count,
        "num_proteins": protein_count,
        "n_hypotheticals": len(hypothetical_products),
        "hypothetical_products": hypothetical_products,
    }


# --------------- 4) PRETTY TABLE FOR THE CONSOLE ---------------
def format_table(rows: List[Dict]) -> str:
    """
    Build a simple aligned text table (so it looks neat in the terminal).
    """
    headers = ["phage_name", "genome_length", "num_genes", "num_trnas", "num_proteins", "n_hypotheticals"]
    # Make rows of strings
    data = [[
        str(r.get("phage_name", "")),
        str(r.get("genome_length", "")),
        str(r.get("num_genes", "")),
        str(r.get("num_trnas", "")),
        str(r.get("num_proteins", "")),
        str(r.get("n_hypotheticals", "")),
    ] for r in rows]

    # Compute column widths
    widths = [max(len(h), *(len(row[i]) for row in data)) for i, h in enumerate(headers)]

    def fmt(cols: Iterable[str]) -> str:
        return " | ".join(c.ljust(widths[i]) for i, c in enumerate(cols))

    sep = "-+-".join("-" * w for w in widths)
    lines = [fmt(headers), sep]
    for row in data:
        lines.append(fmt(row))
    return "\n".join(lines)


# --------------- 5) CACHING HELPERS (SKIP RE-DOWNLOADS) ---------------
def path_for_accession(save_dir: str, acc: str) -> str:
    """Return a standard local filepath for a given accession, ending in .gbk"""
    os.makedirs(save_dir, exist_ok=True)
    return os.path.join(save_dir, f"{acc}.gbk")

def file_contains_genbank_records(path: str) -> bool:
    """
    Try to parse at least one record from a file.
    If we can read one record, we consider it a valid cached file.
    """
    try:
        for _rec in SeqIO.parse(path, "genbank"):
            return True  # success: at least one record found
        return False      # parsed but zero records
    except Exception:
        return False       # parsing failed

def is_cached(path: str) -> bool:
    """A file is 'cached' if it exists, is not empty, and parses as GenBank."""
    return os.path.isfile(path) and os.path.getsize(path) > 0 and file_contains_genbank_records(path)


# --------------- 6) NCBI ENTREZ (DOWNLOAD / SEARCH) ---------------
def entrez_setup(email: str, api_key: str | None = None) -> None:
    """
    Tell Biopython's Entrez your email (required by NCBI) and API key (optional).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

def fetch_by_accession(acc: str, save_dir: str, delay: float = 0.34, force: bool = False) -> str:
    """
    Download one accession from NCBI nuccore in GenBank format and save it.
    - Uses a polite delay between calls (0.34s ~ 3 req/sec without API key).
    - If the file already exists and looks valid, we SKIP downloading unless --force.
    Returns the local file path.
    """
    out_path = path_for_accession(save_dir, acc)

    # CACHE CHECK: if we already have a valid file and not forcing, reuse it
    if not force and is_cached(out_path):
        print(f"[cache] Using cached file for {acc}: {out_path}")
        return out_path

    # Otherwise, fetch from NCBI
    print(f"[fetch] Downloading {acc} from NCBI ...")
    with Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text") as handle:
        text = handle.read()

    # Basic sanity check
    if not text.strip():
        raise RuntimeError(f"No content returned for accession {acc}")

    # Write to disk
    os.makedirs(save_dir, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(text)

    # Polite pause to respect NCBI servers (adjust if you have an API key)
    time.sleep(delay)

    # Final validation: ensure we can parse it
    if not is_cached(out_path):
        raise RuntimeError(f"Downloaded file for {acc} did not parse as GenBank: {out_path}")

    print(f"[saved] {out_path}")
    return out_path

def search_and_fetch(query: str, max_hits: int, save_dir: str, delay: float = 0.34, force: bool = False) -> List[str]:
    """
    Search NCBI nuccore for 'query', then fetch up to max_hits results.
    We fetch each accession individually so cache checking works cleanly.
    Returns a list of local file paths.
    """
    print(f"[search] '{query}' (max {max_hits})")
    os.makedirs(save_dir, exist_ok=True)

    with Entrez.esearch(db="nuccore", term=query, retmax=max_hits) as h:
        rec = Entrez.read(h)

    ids = rec.get("IdList", [])
    if not ids:
        print("[search] No results.")
        return []

    out_paths: List[str] = []
    for acc in ids:
        try:
            p = fetch_by_accession(acc, save_dir, delay=delay, force=force)
            out_paths.append(p)
        except Exception as e:
            sys.stderr.write(f"[WARN] Failed fetching {acc}: {e}\n")

    return out_paths


# --------------- 7) OPTIONAL: LOCAL FILE PICKERS ---------------
def choose_files_from_dir(directory: str) -> List[str]:
    """
    List .gb/.gbk/.genbank files in a directory and ask the user to choose by number.
    Returns absolute paths to the chosen files.
    """
    if not os.path.isdir(directory):
        raise RuntimeError(f"Not a directory: {directory}")
    cand = [f for f in os.listdir(directory) if f.lower().endswith((".gb", ".gbk", ".genbank"))]
    cand.sort()
    if not cand:
        raise RuntimeError(f"No GenBank files found in {directory!r}.")

    print("\nSelect file(s) by typing numbers separated by spaces. Press Enter to select all.\n")
    for i, name in enumerate(cand, 1):
        print(f"{i:2d}) {name}")
    print()

    sel = input("Your choice (e.g., '1 3 5', or just Enter for all): ").strip()
    if not sel:
        chosen = cand
    else:
        idx = []
        for token in sel.split():
            n = int(token)
            if not (1 <= n <= len(cand)):
                raise RuntimeError(f"Choice out of range: {n}")
            idx.append(n - 1)
        chosen = [cand[i] for i in idx]

    return [os.path.abspath(os.path.join(directory, x)) for x in chosen]

def pick_files_gui() -> List[str]:
    """
    GUI file picker using Tkinter (works on most systems; on Linux you may need: sudo apt install python3-tk).
    Returns a list of selected file paths.
    """
    import tkinter as tk
    from tkinter import filedialog
    root = tk.Tk(); root.withdraw()
    paths = filedialog.askopenfilenames(
        title="Select GenBank files",
        filetypes=[("GenBank files", "*.gb *.gbk *.genbank"), ("All files", "*.*")],
    )
    root.update(); root.destroy()
    return list(paths)


# --------------- 8) MAIN: GLUE EVERYTHING TOGETHER ---------------
def main():
    # 8a) Define command-line options for the user
    ap = argparse.ArgumentParser(
        description="Fetch (with caching) and summarize bacteriophage/virome GenBank files."
    )

    # Local inputs (you can pass files directly, or pick interactively, or via GUI)
    ap.add_argument("genbank_files", nargs="*", help="Local GenBank files (.gb/.gbk/.genbank)")
    ap.add_argument("--pick", metavar="DIR", help="Pick local files interactively from a directory.")
    ap.add_argument("--gui", action="store_true", help="Pick local files via a GUI dialog (tkinter).")

    # Remote selection & download (NCBI)
    ap.add_argument("--entrez-email", help="Required for NCBI requests: your contact email.")
    ap.add_argument("--api-key", help="Optional NCBI API key for higher rate limits.")
    ap.add_argument("--acc", nargs="+", help="One or more accessions to fetch (space-separated).")
    ap.add_argument("--search", help="NCBI esearch query (nuccore).")
    ap.add_argument("--max-hits", type=int, default=5, help="Max hits to fetch when using --search (default: 5).")
    ap.add_argument("--save-dir", default="./downloaded_gbk", help="Directory to save fetched files (default: ./downloaded_gbk).")
    ap.add_argument("--force", action="store_true", help="Force re-download even if a cached file exists.")

    # Behavior / outputs
    ap.add_argument("--no-filter", action="store_true", help="Do NOT filter to phage/virome; include all records.")
    ap.add_argument("--csv", help="Write a CSV summary to this path.")
    ap.add_argument("--json", help="Write a JSON summary to this path.")
    ap.add_argument("--list-hypo", help="Write all hypothetical-product lines grouped by file to this TXT path.")

    args = ap.parse_args()

    # 8b) Build a file list from local choices first
    files: List[str] = []
    files.extend(args.genbank_files or [])

    if args.pick:
        try:
            files.extend(choose_files_from_dir(args.pick))
        except Exception as e:
            sys.stderr.write(f"[ERROR] {e}\n")
            sys.exit(1)

    if args.gui:
        try:
            files.extend(pick_files_gui())
        except Exception as e:
            sys.stderr.write(f"[ERROR] GUI selection failed: {e}\n")
            sys.exit(1)

    # 8c) Handle remote fetching (with caching)
    fetched_paths: List[str] = []
    if args.acc or args.search:
        if not args.entrez_email:
            sys.stderr.write("[ERROR] --entrez-email is required when using --acc or --search.\n")
            sys.exit(2)
        entrez_setup(args.entrez_email, args.api_key)

        # For each accession, download (or reuse cached) and collect file paths
        if args.acc:
            for acc in args.acc:
                try:
                    p = fetch_by_accession(acc, args.save_dir, delay=0.12 if args.api_key else 0.34, force=args.force)
                    fetched_paths.append(p)
                except Exception as e:
                    sys.stderr.write(f"[WARN] Failed fetching {acc}: {e}\n")

        # For search queries: find IDs and fetch them (with caching)
        if args.search:
            fetched_paths.extend(
                search_and_fetch(args.search, args.max_hits, args.save_dir, delay=0.12 if args.api_key else 0.34, force=args.force)
            )

    # Add fetched paths to the full list
    files.extend(fetched_paths)

    # 8d) Deduplicate file list while preserving order (avoid analyzing the same file twice)
    seen = set(); ordered: List[str] = []
    for f in files:
        if f not in seen:
            seen.add(f); ordered.append(f)
    files = ordered

    # 8e) If we have no files at all, show a helpful error
    if not files:
        ap.error("No input files. Provide local files or use --acc/--search with --entrez-email.")

    # 8f) Summarize each file and collect rows for output
    rows: List[Dict] = []
    hypo_lines: List[str] = []
    had_error = False
    filter_on = not args.no_filter  # True by default (filter to phage/virome)

    for fp in files:
        try:
            s = summarize_file(fp, filter_to_phage=filter_on)
            rows.append(s)
            if s["hypothetical_products"]:
                hypo_lines.append(f"# {s['phage_name']}")
                for line in s["hypothetical_products"]:
                    hypo_lines.append(f"- {line}")
                hypo_lines.append("")
        except Exception as e:
            had_error = True
            sys.stderr.write(f"[ERROR] {e}\n")

    # 8g) Print the console table
    if rows:
        print(format_table(rows))

    # 8h) Optional file outputs
    if args.csv and rows:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["phage_name", "genome_length", "num_genes", "num_trnas", "num_proteins", "n_hypotheticals"])
            for r in rows:
                w.writerow([
                    r["phage_name"], r["genome_length"], r["num_genes"],
                    r["num_trnas"], r["num_proteins"], r["n_hypotheticals"]
                ])
        print(f"[wrote] {args.csv}")

    if args.json and rows:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(rows, fh, indent=2)
        print(f"[wrote] {args.json}")

    if args.list_hypo and hypo_lines:
        with open(args.list_hypo, "w", encoding="utf-8") as fh:
            fh.write("\n".join(hypo_lines).rstrip() + "\n")
        print(f"[wrote] {args.list_hypo}")

    if had_error and not rows:
        # If *all* files failed, exit with non-zero status for CI/automation use
        sys.exit(1)


# --------------- 9) STANDARD PYTHON ENTRY POINT ---------------
if __name__ == "__main__":
    main()
