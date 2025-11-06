#!/usr/bin/env python3
"""
This script was inspired by and adapted from HISAT-genotype’s internal post-processing scripts.

Key adaptations:
- The original allele-tree flattening logic is simplified as `trim_allele_fields()`.
- The result summarization and CSV export are reimplemented as `summarize_for_table()`
  using pandas for easier table handling.
- Command-line parsing and output generation are redesigned for standalone use
  with HISAT-genotype flat report files (ranked allele abundance format).

All HISAT-genotype–specific dependencies (e.g., `hisatgenotype_typing_common`, `hg_args`)
and hierarchical allele tree traversal have been removed for simplicity.


"""
# hla_parser.py
import argparse
import re
import pandas as pd

def trim_allele_fields(allele_full, fields_to_keep=2):
    if "*" not in allele_full:
        return allele_full
    gene, rest = allele_full.split("*", 1)
    parts = rest.split(":")
    fields_to_keep = max(1, min(fields_to_keep, len(parts)))
    kept = ":".join(parts[:fields_to_keep])
    return f"{gene}*{kept}"

def parse_ranked_lines_to_df(path):
    ranked_re = re.compile(
        r'^\s*\d+\s+ranked\s+([A-Za-z0-9*:+]+)\s+\(abundance:\s+([0-9.]+)%',
        re.IGNORECASE
    )
    rows = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            m = ranked_re.match(line)
            if not m:
                continue

            allele_full_raw = m.group(1)
            abundance = float(m.group(2))

            if "*" not in allele_full_raw:
                continue
            gene = allele_full_raw.split("*", 1)[0]

            flag_match = re.match(r'^(.+?)([NLSCAQ]+)$', allele_full_raw)
            if flag_match and "*" in flag_match.group(1):
                core_allele = flag_match.group(1)
                trailing_flags = flag_match.group(2)
            else:
                core_allele = allele_full_raw
                trailing_flags = ""

            allele_2field = trim_allele_fields(core_allele, 2)

            rows.append({
                "gene": gene,
                "allele_full": core_allele,
                "allele_2field": allele_2field,
                "abundance": abundance,
                "flags": trailing_flags
            })

    return pd.DataFrame(rows)

def summarize_for_table(df, cutoff=10.0):
    keep = df[df["abundance"] >= cutoff].copy()
    keep = keep.sort_values(["gene", "abundance"], ascending=[True, False])
    table = keep[["gene", "allele_full", "allele_2field", "abundance", "flags"]]
    return table

def main():
    ap = argparse.ArgumentParser(
        description="Parse HISAT-genotype report and produce TSV of alleles ≥ cutoff."
    )
    ap.add_argument("infile", help="HISAT-genotype text output")
    ap.add_argument("--cutoff", type=float, default=10.0,
                    help="Abundance cutoff in percent (default 10)")
    ap.add_argument("--tsv-out", default="hla_summary.tsv",
                    help="Output TSV filename (default hla_summary.tsv)")
    args = ap.parse_args()

    df = parse_ranked_lines_to_df(args.infile)
    table = summarize_for_table(df, cutoff=args.cutoff)
    table.to_csv(args.tsv_out, sep="\t", index=False)
    print(table.to_string(index=False))

if __name__ == "__main__":
    main()

## can run the code with command:  python hla_parser.py hisat_output.txt --cutoff 10 --tsv-out donor1.tsv
