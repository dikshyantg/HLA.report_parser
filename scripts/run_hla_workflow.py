#!/usr/bin/env python3
# run_hla_workflow.py
import argparse
import os
from pathlib import Path
from hla_parser import parse_ranked_lines_to_df, summarize_for_table
from hla_plots import make_plots

def main():
    ap = argparse.ArgumentParser(
        description="End-to-end HLA report workflow: TSV + plots."
    )
    ap.add_argument("infile", help="HISAT-genotype text output")
    ap.add_argument(
        "--cutoff", type=float, default=10.0,
        help="Abundance cutoff for both table and plots (default 10)"
    )
    ap.add_argument(
        "--trim-level", type=int, default=2,
        help="Trim alleles to this many fields for plot labels"
    )
    ap.add_argument(
        "--tsv-out", default="hla_summary.tsv",
        help="Output TSV filename (will be placed inside --outdir if given)"
    )
    ap.add_argument(
        "--out-prefix", default="hla_plots",
        help="Output prefix for plots (basename; will live in --outdir)"
    )
    ap.add_argument(
        "--outdir", default=".",
        help="Output directory (default current directory)"
    )
    args = ap.parse_args()

    # ensure outdir exists
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. parse once
    df = parse_ranked_lines_to_df(args.infile)

    # 2. write TSV
    table = summarize_for_table(df, cutoff=args.cutoff)
    tsv_path = outdir / args.tsv_out
    table.to_csv(tsv_path, sep="\t", index=False)
    print(table.to_string(index=False))

    # 3. make plots
    # build full prefix inside outdir
    plot_prefix = outdir / args.out_prefix
    make_plots(
        df,
        cutoff=args.cutoff,
        trim_level=args.trim_level,
        out_prefix=str(plot_prefix)
    )

if __name__ == "__main__":
    main()
