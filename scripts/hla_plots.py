#!/usr/bin/env python3
# hla_plots.py

import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# we import the parser + trimming from the first script
from hla_parser import parse_ranked_lines_to_df, trim_allele_fields


def make_plots(
    df,
    cutoff=10.0,
    trim_level=2,
    out_prefix="hla_plots",
    plot_type="both",
    pdf_out="hla_plots_all.pdf",
):
    """
    Make per-gene plots from parsed HLA DF.

    Parameters
    ----------
    df : pandas.DataFrame
        Must have columns: gene, allele_full, abundance
    cutoff : float
        Only plot alleles with abundance >= cutoff
    trim_level : int
        How many fields to show in plot labels (2 -> A*03:01)
    out_prefix : str
        Prefix for PNG files, e.g. 'hla_plots' -> hla_plots_A.png
    plot_type : str
        'bar', 'pie', or 'both'
    pdf_out : str
        Name of combined PDF file to write
    """

    # we will collect all figures into one PDF
    with PdfPages(pdf_out) as pdf:
        genes = df["gene"].unique()

        for g in genes:
            sub = df[(df["gene"] == g) & (df["abundance"] >= cutoff)].copy()
            if sub.empty:
                continue

            # label alleles in the plot using the requested trim resolution
            sub["allele_label"] = sub["allele_full"].apply(
                lambda a: trim_allele_fields(a, trim_level)
            )

            # ---------------- BAR PLOT ----------------
            if plot_type in ("bar", "both"):
                fig_bar = plt.figure()
                plt.bar(sub["allele_label"], sub["abundance"])
                plt.title(f"HLA-{g} predicted alleles (≥ {cutoff}% abundance)")
                plt.xlabel("Allele")
                plt.ylabel("Abundance (%)")
                plt.xticks(rotation=45, ha="right")
                plt.tight_layout()

                # save individual PNG
                fig_bar.savefig(f"{out_prefix}_{g}.png")
                # also add to combined PDF
                pdf.savefig(fig_bar)
                plt.close(fig_bar)

            # ---------------- PIE CHART ----------------
            if plot_type in ("pie", "both"):
                fig_pie = plt.figure()
                plt.pie(
                    sub["abundance"],
                    labels=sub["allele_label"],
                    autopct="%1.1f%%"
                )
                plt.title(f"HLA-{g} predicted alleles (≥ {cutoff}% abundance)")
                plt.tight_layout()

                # save individual PNG
                fig_pie.savefig(f"{out_prefix}_{g}_pie.png")
                # also add to combined PDF
                pdf.savefig(fig_pie)
                plt.close(fig_pie)


def main():
    ap = argparse.ArgumentParser(
        description="Generate per-gene HLA abundance plots (bar/pie/both) and one combined PDF."
    )
    ap.add_argument("infile", help="HISAT-genotype text output")
    ap.add_argument(
        "--cutoff",
        type=float,
        default=10.0,
        help="Only plot alleles ≥ this abundance (default 10).",
    )
    ap.add_argument(
        "--trim-level",
        type=int,
        default=2,
        help="How many fields to show in labels (default 2 = A*03:01).",
    )
    ap.add_argument(
        "--out-prefix",
        default="hla_plots",
        help="Prefix for per-gene PNGs (default: hla_plots).",
    )
    ap.add_argument(
        "--pdf-out",
        default="hla_plots_all.pdf",
        help="Combined PDF filename (default: hla_plots_all.pdf).",
    )
    ap.add_argument(
        "--plot-type",
        choices=["bar", "pie", "both"],
        default="both",
        help="Which plots to make: bar, pie, or both (default).",
    )
    args = ap.parse_args()

    # 1. parse file into DF
    df = parse_ranked_lines_to_df(args.infile)

    # 2. make plots according to user choices
    make_plots(
        df,
        cutoff=args.cutoff,
        trim_level=args.trim_level,
        out_prefix=args.out_prefix,
        plot_type=args.plot_type,
        pdf_out=args.pdf_out,
    )


if __name__ == "__main__":
    main()


# can run the code with command: python hla_plots.py hisat_output.txt --cutoff 10 --trim-level 2 --out-prefix donor1_plots