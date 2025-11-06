This directory contains a Python-based workflow for parsing, summarizing, and visualizing HLA typing results inspired by HISAT-genotype.
It can be called directly from the command line or integrated into workflow managers like Nextflow or Snakemake.
#Files
•	hla_parser.py – parses HISAT-genotype–style flat reports, trims allele names to the specified resolution (default = 2nd field), and summarizes results to a TSV.
•	hla_plots.py – generates generate pie charts or bar plots for the predicted alleles of each gene.
•	run_hla_workflow.py – main entrypoint; orchestrates parsing, summarization, and visualization for one or more input files.
Inputs
1.	Flat HLA Report (NA12877.hla.report.txt)
o	HISAT-genotype–style file listing ranked alleles with abundance values
(e.g., HLA-A*02:01 (abundance: 45.3%))
o	Required columns or fields:
	Locus name (e.g., HLA-A, -B, -C)
	Allele name (e.g., HLA-A*02:01:01:01)
	Abundance or score (%)
2.	Optional Configuration
o	-- --trim-level (1, 2, or allele-level)
o	--outdir for output tables and plots
## Run
python run_hla_workflow.py \
  NA12877.hla.report.txt\
  --outdir results/ \
  --trim-level 2

This will create:
results/
├── sample1_hla_summary.tsv
├── allele_abundance_plots/
│   ├── sample1_abundance.png
│   └── sample2_abundance.png
└── combined_heatmap.pdf
<img width="468" height="639" alt="image" src="https://github.com/user-attachments/assets/e49199de-78b6-455f-aef8-ec35d52fcf1e" />
