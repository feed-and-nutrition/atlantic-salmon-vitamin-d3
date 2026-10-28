# Vitamin D3-related genes in Atlantic salmon

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23711905.svg)](https://doi.org/10.5281/zenodo.23711905)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository holds the analysis code behind the manuscript **"Genomic identification and multi-tissue
transcriptomic profiling of vitamin D3-related genes in Atlantic salmon (*Salmo salar*)"** by Takaya Saito,
Øystein Sæle, Pauline Wischhusen and Anne-Catrin Adam (Institute of Marine Research, Bergen, Norway). An earlier
version is available as a preprint on bioRxiv ([doi:10.64898/2026.03.03.709343](https://doi.org/10.64898/2026.03.03.709343)).

The scripts:

- **identify** 16 vitamin D3-related genes in the Atlantic salmon genome (RefSeq annotation RS_2025_05),
- **build protein trees** to tell orthologs and paralogs apart,
- **search six fish genomes** for the vitamin D-binding protein gene (*gc*),
- **profile expression** in a 14-tissue atlas and seven independent RNA-seq studies,
- and **draw the figures and write the supplementary tables** of the manuscript.

## Quick start

Most figures and tables can be rebuilt from the files in this repository, without downloading any sequencing
data.

```bash
git clone https://github.com/feed-and-nutrition/atlantic-salmon-vitamin-d3.git
cd atlantic-salmon-vitamin-d3
Rscript scripts/run_all.R
```

`run_all.R` runs each step in turn and prints a summary at the end. Steps that need the genome-wide expression
tables (see [What needs the large data](#what-needs-the-large-data)) are reported as failed and skipped, and the
others still run. Figures and Supplementary Tables (Excel files) are written to `manuscript/v6/`.

To write them somewhere else:

```bash
VD3_FIG_DIR=results Rscript scripts/run_all.R
```

## Installation

### R and packages

R 4.6.1 was used. Install the packages once:

```r
install.packages(c("tidyverse", "cowplot", "systemfonts", "ragg", "ggrepel", "pheatmap", "ape",
                   "xml2", "jsonlite", "BiocManager"))
BiocManager::install(c("ggtree", "Biostrings"))
```

### Command-line tools

Only some steps need these. The versions used for the manuscript are given in brackets.

| Tool | Used by | Needed for |
|---|---|---|
| python3 (standard library only) | most table scripts | writing the Excel tables |
| `muscle5` (5.3) and `iqtree2` (2.0.7) | 01_0, 06_4 (muscle5 only) | rebuilding the protein trees, codon alignments of the paralog pairs |
| BLAST+ (2.12.0) | 05_8, 06_1, 07_1, 07_3, 07_4, 07_5 | sequence comparisons and genome searches |
| `wgsim` (1.19.2) and `salmon` (1.10.2) | 06_2, 06_3 | the paralog read simulation and the pike expression |
| `fastq-dump` (sra-tools) | 06_3 | downloading the pike reads |
| `yn00` (PAML 4.10.10) | 06_4 | dN/dS of the paralog pairs (only if `divergence.tsv` is missing; set `PAML_YN00` to its path if it is not on the `PATH`) |

The tools must be on your `PATH`. With conda, for example:

```bash
conda install -c conda-forge -c bioconda "muscle>=5" iqtree blast wgsim salmon sra-tools paml
ln -s "$(command -v muscle)" "$CONDA_PREFIX/bin/muscle5"   # the scripts call MUSCLE 5 as muscle5
```

### Fonts (optional)

The figures use the fonts Arimo and Nimbus Roman, so that they look the same on every computer. If they are not
installed, a warning is printed and the system default fonts are used instead.

## What is in the repository

```
scripts/        all analysis scripts, run from the repository root (or from scripts/)
ncbi_samples/   list of the RNA-seq studies and cached SRA run information
salmon_quant/   sample sheets and sample tables of each RNA-seq study
msa/            protein sequences, alignments and IQ-TREE results for each tree
vd3_genes/      the curated gene list and the expression, annotation and gc search results
hk_genes/       housekeeping genes used as an expression reference
gc_blast/       BLASTp results of the first GC search (NCBI web BLAST)
d3pathway/      Figure 1 (and the draw.io drawing of earlier versions)
```

The curated gene list is `vd3_genes/salmon_vd3_orthologs_v1.tsv`. Expression values of the 16 genes are in
`vd3_genes/tpm_vd3_genes.tsv` (atlas) and `vd3_genes/tpm_vd3_genes_studies.tsv` (all studies).

Intermediate results are included, so any step can be rerun on its own. Large files are not: the genome-wide
expression tables (`*.tsv.gz`) and downloaded genomes and BLAST databases (`annotation/`).

## Which script makes which figure or table

| Output | Script | Runs from this repository |
|---|---|---|
| Figure 1 | `09_1_pathway_figure.R` | yes |
| Figure 2, Figure S4 | `01_phylogenetic_tree.R` (trees built by `01_0_build_trees.R`) | yes |
| Figure 3, Figure S6 | `04_2_vd3_genes_heatmap_split.R` | yes |
| Figure 4, Figure S5 | `05_5_vd3_genes_studies_fig.R` | yes |
| Figure S1 | `03_2_salmon_quant_pca.R` | needs the large data |
| Figure S2 | `05_4_atlas_heart_liver.R` | needs the large data |
| Figure S3 | `02_housekeeping_heatmap.R` | yes |
| Table S3 | `07_3_gc_genome_search.R` | yes, but downloads about 3.5 GB |
| Table S6 | `07_2_synteny_table.R` (from the results of `07_1_annotation_check.R` and `07_4_vdr_synteny.R`) | yes |
| Table S7 | `08_1_ligand_contacts.R` | yes |
| Table S8 | `04_1_vd3_genes_heatmap.R` | yes |
| Table S9 | `06_4_paralog_fate.R` (with the results of `06_1_paralog_similarity.R`, `06_2_paralog_simulation.R` and `06_3_pike_expression.R`) | yes |
| Table S10 | `05_6_study_table.R` | yes |
| Table S11 | `05_7_kidney_detection.R` (with the results of `07_5_cyp27b1_genome_search.R`) | needs the large data |
| Table S12 | `05_8_related_genes.R` | needs the large data |

Tables S1, S2, S4 and S5 were compiled by hand from the gene search and are provided with the manuscript.

Six scripts are not part of `run_all.R`, because they download large files to `annotation/` and can take an hour
or more: `06_2_paralog_simulation.R`, `06_3_pike_expression.R` (about 9 GB of pike reads), `07_1_annotation_check.R`,
`07_3_gc_genome_search.R`, `07_4_vdr_synteny.R` and `07_5_cyp27b1_genome_search.R` (both use the downloads of
`07_3`). Their results are included, so the tables that use them can still be rebuilt.

`01_0_build_trees.R` reuses the trees in `msa/` and only rebuilds a tree when its sequence file is newer than the
tree. Without `muscle5` and `iqtree2` it may fail, and `01_phylogenetic_tree.R` still draws the included trees.

## What needs the large data

The genome-wide expression tables are not included because of their size. They are made from the public reads
in three steps.

1. **Quantify the reads.** Download the reads of each study from the SRA (run accessions in Supplementary
   Table S10 and `ncbi_samples/runinfo/`). In the manuscript, reads were checked with FastQC (0.11.9) and MultiQC
   (1.17), trimmed with TrimGalore! (0.6.7), quantified against the RefSeq RS_2025_05 transcripts with Salmon
   (1.4.0, options `-l A --seqBias --gcBias`) and summarized to genes with the tximport R package (1.40.0). Save
   the result of each study as `txi_gene.rds`, one directory per study.
2. **Export the gene tables.** Run `05_1_export_studies_server.R` on the computer that holds these files. It
   needs base R only. Each study's directory name is in the `location` column of `ncbi_samples/sample_info.txt`:

   ```bash
   VD3_MAPPED_DIR=/path/to/study/directories Rscript scripts/05_1_export_studies_server.R
   ```

   It writes `salmon_quant/<study>/tpm_gene.tsv.gz`, `counts_gene.tsv.gz` and `samples.tsv`.
3. **Run everything.** Copy `salmon_quant/` back into the repository if needed and run `Rscript scripts/run_all.R`
   again. All steps should now succeed.

## Settings

Paths, tissue and study lists, excluded samples and genes, and plot helpers are in `scripts/00_config.R`. These
environment variables change the locations:

| Variable | Meaning | Default |
|---|---|---|
| `VD3_ROOT` | repository root | the working directory |
| `VD3_FIG_DIR` | where figures and Supplementary Tables are written | `manuscript/v6` |
| `VD3_TXI_DIR` | directory with the atlas `txi_gene.rds` (used by 00_1 only when `tpm_gene.tsv.gz` is missing) | not set |
| `VD3_MAPPED_DIR` | directory that holds the study directories (used by 05_1 only) | the working directory |

## Troubleshooting

- **"Font 'Arimo' not found"**: a warning only. Install the fonts or ignore it.
- **A step is reported as failed by `run_all.R`**: the summary at the end gives the reason. A missing
  `tpm_gene.tsv.gz` means that the step needs the large data.
- **A download fails**: some scripts fetch sequences and run information from NCBI the first time. Check the
  internet connection and rerun. Downloaded files are cached and not fetched again.

## Citation

If you use this code, please cite the manuscript and this repository:

> Saito T, Sæle Ø, Wischhusen P, Adam A-C. Analysis scripts for vitamin D3-related genes in Atlantic salmon.
> Zenodo; 2026. [doi:10.5281/zenodo.23711905](https://doi.org/10.5281/zenodo.23711905)

GitHub's "Cite this repository" button gives the same reference in other formats.

## License

The code is released under the MIT License (see [`LICENSE`](LICENSE)).
