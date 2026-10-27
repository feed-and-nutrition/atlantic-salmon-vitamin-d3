# Vitamin D3-related genes in Atlantic salmon

Analysis scripts for Saito T, Sæle Ø, Wischhusen P, Adam A-C, "Genomic identification and multi-tissue
transcriptomic profiling of vitamin D3-related genes in Atlantic salmon (*Salmo salar*)" (preprint of an earlier
version: bioRxiv, doi:10.64898/2026.03.03.709343). They identify the 16 genes, build the protein trees, profile expression in a
14-tissue atlas and seven replicate RNA-seq studies, and write Figures 2 to 4, Figures S1 to S6 and the
Supplementary Tables.

## Requirements

- R (4.6.1 was used) with tidyverse, cowplot, systemfonts, ragg, ggrepel, pheatmap, ape, ggtree, Biostrings
  (Bioconductor), xml2 and jsonlite.
- python3 (standard library only), for the Excel tables (`scripts/tsv_to_xlsx.py`).
- On the PATH: `muscle5` (5.3) and `iqtree2` (2.0.7) for the trees (01_0), BLAST+ (2.12.0) for 06_1, 07_1 and
  07_3, and `wgsim` (1.19.2) and `salmon` (1.10.2 was used) for the paralog simulation (06_2).
- The fonts Arimo and Nimbus Roman. Without them the figures use the system default fonts.
- Internet access the first time a script downloads sequences, genomes or SRA run info from NCBI.

## Layout

```
scripts/          all scripts (run from the repository root or from scripts/)
ncbi_samples/     sample_info.txt (study, sample sheet, location of the quantification) and SRA run info
salmon_quant/     samp_info/ sample sheets, per-study sample tables and the combined samples_all.tsv
msa/              protein sequences, accession lists, alignments and IQ-TREE output per tree
vd3_genes/        the curated gene list (salmon_vd3_orthologs_v1.tsv) and the expression and annotation results
hk_genes/         housekeeping gene candidates (hk_genes.tsv) and their expression
gc_blast/         BLASTp results of the GC search (web BLAST)
d3pathway/        Figure 1, drawn by hand
```

Results tables are included, so each step can be run on its own. Large files are not: gene-level TPM and count
tables (`*.tsv.gz`), genomes and BLAST databases (`annotation/`).

## Workflow

### 1. Read quantification

Reads were downloaded from the SRA (run accessions in Supplementary Table S10 and `ncbi_samples/runinfo/`),
checked with FastQC (0.11.9) and MultiQC (1.17), trimmed with TrimGalore! (0.6.7), quantified against the RefSeq
transcripts of RS_2025_05 with Salmon (1.4.0, `-l A --seqBias --gcBias`) and summarized to genes with tximport
(1.40.0), which gives one `txi_gene.rds` per study.

### 2. Export of the gene tables

`05_1_export_studies_server.R` runs where the `txi_gene.rds` files are. It needs base R only. The `location`
column of `ncbi_samples/sample_info.txt` names each study's directory under `VD3_MAPPED_DIR`:

```
VD3_MAPPED_DIR=/path/to/quant Rscript scripts/05_1_export_studies_server.R
```

It writes `salmon_quant/<study>/{samples.tsv,tpm_gene.tsv.gz,counts_gene.tsv.gz}`, which the other scripts read.

### 3. Tables and figures

```
Rscript scripts/run_all.R
```

runs the steps below in order, reports any step that fails and copies the figure PNGs to `manuscript/v6/`
(set `VD3_FIG_DIR` to write them and the Supplementary Tables elsewhere). A figure whose step failed is not copied.

| Script | Output |
|---|---|
| 00_1_export_tpm.R | atlas TPM of all genes, `tpm_vd3_genes.tsv`, `tpm_hk_genes.tsv` |
| 01_0_build_trees.R | downloads, aligns and builds the trees listed in `msa/<tree>/<tree>_accessions.tsv` |
| 01_phylogenetic_tree.R | Figure 2 (`msa/tree2.png`) and Figure S4 (`msa/tree1.png`) |
| 02_housekeeping_heatmap.R | Figure S3 |
| 03_2_salmon_quant_pca.R | Figure S1 (uses 03_1) |
| 04_1_vd3_genes_heatmap.R | Table S8, expression level, SD and tau per gene |
| 04_2_vd3_genes_heatmap_split.R | Figure 3 and Figure S6 |
| 05_2_combine_studies.R | `salmon_quant/samples_all.tsv` and the per-gene tables of all studies |
| 05_3_vd3_genes_studies.R | exploratory per-study plots in `vd3_genes/studies/` |
| 05_4_atlas_heart_liver.R | sample identity check, Figure S2 |
| 05_5_vd3_genes_studies_fig.R | Figure 4 and Figure S5 |
| 05_6_study_table.R | Table S10 (studies and samples with their SRA runs) |
| 05_7_kidney_detection.R | Table S11 (*cyp27b1* in the kidney samples) |
| 06_1_paralog_similarity.R | Table S9 (paralog similarity, with the 06_2 simulation) |
| 07_2_synteny_table.R | Table S6 (synteny of the paralog pairs, from 07_1) |
| 08_1_ligand_contacts.R | Table S7 (ligand-contact residues from PDB structures) |

### 4. Steps not in run_all

These download large files to `annotation/` and take an hour or more.

| Script | Output |
|---|---|
| 06_2_paralog_simulation.R | simulated paralog mixtures quantified with Salmon (`vd3_genes/paralogs/simulation.tsv`) |
| 07_1_annotation_check.R | gene map to RS_2026_09, reciprocal-BLAST paralog search and synteny (`vd3_genes/annotation_check/`) |
| 07_3_gc_genome_search.R | Table S3 (search for *gc* in six fish genomes, about 3.5 GB of downloads) |

## Settings

`scripts/00_config.R` holds the paths, the tissue and study lists, the excluded samples and genes, and the plot
helpers. Environment variables:

- `VD3_ROOT`: repository root (default: the working directory, or its parent when run from `scripts/`).
- `VD3_FIG_DIR`: where figures and Supplementary Tables are written (default `manuscript/v6`).
- `VD3_TXI_DIR`: directory with the atlas `txi_gene.rds`, used by 00_1 only when
  `salmon_quant/PRJNA72713/tpm_gene.tsv.gz` is missing.
- `VD3_MAPPED_DIR`: parent directory of the study directories, for 05_1.

Figure 1 (`d3pathway/vd3_pathwa_colors.svg`) was drawn in draw.io and exported to PNG with Inkscape.

## License

MIT (see `LICENSE`).
