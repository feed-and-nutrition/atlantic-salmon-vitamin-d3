# Rebuild all tables and figures, then copy the figure PNGs into manuscript/v6 (fig_dir in 00_config.R).
# Run from the project root or scripts/:  Rscript scripts/run_all.R
# A step that fails (for example 00_1 without the 15-tissue txi) is reported and skipped.

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

steps <- c("00_1_export_tpm.R", "01_0_build_trees.R", "01_phylogenetic_tree.R", "02_housekeeping_heatmap.R",
           "03_2_salmon_quant_pca.R", "04_1_vd3_genes_heatmap.R", "04_2_vd3_genes_heatmap_split.R",
           "05_2_combine_studies.R", "05_3_vd3_genes_studies.R", "05_4_atlas_heart_liver.R",
           "05_5_vd3_genes_studies_fig.R", "05_6_study_table.R", "05_7_kidney_detection.R", "06_1_paralog_similarity.R",
           "07_2_synteny_table.R", "08_1_ligand_contacts.R")

old_wd <- setwd(file.path(proj_dir, "scripts"))
status <- map_chr(steps, function(f) {
  message("== ", f)
  tryCatch({
    source(f, local = new.env())
    "ok"
  }, error = function(e) {
    message("   FAILED: ", conditionMessage(e))
    paste("failed:", conditionMessage(e))
  })
})
setwd(old_wd)

# Figure number -> source PNG (Fig 1 is drawn by hand in d3pathway/vd3_pathwa_colors.svg, see CLAUDE.md)
figures <- c(fig1 = "d3pathway/pathway.png",
             fig2 = "msa/tree2.png",
             fig3 = "vd3_genes/vd3_genes_all_expr.png",
             fig4 = "vd3_genes/studies/vd3_genes_tissues.png",
             fig_s1 = "salmon_quant/PRJNA72713/salmon_stats.png",
             fig_s2 = "salmon_quant/PRJNA72713/atlas_heart_liver.png",
             fig_s3 = "hk_genes/hk_genes.png",
             fig_s4 = "msa/tree1.png",
             fig_s5 = "vd3_genes/studies/vd3_genes_per_study.png",
             fig_s6 = "vd3_genes/vd3_genes_all_expr_z.png")
made_by <- c(fig1 = NA, fig2 = "01_phylogenetic_tree.R",
             fig3 = "04_2_vd3_genes_heatmap_split.R", fig4 = "05_5_vd3_genes_studies_fig.R",
             fig_s1 = "03_2_salmon_quant_pca.R", fig_s2 = "05_4_atlas_heart_liver.R",
             fig_s3 = "02_housekeeping_heatmap.R", fig_s4 = "01_phylogenetic_tree.R",
             fig_s5 = "05_5_vd3_genes_studies_fig.R", fig_s6 = "04_2_vd3_genes_heatmap_split.R")
# Copy only figures whose step succeeded, so a stale PNG never replaces a manuscript figure
ok_steps <- steps[status == "ok"]
to_copy <- names(figures)[is.na(made_by) | made_by %in% ok_steps]
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
copied <- set_names(map_lgl(to_copy, ~ file.copy(file.path(proj_dir, figures[[.x]]),
                                                  file.path(fig_dir, paste0(.x, ".png")), overwrite = TRUE)), to_copy)

message("\nSteps:")
walk2(steps, status, ~ message(sprintf("  %-34s %s", .x, .y)))
message("Figures copied to ", fig_dir, ": ", paste(names(copied)[copied], collapse = ", "))
skipped <- setdiff(names(figures), names(copied)[copied])
if (length(skipped)) message("Figures NOT updated: ", paste(skipped, collapse = ", "))
