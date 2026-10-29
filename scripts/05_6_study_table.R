# Table S10: the RNA-seq studies, with tissues and sample counts as analyzed (atlas and replicate studies),
# and one row per sample with its SRA runs, layout, instrument and read length.
# Counts come from samples_all.tsv (05_2). Excluded samples (00_config.R) are listed but not counted, and
# relabelled tissues (tissue_raw differs from tissue) are noted. Design and references are written by hand
# below and follow the Methods. SRA run info is downloaded once per BioProject to ncbi_samples/runinfo/.
# A sample is matched to its runs by the run accession in its sample sheet, or else by the SRA library name.
# Needs python3 for the Excel file, and internet access on the first run.
# Input:  salmon_quant/samples_all.tsv, salmon_quant/<study>/samples.tsv
# Output: vd3_genes/studies/{study_table,sample_runs}.tsv, manuscript/v7/Table S10 ... .xlsx

source(if (file.exists("00_config.R")) "00_config.R" else "scripts/00_config.R")

out_file <- file.path(vd3_dir, "studies", "study_table.tsv")
table_file <- file.path(fig_dir, "Table S10 RNA-seq studies and samples.xlsx")

design <- tribble(
  ~study, ~bioproject, ~reference, ~design,
  "PRJNA72713", "PRJNA72713", "Lien et al. 2016", "Multi-tissue atlas, one sample per tissue",
  "PRJNA359635", "PRJNA359635", "Robertson et al. 2017", "Females and males in freshwater and seawater, one pool of three fish per group",
  "PRJEB47409", "PRJEB47409", "Nguyen et al. 2025", "Immature and mature females and males, three fish per group",
  "ROS_1", "PRJNA765504, PRJNA765887, PRJNA766174, PRJNA777794", "Yin et al. 2022",
  "18 fish fed a low-antioxidant diet (Low and LowSham groups)",
  "PRJNA647285", "PRJNA647285", "Gervais et al. 2021", "Uninfected parr from resistant and susceptible families, before a virus challenge",
  "PRJNA543940", "PRJNA543940", "Robinson et al. 2020", "Parr from families with high and low genomic breeding values for pancreas disease resistance",
  "PRJNA380580", "PRJNA380580", "Skaftnesmo et al. 2017", "Males at three maturation stages, three fish per stage",
  "PRJNA626461", "PRJNA626461", "Andrew et al. 2021", "Immature parr homozygous for the early (EE) or late (LL) vgll3 allele"
)
stopifnot(setequal(design$study, c("PRJNA72713", replicate_studies)))

samples <- read_tsv(file.path(proj_dir, "salmon_quant", "samples_all.tsv"), show_col_types = FALSE) %>%
  filter(study %in% design$study)

tissues <- samples %>%
  filter(!excluded) %>%
  count(study, tissue) %>%
  group_by(study) %>%
  summarise(tissues = paste0(tissue_label(tissue), " (", n, ")", collapse = ", "),
            n_samples = sum(n), .groups = "drop")

# Spelling differences (Head Kidney, Pyloric caecum) and gonads split by sex are not relabels
relabelled <- function(raw, tissue) {
  str_replace_all(str_to_lower(raw), " ", "_") != tissue & str_to_lower(raw) != "gonad"
}
notes <- samples %>%
  mutate(note = case_when(
    excluded ~ paste0("Sample ", sample, " excluded (", excluded_samples[id], ")"),
    relabelled(tissue_raw, tissue) ~ paste0(tissue_raw, " samples analyzed as ", str_to_lower(tissue_label(tissue))))) %>%
  filter(!is.na(note)) %>%
  distinct(study, note) %>%
  group_by(study) %>%
  summarise(notes = paste(note, collapse = ". "), .groups = "drop")

study_table <- design %>%
  left_join(tissues, by = "study") %>%
  left_join(notes, by = "study") %>%
  transmute(Study = if_else(study == "ROS_1", "ROS_1 (ROS Work Package 1)", study),
            BioProject = bioproject, Reference = reference, Design = design,
            `Tissues (samples)` = tissues, `Samples analyzed` = n_samples, Notes = replace_na(notes, ""))

write_tsv(study_table, out_file)

# SRA runs of each analyzed sample
runinfo_dir <- file.path(proj_dir, "ncbi_samples", "runinfo")
dir.create(runinfo_dir, showWarnings = FALSE)
runinfo <- map_dfr(unlist(str_split(design$bioproject, ", ")), function(p) {
  f <- file.path(runinfo_dir, paste0(p, ".csv"))
  if (!file.exists(f)) download.file(paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=", p),
                                     f, quiet = TRUE)
  read_csv(f, col_types = cols(.default = col_character()))
}) %>%
  filter(!is.na(Run), Run != "Run") %>%
  transmute(run = Run, library = str_remove(LibraryName, "_L00\\d$"), layout = LibraryLayout, instrument = Model,
            spots = as.numeric(spots), bases = as.numeric(bases))

sheets <- map_dfr(design$study, function(s) {
  read_tsv(file.path(proj_dir, "salmon_quant", s, "samples.tsv"), col_types = cols(.default = col_character())) %>%
    transmute(study, sample, key = coalesce(!!!syms(intersect(c("read_group", "read_file"), names(.)))),
              lims = if ("lims" %in% names(.)) lims else NA_character_)
})
by_run <- inner_join(sheets, runinfo, by = c(key = "run"), keep = TRUE)
by_library <- sheets %>%
  anti_join(by_run, by = c("study", "sample")) %>%
  mutate(library = coalesce(na_if(lims, ""), key)) %>%
  inner_join(runinfo, by = "library")
sample_runs <- bind_rows(by_run, by_library) %>%
  group_by(study, sample) %>%
  summarise(runs = paste(sort(run), collapse = ", "), layout = paste(unique(layout), collapse = ", "),
            instrument = paste(unique(instrument), collapse = ", "), m_spots = sum(spots) / 1e6,
            read_length = sum(bases) / sum(spots) / if_else(dplyr::first(layout) == "PAIRED", 2, 1), .groups = "drop")
analyzed <- samples %>% filter(!excluded)
missing <- anti_join(analyzed, sample_runs, by = c("study", "sample"))
if (nrow(missing)) stop("No SRA run for ", paste(missing$id, collapse = ", "))

# Group: stage and treatment, without "normal" and the read-group labels that some sheets keep as treatment
sample_table <- analyzed %>%
  inner_join(sample_runs, by = c("study", "sample")) %>%
  mutate(treatment = if_else(treatment %in% c("normal", stage) | str_detect(treatment, "_"), NA, treatment),
         group = str_squish(paste(replace_na(stage, ""), replace_na(treatment, "")))) %>%
  mutate(study = factor(study, levels = design$study)) %>%
  arrange(study, tissue, sample) %>%
  transmute(Study = as.character(study), Sample = sample, Tissue = tissue_label(tissue),
            Sex = replace_na(na_if(sex, "NA"), ""), Group = group,
            `SRA runs` = runs, Layout = str_to_lower(layout), Instrument = instrument,
            `Read length (nt)` = round(read_length), `Spots (M)` = round(m_spots, 1),
            `Salmon reads (M)` = m_counts)
write_tsv(sample_table, file.path(vd3_dir, "studies", "sample_runs.tsv"))

tmp <- tempfile()
dir.create(tmp)
write_tsv(study_table, file.path(tmp, "Studies.tsv"))
write_tsv(sample_table, file.path(tmp, "Samples.tsv"), na = "")
cat("\nRead length is the mean length of one read (one mate for paired runs), from the SRA run info. ",
    "Salmon reads are the estimated counts summed over all genes.\n", file = file.path(tmp, "Samples.tsv"),
    append = TRUE, sep = "")
py <- file.path(proj_dir, "scripts", "tsv_to_xlsx.py")
if (system2("python3", c(shQuote(py), shQuote(file.path(tmp, c("Studies.tsv", "Samples.tsv"))), shQuote(table_file))) != 0)
  stop("tsv_to_xlsx.py failed")
message("Wrote ", out_file, " and ", table_file)
