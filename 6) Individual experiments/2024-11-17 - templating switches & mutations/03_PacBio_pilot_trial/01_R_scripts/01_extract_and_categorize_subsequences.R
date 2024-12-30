## 2024-11-17


# Load packages and source code -------------------------------------------

root_dir <- "~/CRISPR_4sgRNA"
exper_dir <- file.path(root_dir, "6) Individual experiments")

first_QC_dir <- file.path(exper_dir, "2020-08-29 - PacBio - first 384-well plate")
nanopore_dir <- file.path(exper_dir, "2022-01-05 - first nanopore sequencing run")
project_dir  <- file.path(exper_dir, "2024-11-17 - templating switches & mutations")

source(file.path(first_QC_dir, "1) R functions", "07) Categorizing subsequences of reads aligned to the reference.R"))
source(file.path(nanopore_dir, "01_R_scripts", "1_R_functions", "03_extracting_aligned_sgRNAs.R"))
source(file.path(project_dir, "01_R_functions", "01_extracting_and_categorizing_subsequences.R"))



# Define paths ------------------------------------------------------------

rdata_dir <- file.path(project_dir, "03_PacBio_pilot_trial", "02_R_objects")
first_rdata_dir <- file.path(exper_dir, "2022-04-06 - PacBio pooled 4sg - first trial", "03_R_objects")



# Read in data ------------------------------------------------------------

amplicon_ref <- read.table(file.path(nanopore_dir, "02_input_data", "amplicon_4sg.txt"),
                           quote = "", stringsAsFactors = FALSE
                           )[, 1]


# Define maps -------------------------------------------------------------

features_list <- list(
  "promoter1_hU6"  = c(177, 426),
  "sg1"            = c(427, 446),
  "sg1_cr1"        = c(427, 532),
  "tracrRNA1"      = c(447, 532),
  "polyT_1"        = c(533, 539),

  "EM7_promoter"   = c(540, 587),
  "pre_TpR"        = c(588, 605),
  "TpR_DHFR"       = c(606, 842),
  "polyT_TpR"      = c(843, 849),

  "promoter2_mU6"  = c(850, 1165),
  "sg2"            = c(1166, 1185),
  "sg2_cr2"        = c(1166, 1273),
  "tracrRNA2"      = c(1186, 1273),
  "polyT_2"        = c(1274, 1280),

  "promoter3_hH1"  = c(1281, 1504),
  "sg3"            = c(1505, 1524),
  "sg3_cr3"        = c(1505, 1612),
  "tracrRNA3"      = c(1525, 1612),
  "polyT_3"        = c(1613, 1619),

  "promoter4_h7SK" = c(1620, 1863),
  "sg4"            = c(1864, 1883),
  "sg4_cr4"        = c(1864, 1969),
  "tracrRNA4"      = c(1884, 1969),
  "polyT_4"        = c(1970, 1976)
)



# Load data ---------------------------------------------------------------

load(file.path(first_rdata_dir, "02_align_reads.RData"))
load(file.path(first_rdata_dir, "07_assign_sgRNAs_to_plasmids.RData"))



# Prepare features_df -----------------------------------------------------

features_df <- TweakFeaturesDf(FeaturesListToDf(features_list))
features_indices_list <- lapply(seq_len(nrow(features_df)),
                                function(y) features_df[y, "Start"]:features_df[y, "End"]
                                )


# Extract aligned sequences -----------------------------------------------

extracted_df <- ExtractAllFeatures(FilterAlignmentsDf(alignments_df, pb_df))



# Categorize extracted sequences ------------------------------------------

categorized_df <- AddThreeCategories(extracted_df, features_df)



# Incorporate data on sgRNAs ----------------------------------------------

categorized_df <- AddGuideData(categorized_df, pb_df)



# Save data ---------------------------------------------------------------

save(categorized_df, file = file.path(rdata_dir, "01_extract_and_categorize_sequences__categorized_df.RData"))
save(features_df, file = file.path(rdata_dir, "01_extract_and_categorize_sequences__features_df.RData"))


