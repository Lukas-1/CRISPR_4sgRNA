### 18th July 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
source(file.path(general_functions_directory, "05) Mapping sequences to the human genome.R"))
source(file.path(general_functions_directory, "07) Annotating mapped sequences with additional information.R"))
source(file.path(general_functions_directory, "31) Finding 0MM hits for variable-length sgRNAs.R"))



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")




# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "08) Categorize subsequences of reads aligned to the reference.RData"))
load(file.path(sql2_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))




# Define functions --------------------------------------------------------

FilterMutations <- function(extract_df, use_zmws, min_length = 19, max_length = 42) {

  are_sg <- extract_df[["Feature"]] %in% paste0("sg", 1:4)
  are_mutated <- extract_df[["Category"]] %in% "Mutation"
  are_eligible <- extract_df[["ZMW"]] %in% use_zmws

  mutations_df <- extract_df[are_sg & are_mutated & are_eligible, ]

  mutations_df[["Read_without_gaps"]] <- gsub("-", "", mutations_df[["Aligned_read"]], fixed = TRUE)

  sequence_lengths <- nchar(mutations_df[["Read_without_gaps"]])
  are_too_short <- sequence_lengths < min_length
  are_too_long <- sequence_lengths > max_length
  mutations_df <- mutations_df[!(are_too_short | are_too_long), ]
  row.names(mutations_df) <- NULL
  return(mutations_df)
}




# Filter reads ------------------------------------------------------------

pass_filters <- ccs7_df_list[["individual_reads_df"]][["Passes_filters"]] == 1
ccs7_filtered_zmws <- ccs7_df_list[["individual_reads_df"]][["ZMW"]][pass_filters]




# Define all "mutated" sequences to be considered -------------------------

mutations_df <- FilterMutations(extracted_df, ccs7_filtered_zmws)




# Annotate mutated sequences with 0MM off-target sites --------------------

groooo
mutations_df <- Add0MMHits(mutations_df, "Read_without_gaps")




# Save data ---------------------------------------------------------------

save(list = "mutations_df",
     file = file.path(sql2_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData")
     )



