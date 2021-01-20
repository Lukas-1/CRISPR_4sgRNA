### 20th November 2019 ###




# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRko_RData_directory  <- file.path(RData_directory, "3) CRISPRko")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRko", "GPP sgRNA designer", "1) Input files")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))
load(file.path(CRISPRko_RData_directory, "12) Create a gene-based summary of the human genome - sgRNAs_overview_df.RData"))





# Collect Entrez IDs for submission to the GPP sgRNA designer -------------

problematic_entrezs <- FindProblematicEntrezs(merged_CRISPRko_df, sgRNAs_overview_df)

top_4_df <- merged_CRISPRko_df[(merged_CRISPRko_df[["Rank"]] %in% 1:4), ]
table(top_4_df[(top_4_df[["Source"]] %in% "GPP"), "Entrez_ID"] %in% problematic_entrezs)





# Separate the Entrez IDs into chunks -------------------------------------

submit_entrez_chunks_list   <- lapply(entrez_chunks_list, function(x) intersect(x, problematic_entrezs))
optional_entrez_chunks_list <- lapply(entrez_chunks_list, function(x) setdiff(x, problematic_entrezs))





# Build data frames for submission to the GPP sgRNA designer --------------

submit_df_list <- lapply(submit_entrez_chunks_list, BuildDfForGPP)
optional_df_list <- lapply(optional_entrez_chunks_list, BuildDfForGPP)

combined_submit_df_list <- CombineDfChunks(submit_df_list, max_num_per_chunk = 200L)





# Write input files for the GPP sgRNA designer to disk --------------------

for (chunk_ID in names(submit_df_list)) {
  WriteGPPInputDf(submit_df_list[[chunk_ID]],
                  chunk_ID,
                  file.path(GPP_input_files_directory, "1) High-priority")
                  )
}

for (chunk_ID in names(optional_df_list)) {
  WriteGPPInputDf(optional_df_list[[chunk_ID]],
                  chunk_ID,
                  file.path(GPP_input_files_directory, "2) Optional"),
                  input_prefix = "optional_"
                  )
}

for (chunk_ID in names(combined_submit_df_list)) {
  WriteGPPInputDf(combined_submit_df_list[[chunk_ID]],
                  sub("chunk_", "", chunk_ID),
                  file.path(GPP_input_files_directory, "3) Combined"),
                  input_prefix = "combined_"
                  )
}





# Save data ---------------------------------------------------------------

save(list = "problematic_entrezs",
     file = file.path(CRISPRko_RData_directory, "18) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData")
     )




