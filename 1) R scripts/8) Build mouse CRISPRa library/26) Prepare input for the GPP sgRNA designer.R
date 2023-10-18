### 7th October 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR_4sgRNA"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "6) Mouse - General")
CRISPRa_RData_directory   <- file.path(RData_directory, "7) Mouse - CRISPRa")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "Mouse - CRISPRa", "GPP sgRNA designer", "1) Input files")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "04) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRa_RData_directory, "20) Create a gene-based summary of the mouse genome - sgRNAs_overview_df.RData"))





# Identify problematic Entrez IDs -----------------------------------------

problematic_entrezs <- FindProblematicEntrezs(merged_replaced_CRISPRa_df, sgRNAs_overview_df)






# Build data frames for submission to the GPP sgRNA designer --------------

all_genes_df_list <- lapply(entrez_chunks_list, BuildDfForGPP)







# Write input files for the GPP sgRNA designer to disk --------------------

for (chunk_ID in names(all_genes_df_list)) {
  WriteGPPInputDf(all_genes_df_list[[chunk_ID]],
                  chunk_ID,
                  file.path(GPP_input_files_directory, "All genes")
                  )
}





# Save data ---------------------------------------------------------------

save(list = "problematic_entrezs",
     file = file.path(CRISPRa_RData_directory, "26) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData")
     )






