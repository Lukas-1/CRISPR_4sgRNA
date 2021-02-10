### 9th April 2020 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "18) Using the Broad Institute's GPP sgRNA designer.R"))
source(file.path(general_functions_directory, "21) Splitting sgRNAs into chunks for parallel analysis.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory     <- "~/CRISPR"
RData_directory           <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory   <- file.path(RData_directory, "1) General")
CRISPRi_RData_directory   <- file.path(RData_directory, "4) CRISPRi")
GPP_input_files_directory <- file.path(CRISPR_root_directory, "4) Intermediate files", "CRISPRi", "GPP sgRNA designer", "1) Input files")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "09) Divide the entire set of protein-coding genes into chunks - entrez_chunks_list.RData"))
load(file.path(CRISPRi_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRi_RData_directory, "20) Create a gene-based summary of the human genome.RData"))
load(file.path(CRISPRi_RData_directory, "20) Create a gene-based summary of the human genome - vacuolation_entrezs.RData"))





# Identify problematic Entrez IDs -----------------------------------------

problematic_entrezs <- FindProblematicEntrezs(merged_replaced_CRISPRi_df, sgRNAs_overview_df)

vacuolation_only_entrezs <- setdiff(vacuolation_entrezs, sgRNAs_overview_df[["Entrez_ID"]])

problematic_entrezs <- c(vacuolation_only_entrezs, problematic_entrezs)





# Build data frames for submission to the GPP sgRNA designer --------------

all_genes_df_list <- lapply(entrez_chunks_list, BuildDfForGPP)

vacuolation_df <- BuildDfForGPP(vacuolation_entrezs)
vacuolation_only_df <- BuildDfForGPP(setdiff(vacuolation_entrezs, unlist(entrez_chunks_list, use.names = FALSE)))

all_genes_df_list <- c(all_genes_df_list, list("Vacuolation_only" = vacuolation_only_df))





# Write input files for the GPP sgRNA designer to disk --------------------

for (chunk_ID in names(all_genes_df_list)) {
  WriteGPPInputDf(all_genes_df_list[[chunk_ID]],
                  chunk_ID,
                  file.path(GPP_input_files_directory, "All genes")
                  )
}

WriteGPPInputDf(vacuolation_df,
                "vacuolation",
                file.path(GPP_input_files_directory, "Vacuolation genes"),
                input_prefix = ""
                )




# Save data ---------------------------------------------------------------

save(list = "problematic_entrezs",
     file = file.path(CRISPRi_RData_directory, "26) Prepare input for the GPP sgRNA designer - problematic_entrezs.RData")
     )




