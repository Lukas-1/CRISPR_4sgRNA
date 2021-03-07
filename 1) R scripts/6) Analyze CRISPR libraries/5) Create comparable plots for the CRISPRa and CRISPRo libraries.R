### 6th March 2021 ##





# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Bar charts")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))




# Define functions --------------------------------------------------------

CommonGenesBarplotMat <- function(var_name) {

  stopifnot(all(c("merged_CRISPRko_df", "merged_replaced_CRISPRa_df") %in% ls(envir = globalenv())))

  CRISPRko_var_df <- BarPlot_Sources(merged_CRISPRko_df,
                                     var_name,
                                     filter_top4           = TRUE,
                                     show_sublibraries     = FALSE,
                                     filter_complete_genes = TRUE
                                     )[["plot_df"]]

  CRISPRa_var_df  <- BarPlot_Sources(merged_replaced_CRISPRa_df,
                                     var_name,
                                     filter_top4           = TRUE,
                                     show_sublibraries     = FALSE,
                                     filter_complete_genes = TRUE
                                     )[["plot_df"]]


  shared_entrezs <- intersect(CRISPRko_var_df[["Entrez_ID"]],
                              CRISPRa_var_df[["Entrez_ID"]]
                              )

  CRISPRko_var_list <- BarPlot_Sources(merged_CRISPRko_df[merged_CRISPRko_df[["Entrez_ID"]] %in% shared_entrezs, ],
                                       var_name,
                                       filter_top4           = TRUE,
                                       show_sublibraries     = FALSE,
                                       filter_complete_genes = TRUE
                                       )

  CRISPRa_var_list  <- BarPlot_Sources(merged_replaced_CRISPRa_df[merged_replaced_CRISPRa_df[["Entrez_ID"]] %in% shared_entrezs, ],
                                       var_name,
                                       filter_top4           = TRUE,
                                       show_sublibraries     = FALSE,
                                       filter_complete_genes = TRUE
                                       )

  results_list <- list(
    "CRISPRa_mat"  = CRISPRa_var_list[["counts_mat"]],
    "CRISPRko_mat" = CRISPRko_var_list[["counts_mat"]]
  )
  return(results_list)
}





# Plot sgRNAs affected by SNPs using shared genes (CRISPRo & a) -----------

shared_barplot_vars <- c(
  "all22_SNP_AF_max_Kaviar",
  "Expected_all22_SNP_AF_max_Kaviar"
)

bar_mat_list_list <- sapply(shared_barplot_vars, CommonGenesBarplotMat, simplify = FALSE)


ManuscriptBars(bar_mat_list_list[["all22_SNP_AF_max_Kaviar"]][["CRISPRa_mat"]])
ManuscriptBars(bar_mat_list_list[["all22_SNP_AF_max_Kaviar"]][["CRISPRko_mat"]])

ManuscriptBars(bar_mat_list_list[["Expected_all22_SNP_AF_max_Kaviar"]][["CRISPRa_mat"]],
               expected_SNP_percent = FALSE
               )
ManuscriptBars(bar_mat_list_list[["Expected_all22_SNP_AF_max_Kaviar"]][["CRISPRko_mat"]],
               expected_SNP_percent = FALSE
               )


labels_list <- list(
  "all22_SNP_AF_max_Kaviar"          = "Target polymorphism > 0.1%",
  "Expected_all22_SNP_AF_max_Kaviar" = "Expected to hit alternate allele"
)


for (SNP_as_percent in c(TRUE, FALSE)) {
  for (var_name in names(bar_mat_list_list)) {
    for (modality in c("CRISPRa", "CRISPRko")) {
      file_name <- paste0(var_name, " - ", modality)
      if (var_name == "Expected_all22_SNP_AF_max_Kaviar") {
        if (SNP_as_percent) {
          file_name <- paste0(file_name, " - percentage")
          use_numeric_limits <- c(0, 2)
        } else {
          file_name <- paste0(file_name, " - absolute number")
          use_numeric_limits <- c(0, 1200)
        }
      } else if (SNP_as_percent) {
        next
      } else {
        use_numeric_limits <- NULL
      }
      pdf(file.path(file_output_directory,
                    paste0("Comparison - ", file_name, ".pdf")
                    ),
          width  = horizontal_width,
          height = horizontal_height
          )
      ManuscriptBars(bar_mat_list_list[[var_name]][[paste0(modality, "_mat")]],
                     axis_label           = labels_list[[var_name]],
                     numeric_limits       = use_numeric_limits,
                     expected_SNP_percent = SNP_as_percent
                     )
      dev.off()
    }
  }
}








