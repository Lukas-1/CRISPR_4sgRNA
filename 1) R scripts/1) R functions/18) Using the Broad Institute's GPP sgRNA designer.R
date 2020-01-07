### 18th December 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))            # For EntrezIDsToSymbols (for the status message)
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))         # For FormatFixedWidthInteger




# Functions for producing input files for the GPP sgRNA designer ----------

FindProblematicEntrezs <- function(CRISPR_df, overview_df) {

  are_top_four <- CRISPR_df[, "Rank"] %in% 1:4

  have_overlaps <- !(CRISPR_df[, "Num_overlaps"] %in% 0) | is.na(CRISPR_df[, "Num_overlaps"])
  meet_criteria <- MeetCriteria(CRISPR_df)

  are_problematic_sgRNAs <- are_top_four & (have_overlaps | !(meet_criteria))
  submit_entrezs <- unlist(strsplit(CRISPR_df[are_problematic_sgRNAs, "Entrez_ID"], ", ", fixed = TRUE))

  are_problematic_genes <- !(overview_df[, "Spacing"] %in% paste0(seq(4, 100, by = 4), "*50"))
  submit_entrezs_genes <- overview_df[are_problematic_genes, "Entrez_ID"]

  already_GPP_entrezs <- unique(CRISPR_df[grepl("GPP", CRISPR_df[, "Source"], fixed = TRUE), "Entrez_ID"])

  submit_entrezs <- unique(c(submit_entrezs, submit_entrezs_genes, already_GPP_entrezs))
  submit_entrezs <- submit_entrezs[!(is.na(submit_entrezs))]

  # Re-order the Entrez IDs
  submit_entrezs <- submit_entrezs[order(as.integer(submit_entrezs))]

  num_all <- length(submit_entrezs)
  num_previous <- length(already_GPP_entrezs)
  num_new <- num_all - num_previous

  show_message <- paste0(num_all, " genes (Entrez IDs) ought to be submitted to the GPP sgRNA designer.")
  if (num_new == 0) {
    show_message <- paste0(show_message, " SgRNAs from GPP were previously available for all of them.")
  } else {
    show_message <- paste0(show_message, " Of these, previous data from GPP were not found for ", num_new, " of them")
    if (num_new < 10) {
      new_entrezs <- submit_entrezs[!(submit_entrezs %in% already_GPP_entrezs)]
      symbols_vec <- EntrezIDsToSymbols(new_entrezs)
      show_message <- paste0(show_message, ": ", paste0(vapply(seq_along(new_entrezs), function(x) paste0(new_entrezs[[x]], " (", symbols_vec[[x]], ")"), ""), collapse = ", "))
    } else {
      show_message <- paste0(show_message, ".")
    }

  }
  message(show_message)

  return(submit_entrezs)
}



BuildDfForGPP <- function(submit_entrezs) {

  num_genes_per_file <- 200L
  num_genes <- length(submit_entrezs)
  num_files <- ceiling(num_genes / num_genes_per_file)

  file_sequence <- rep(seq_len(num_files), each = num_genes_per_file)
  file_sequence <- file_sequence[seq_len(num_genes)]

  submit_df <- data.frame("Entrez_ID"      = submit_entrezs,
                          "File_number"    = file_sequence,
                          "File_name"      = FormatFixedWidthInteger(file_sequence),
                          stringsAsFactors = FALSE,
                          row.names        = NULL
                          )
  return(submit_df)
}




WriteGPPInputDf <- function(submit_df, chunk_ID, GPP_input_directory, input_prefix = "") {
  num_files <- length(unique(submit_df[, "File_number"]))
  for (i in seq_len(num_files)) {
    are_this_file <- submit_df[, "File_number"] == i
    file_name <- paste0("GPPsg_", input_prefix, "input__chunk_", chunk_ID,
                        "__file_", submit_df[are_this_file, "File_name"][[1]],
                        ".txt"
                        )
    write.table(submit_df[are_this_file, "Entrez_ID"],
                file = file.path(GPP_input_directory, file_name),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )
  }
  return(invisible(NULL))
}







# Functions for processing output from the GPP sgRNA designer -------------

# ReadGPPOutputFiles <- function(output_file_names, GPP_path) {
#   GPP_output_df_list <- sapply(output_file_names, function (x) {
#     GPP_output_colnames <- scan(file = file.path(GPP_path, x), nlines = 1, what = "character", sep = "\t")
#     GPP_output_df <- read.table(file.path(GPP_path, x),
#                                 sep = "\t", quote = "", header = FALSE,
#                                 fill = TRUE, skip = 5, stringsAsFactors = FALSE
#                                 )
#     colnames(GPP_output_df) <- GPP_output_colnames
#     return(GPP_output_df)
#   }, simplify = FALSE)
#   results_df <- do.call(rbind.data.frame, c(GPP_output_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
#   return(results_df)
# }

ReadGPPOutputFiles <- function(output_file_names, GPP_path) {
  GPP_output_df_list <- sapply(output_file_names, function (x) {
    GPP_output_df <- read.table(file.path(GPP_path, x),
                                sep = "\t", quote = "", header = TRUE,
                                fill = TRUE, stringsAsFactors = FALSE
                                )
    return(GPP_output_df)
  }, simplify = FALSE)
  results_df <- do.call(rbind.data.frame, c(GPP_output_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}



CRISPRa_GPP_output_columns <- c(
  "Input", "Target Gene ID", "Target Gene Symbol",
  "Reference Sequence", "Strand of Target", "TSS Position",
  "Strand of sgRNA", "sgRNA Cut Position",
  "sgRNA Sequence", "sgRNA Context Sequence", "PAM Sequence",
  "sgRNA Cut Site TSS Offset",
  "# Off-Target Tier I Match Bin I Matches",   "# Off-Target Tier II Match Bin I Matches",   "# Off-Target Tier III Match Bin I Matches",
  "# Off-Target Tier I Match Bin II Matches",  "# Off-Target Tier II Match Bin II Matches",  "# Off-Target Tier III Match Bin II Matches",
  "# Off-Target Tier I Match Bin III Matches", "# Off-Target Tier II Match Bin III Matches", "# Off-Target Tier III Match Bin III Matches",
  "# Off-Target Tier I Match Bin IV Matches",  "# Off-Target Tier II Match Bin IV Matches",  "# Off-Target Tier III Match Bin IV Matches",
  "On-Target Efficacy Score", "DHS Score",
  "On-Target Rank", "Off-Target Rank", "Combined Rank",
  "Pick Order", "Picking Round", "Picking Notes"
)

CRISPRko_GPP_output_columns <- c(
  "Input", "Target Gene ID", "Target Gene Symbol",
  "Target Transcript",
  "Reference Sequence", "Strand of Target",
  "Strand of sgRNA", "Orientation", "sgRNA Cut Position (1-based)",
  "sgRNA Sequence", "sgRNA Context Sequence", "PAM Sequence",
  "Exon Number",
  "Target Cut Length", "Target Total Length", "Target Cut %",
  "# Off-Target Tier I Match Bin I Matches",   "# Off-Target Tier II Match Bin I Matches",   "# Off-Target Tier III Match Bin I Matches",   "# Off-Target Tier IV Match Bin I Matches",
  "# Off-Target Tier I Match Bin II Matches",  "# Off-Target Tier II Match Bin II Matches",  "# Off-Target Tier III Match Bin II Matches",  "# Off-Target Tier IV Match Bin II Matches",
  "# Off-Target Tier I Match Bin III Matches", "# Off-Target Tier II Match Bin III Matches", "# Off-Target Tier III Match Bin III Matches", "# Off-Target Tier IV Match Bin III Matches",
  "# Off-Target Tier I Match Bin IV Matches",  "# Off-Target Tier II Match Bin IV Matches",  "# Off-Target Tier III Match Bin IV Matches",  "# Off-Target Tier IV Match Bin IV Matches",
  "On-Target Efficacy Score",
  "On-Target Rank", "Off-Target Rank", "Combined Rank",
  "Pick Order", "Picking Round", "Picking Notes"
)





TidyGPPOutputDf <- function(GPP_output_df, choose_columns) {

  were_not_found <- grepl("^ERROR: Gene .+ not found", GPP_output_df[, "Picking Notes"]) & is.na(GPP_output_df[, "Quota"])

  message(paste0(sum(were_not_found), " genes were not found by the Broad Institute's Genetic Perturbation Platform (GPP) sgRNA picker tool and were omitted from the data frame!"))

  results_df <- GPP_output_df[!(were_not_found), choose_columns]
  rownames(results_df) <- NULL
  return(results_df)
}












