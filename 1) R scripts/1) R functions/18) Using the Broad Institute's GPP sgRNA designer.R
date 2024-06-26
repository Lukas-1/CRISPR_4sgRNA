### 18th December 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR_4sgRNA/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "01) Retrieving annotation data for a gene.R"))            # For EntrezIDsToSymbols (for the status message)
source(file.path(general_functions_directory, "16) Producing per-gene summaries of CRISPR libraries.R")) # For MeetCriteria
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))         # For FormatFixedWidthInteger




# Functions for producing input files for the GPP sgRNA designer ----------

PrettyGeneString <- function(entrezs_vec) {
  symbols_vec <- EntrezIDsToSymbols(new_entrezs)
  result_vec <- mapply(function(x, y) paste0(x, " (", y, ")"), entrezs_vec, symbols_vec)
  result_string <- paste0(result_vec, collapse = ", ")
  return(result_string)
}



FindProblematicEntrezs <- function(CRISPR_df, overview_df) {
  are_top_four <- CRISPR_df[["Rank"]] %in% 1:4

  have_overlaps <- !(CRISPR_df[["Num_overlaps"]] %in% 0) | is.na(CRISPR_df[["Num_overlaps"]])
  meet_criteria <- MeetCriteria(CRISPR_df)
  have_multiple_locations <- (CRISPR_df[["Num_0MM"]] == 1) %in% FALSE

  are_problematic_sgRNAs <- are_top_four & (have_overlaps | (!(meet_criteria)) | have_multiple_locations)
  submit_entrezs <- unlist(strsplit(CRISPR_df[["Entrez_ID"]][are_problematic_sgRNAs], ", ", fixed = TRUE))

  are_problematic_genes <- !(overview_df[["Spacing"]] %in% paste0(seq(4, 100, by = 4), "*50"))
  submit_entrezs_genes <- overview_df[["Entrez_ID"]][are_problematic_genes]

  are_GPP <- CRISPR_df[["Source"]] %in% c("GPP", "Curated, GPP")
  are_top4 <- CRISPR_df[["Rank"]] %in% 1:4
  top4_GPP_entrezs <- unique(CRISPR_df[["Entrez_ID"]][are_GPP & are_top4])
  have_GPP_entrezs <- unique(CRISPR_df[["Entrez_ID"]][are_GPP])

  required_entrezs <- setdiff(c(submit_entrezs, submit_entrezs_genes, top4_GPP_entrezs),
                              NA_character_
                              )

  # Re-order the Entrez IDs
  required_entrezs <- required_entrezs[order(as.integer(required_entrezs))]
  not_found_entrezs <- setdiff(required_entrezs, have_GPP_entrezs)
  num_not_found <- length(not_found_entrezs)

  show_message <- paste0(length(required_entrezs), " genes (Entrez IDs) ought to be submitted to the GPP sgRNA designer.")
  if (num_not_found == 0) {
    show_message <- paste0(show_message, " SgRNAs from GPP were previously available for all of them.")
  } else {
    show_message <- paste0(show_message, " Of these, previous data from GPP were not found for ",
                           num_not_found, " of them"
                           )
    if (num_not_found < 10) {
      show_message <- paste0(show_message, ": ", PrettyGeneString(not_found_entrezs))
    } else {
      show_message <- paste0(show_message, ".")
    }
  }
  message(show_message)
  return(required_entrezs)
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
  num_files <- length(unique(submit_df[["File_number"]]))
  for (i in seq_len(num_files)) {
    are_this_file <- submit_df[["File_number"]] == i
    file_name <- paste0("GPPsg_", input_prefix, "input__chunk_", chunk_ID,
                        "__file_", submit_df[["File_name"]][are_this_file][[1]],
                        ".txt"
                        )
    write.table(submit_df[["Entrez_ID"]][are_this_file],
                file = file.path(GPP_input_directory, file_name),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )
  }
  return(invisible(NULL))
}







# Functions for processing output from the GPP sgRNA designer -------------

ReadGPPOutputFiles <- function(GPP_folder_paths) {
  file_paths_vec <- unlist(lapply(GPP_folder_paths, function(x) file.path(x, list.files(x))))
  GPP_output_df_list <- sapply(file_paths_vec, function (x) {
    GPP_output_df <- read.table(x, sep = "\t", quote = "", comment.char = "",
                                header = TRUE, check.names = FALSE, fill = TRUE,
                                stringsAsFactors = FALSE
                                )
    return(GPP_output_df)
  }, simplify = FALSE)
  results_df <- do.call(rbind.data.frame, c(GPP_output_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  are_duplicated <- duplicated(results_df)
  if (any(are_duplicated)) {
    message(paste0(sum(are_duplicated), " duplicate entries were omitted from the collected GPP output."))
    results_df <- results_df[!(are_duplicated), ]
    rownames(results_df) <- NULL
  }
  return(results_df)
}


CRISPRa_GPP_output_columns <- c(
  "Input", "Target Gene ID", "Target Gene Symbol",
  "Reference Sequence", "Strand of Target", "TSS Position",
  "Strand of sgRNA", "sgRNA 'Cut' Position",
  "sgRNA Sequence", "sgRNA Context Sequence", "PAM Sequence",
  "sgRNA 'Cut' Site TSS Offset",
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

  are_errors <- grepl("ERROR", GPP_output_df[["Picking Notes"]], fixed = TRUE)

  were_not_found <- (grepl("^ERROR: Gene .+ not found", GPP_output_df[["Picking Notes"]]) |
                     grepl("ERROR: Unable to locate a current transcript for ", GPP_output_df[["Picking Notes"]], fixed = TRUE)
                     ) & is.na(GPP_output_df[["Quota"]])

  message(paste0(sum(were_not_found), " genes were not found by the Broad Institute's ",
                 "Genetic Perturbation Platform (GPP) sgRNA picker tool and were omitted from the data frame!"
                 )
          )
  are_other_error <- are_errors & !(were_not_found)
  if (any(are_other_error)) {
    message("")
    message(paste0(sum(are_other_error), " other error messages were found and are shown below:"))
    print(GPP_output_df[are_other_error, ])
    message("")
  }
  results_df <- GPP_output_df[!(are_errors), choose_columns]
  row.names(results_df) <- NULL
  return(results_df)
}




FilterGPPOutputDf <- function(GPP_output_df, problematic_genes, n_unproblematic, n_problematic = NULL) {
  split_df_list <- split(GPP_output_df, GPP_output_df[["Target Gene ID"]])
  filtered_df_list <- lapply(split_df_list, function(x) {
    is_problematic <- x[["Target Gene ID"]][[1]] %in% problematic_genes
    if (is_problematic) {
      if (is.null(n_problematic)) {
        return(x)
      } else {
        use_n <- n_problematic
      }
    } else {
      use_n <- n_unproblematic
    }
    use_seq <- seq_len(use_n)
    use_indices <- which(x[["Pick Order"]] %in% use_seq)
    if (length(use_indices) < use_n) {
      use_indices <- c(use_indices, which(x[["Combined Rank"]] %in% use_seq))
      if (length(use_indices) > use_n) {
        use_indices <- use_indices[use_seq]
      }
    }
    return(x[use_indices, ])
  })
  results_df <- do.call(rbind.data.frame, c(filtered_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}








