### 3rd May 2021 ###



# Import packages and source code -----------------------------------------

library("Biostrings")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
sql2_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-04-03 - PacBio - first Sequel-II run")
sql2_R_objects_directory <- file.path(sql2_directory, "3) R objects")



# Load data ---------------------------------------------------------------

load(file.path(sql2_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(sql2_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(sql2_R_objects_directory, "05) Read in PacBio data.RData"))
load(file.path(sql2_R_objects_directory, "07) Extract barcode sequences and quality scores.RData"))
load(file.path(sql2_R_objects_directory, "09) Process demultiplexed PacBio reads.RData"))







# Try stuff ---------------------------------------------------------------

sg_sequences_mat <- as.matrix(library_df[, paste0("sg_sequence_", 1:4)])











# Try stuff ---------------------------------------------------------------

use_df <- ccs7_df_list[["individual_reads_df"]]
pass_filters <- use_df[["Passes_filters"]] == 1
use_zmws <- use_df[["ZMW"]][pass_filters]

zmw_matches <- match(use_zmws, ccs_df[["ZMW"]])

ccs_matched_df <- ccs_df[zmw_matches, ]
row.names(ccs_matched_df) <- NULL

TabulateBarcodes <- function(input_vec) {
  input_table <- sort(table(input_vec), decreasing = TRUE)
  num_total <- length(input_vec)
  only_one <- length(input_table) == 1
  results_list <- list(
    "first_item"      = names(input_table)[[1]],
    "first_fraction"  = input_table[[1]] / num_total,
    "second_item"     = if (only_one) NA_character_ else names(input_table)[[2]],
    "second_fraction" = if (only_one) 0 else input_table[[2]] / num_total

  )
  return(results_list)
}


plates_BC_list <- lapply(seq_len(nrow(library_df)), function(x) {
  print(library_df[["Combined_ID"]][[x]])
  are_this_ID <- ccs_matched_df[["Combined_ID"]] == library_df[["Combined_ID"]][[x]]
  barcode_1_results <- TabulateBarcodes(ccs_matched_df[["Plate_barcode_1_sequence"]][are_this_ID])
  barcode_2_results <- TabulateBarcodes(ccs_matched_df[["Plate_barcode_2_sequence"]][are_this_ID])
  names(barcode_1_results) <- paste0("BC1_", names(barcode_1_results))
  names(barcode_2_results) <- paste0("BC2_", names(barcode_2_results))
  results_list <- c(
    list("Combined_ID"  = library_df[["Combined_ID"]][[x]],
         "Plate_number" = library_df[["Plate_number"]][[x]],
         "Well_number"  = library_df[["Well"]][[x]],
         "Num_reads"    = NA
         ),
    barcode_1_results, barcode_2_results
  )
  return(results_list)
})

plates_BC_df <- do.call(rbind.data.frame,
                        c(plates_BC_list, list(stringsAsFactors = FALSE,
                                               make.row.names = FALSE
                                               )
                        ))

use_summary_df <- ccs7_df_list[["filtered_summary_df"]]
summary_matches <- match(library_df[["Combined_ID"]], use_summary_df[["Combined_ID"]])
plates_BC_df[["Num_reads"]] <- use_summary_df[["Count_total"]][summary_matches]


plates_BC_df[plates_BC_df[["filtered_summary_df"]] == "CATAGAGAGATAGTAT", ]

plates_BC1 <- tapply(plates_BC_df[["BC1_first_item"]],
                     plates_BC_df[["Plate_number"]],
                     function(x) sort(table(x, dnn = NULL), decreasing = TRUE)
                     )

fwd_barcodes <- paste0(plates_df[["Barcode_sequence"]], "T")
BC1_vec <- vapply(plates_BC1, function(x) names(x)[[1]], "")

BC1_vec == fwd_barcodes


plates_BC2 <- tapply(plates_BC_df[["BC2_first_item"]],
                     plates_BC_df[["Plate_number"]],
                     function(x) sort(table(x, dnn = NULL), decreasing = TRUE),
                     simplify = FALSE
                     )

rev_barcodes <- as.character(reverseComplement(DNAStringSet(fwd_barcodes)))
BC2_vec <- vapply(plates_BC2, function(x) names(x)[[1]], "")

BC2_vec == rev_barcodes


bc_df <- data.frame(
  "Plate_number"            = plates_df[["Plate_number"]],
  "BC1_original"            = BC1_vec,
  "BC1_complement"          = as.character(reverseComplement(DNAStringSet(BC1_vec))),
  "BC1_reverted"            = reverse(BC1_vec),
  "BC1_reverted_complement" = as.character(reverseComplement(DNAStringSet(reverse(BC1_vec)))),
  "BC1_complement_reverted" = reverse(as.character(reverseComplement(DNAStringSet(BC1_vec)))),
  "BC2_original"            = BC2_vec,
  "BC2_complement"          = as.character(reverseComplement(DNAStringSet(BC2_vec))),
  "BC2_reverted"            = reverse(BC2_vec),
  "BC2_reverted_complement" = as.character(reverseComplement(DNAStringSet(reverse(BC2_vec)))),
  stringsAsFactors = FALSE
)











