### 19th December 2019 ###



# Import packages and source code -----------------------------------------

library("Biostrings")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
exchange_with_tools_directory    <- file.path(file_directory, "3) Exchange with other tools")





# Load data ---------------------------------------------------------------

load(file.path(intermediate_R_objects_directory, "1) Import the data on the accuracy of Gibson assembly.RData"))





# Functions for using DINAMelt --------------------------------------------

FwdAndRevCombinations <- function(four_sg_vec) {
  stopifnot(length(four_sg_vec) == 4)
  rev_sequences <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(four_sg_vec)))
  mat_list <- lapply(1:4, function(x) {
    second_sequence_vec <- c(four_sg_vec, rev_sequences[(1:4) != x])
    sub_mat <- cbind(rep(four_sg_vec[[x]], length(second_sequence_vec)), second_sequence_vec)
    return(sub_mat)
  })
  results_mat <- do.call(rbind, mat_list)
  colnames(results_mat) <- paste0("Sequence_", 1:2)
  rownames(results_mat) <- NULL
  return(results_mat)
}


OnlyFwdCombinations <- function(four_sg_vec) {
  stopifnot(length(four_sg_vec) == 4)
  results_mat <- as.matrix(expand.grid(four_sg_vec,
                                       four_sg_vec,
                                       stringsAsFactors = FALSE,
                                       KEEP.OUT.ATTRS = FALSE
                                       )
                           )[, 2:1]
  colnames(results_mat) <- paste0("Sequence_", 1:2)
  return(results_mat)
}


MatListToDf <- function(mat_list) {
  df_list <- lapply(1:384, function(x) data.frame("Gene_symbol" = assembly_df[["Gene_symbol"]][[x]],
                                                  "Well_number" = x,
                                                  mat_list[[x]],
                                                  stringsAsFactors = FALSE
                                                  )
                    )
  results_df <- do.call(rbind.data.frame, c(df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}


ExportToTwoStateMelting <- function(sequence_vec, file_prefix) {

  num_sequences_per_file <- 1000L
  num_sequences <- length(sequence_vec)
  num_files <- ceiling(num_sequences / num_sequences_per_file)

  file_sequence <- rep(seq_len(num_files), each = num_sequences_per_file)
  file_sequence <- file_sequence[seq_len(num_sequences)]

  for (i in seq_len(num_files)) {
    sequence_sub_vec <- sequence_vec[file_sequence == i]
    sequence_sub_vec <- paste0(sequence_sub_vec, c(rep(";", times = length(sequence_sub_vec) - 1), ""))
    write.table(sequence_sub_vec,
                file = file.path(exchange_with_tools_directory, "1) Input for DINAMelt", paste0(file_prefix, " - file ", i, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE
                )

  }
}






# Generate input for the DINAMelt two-state melting tool ------------------
# http://unafold.rna.albany.edu/?q=DINAMelt/Two-state-melting

sequences_mat <- as.matrix(assembly_df[, paste0("sg_", 1:4)])

only_fwd_df    <- MatListToDf(lapply(1:384, function(x) OnlyFwdCombinations(sequences_mat[x, ])))
fwd_and_rev_df <- MatListToDf(lapply(1:384, function(x) FwdAndRevCombinations(sequences_mat[x, ])))

ExportToTwoStateMelting(only_fwd_df[["Sequence_1"]], "1a) Only forward sequences")
ExportToTwoStateMelting(only_fwd_df[["Sequence_2"]], "1b) Only forward sequences")
ExportToTwoStateMelting(fwd_and_rev_df[["Sequence_1"]], "2a) Forward and reverse sequences")
ExportToTwoStateMelting(fwd_and_rev_df[["Sequence_2"]], "2b) Forward and reverse sequences")






# Export sequences for viennaRNA self-annealing prediction ----------------

single_sequences <- unlist(lapply(paste0("sg_", 1:4), function(x) assembly_df[[x]]))

write.table(single_sequences,
            file = file.path(exchange_with_tools_directory, "Single sequences for self-annealing prediction.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE
            )






# Save data ---------------------------------------------------------------

save(list = c("only_fwd_df", "fwd_and_rev_df"),
     file = file.path(intermediate_R_objects_directory,
                      "2) Export data for external tools (e.g. melting temperature).RData"
                      )
     )




