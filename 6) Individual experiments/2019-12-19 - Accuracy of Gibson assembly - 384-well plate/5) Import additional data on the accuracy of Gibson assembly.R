### 4th September 2020 ###



# Import packages and source code -----------------------------------------




# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR_4sgRNA"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
file_input_directory             <- file.path(file_directory, "1) Input")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")





# Define functions --------------------------------------------------------

ReadTable <- function(sub_path, header = TRUE, sep = "\t") {
  read.table(file.path(file_input_directory, sub_path),
             sep = sep, quote = "", stringsAsFactors = FALSE,
             header = header, row.names = NULL,
             check.names = FALSE
             )
}




# Read in data ------------------------------------------------------------

CCS3_mean_accuracies_df <- ReadTable(file.path("CCS3", "CCS3_identical.perc.txt"))
CCS5_mean_accuracies_df <- ReadTable(file.path("CCS5", "CCS5_identical.perc.txt"))

CCS3_combos_df <- ReadTable(file.path("CCS3", "CCS3_all_comb.count.txt"))
CCS5_combos_df <- ReadTable(file.path("CCS5", "CCS5_all_comb.count.txt"))

CCS3_counts_df <- ReadTable(file.path("CCS3", "CCS3_tot.count.txt"), header = FALSE, sep = "")
CCS5_counts_df <- ReadTable(file.path("CCS5", "CCS5_tot.count.txt"), header = FALSE, sep = "")

CCS5_paths <- list.files(file.path(file_input_directory,
                                   "CCS5",
                                   "Processed data FGCZ"
                                   )
                         )
CCS5_paths <- file.path("CCS5", "Processed data FGCZ", CCS5_paths)

CCS5_individual_guides_list <- lapply(CCS5_paths, ReadTable)


CC5_blastn_files <- list.files(file.path(file_input_directory,
                                         "CCS5",
                                         "blastn"
                                         )
                               )

CCS5_blastn_sg1 <- lapply(grep("db1.best.bln", CC5_blastn_files, value = TRUE, fixed = TRUE),
                          function(x) {
                            ReadTable(file.path("CCS5", "blastn", x), header = FALSE)
                          })


CC3_blastn_files <- list.files(file.path(file_input_directory,
                                         "CCS3",
                                         "blastn"
                                         )
                               )

CCS3_blastn_sg1 <- lapply(grep("db1.best.bln", CC3_blastn_files, value = TRUE, fixed = TRUE),
                          function(x) {
                            ReadTable(file.path("CCS3", "blastn", x), header = FALSE)
                          })





# Check data --------------------------------------------------------------

CCS5_all_4_correct <- vapply(CCS5_individual_guides_list,
                             function(x) sum(rowSums(as.matrix(x[, 2:5])) == 4),
                             integer(1)
                             )

files_sg1_mean_accuracy_vec <- vapply(CCS5_individual_guides_list,
                                      function(x) mean(x[[2]]),
                                      numeric(1)
                                      )
CCS5_mean_accuracies_df[["sg1_cr1_perc"]] == round(files_sg1_mean_accuracy_vec * 100, 2)


CCS5_num_reads_combos <- rowSums(CCS5_combos_df[, 2:ncol(CCS5_combos_df)])
CCS5_num_reads_individual_guides <- vapply(CCS5_individual_guides_list, nrow, integer(1))

table(CCS5_num_reads_combos > CCS5_counts_df[[2]])
table(CCS5_num_reads_individual_guides > CCS5_counts_df[[2]])


data.frame("combos_all_correct_count" = CCS5_combos_df[["sg1_cr1:sg2_cr2:sg3_cr3:sg4_cr4"]],
           "files_all_correct_count"  = CCS5_all_4_correct,
           "combos_all_correct_perc"  = CCS5_combos_df[["sg1_cr1:sg2_cr2:sg3_cr3:sg4_cr4"]] / CCS5_num_reads_combos,
           "files_all_correct_perc"   = CCS5_all_4_correct / CCS5_num_reads_individual_guides
           )




table(vapply(CCS5_blastn_sg1, nrow, integer(1)) > CCS3_counts_df[[2]])

table(vapply(CCS5_blastn_sg1, nrow, integer(1)) > CCS5_counts_df[[2]])



data.frame(CCS5_num_reads_individual_guides,
           CCS5_mean_accuracies_df[["tot"]]
           )

table(CCS5_num_reads_individual_guides < CCS5_mean_accuracies_df[["tot"]])
goo


