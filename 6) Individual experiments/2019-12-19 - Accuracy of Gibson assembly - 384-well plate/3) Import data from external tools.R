### 19th December 2019 ###





# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory            <- "~/CRISPR"
file_directory                   <- file.path(CRISPR_root_directory, "6) Individual experiments/2019-12-19 - Accuracy of Gibson assembly - 384-well plate")
intermediate_R_objects_directory <- file.path(file_directory, "2) Intermediate R objects")
exchange_with_tools_directory    <- file.path(file_directory, "3) Exchange with other tools")




# Load data ---------------------------------------------------------------

load(file.path(intermediate_R_objects_directory, "1) Import the data on the accuracy of Gibson assembly.RData"))
load(file.path(intermediate_R_objects_directory, "2) Export data for external tools (e.g. melting temperature).RData"))


# This RData object comes the R script
# "Compute free energies for single sgRNA sequences",
# run on a Ubuntu virtual machine
load(file.path(exchange_with_tools_directory, "ViennaRNA_minimum_free_energies.RData"))




# Define functions --------------------------------------------------------

ReadTwoStateMeltingOutput <- function(file_prefix) {
  # Read in the output from the DINAMelt two-state melting tool

  all_files <- list.files(file.path(exchange_with_tools_directory, "2) Output from DINAMelt"))
  selected_files <- grep(file_prefix, all_files, fixed = TRUE, value = TRUE)

  file_numbers <- as.integer(sub("file ", "", sapply(strsplit(selected_files, " - ", fixed = TRUE), "[[", 2), fixed = TRUE))
  selected_files <- selected_files[order(file_numbers)]

  output_mat_list <- lapply(selected_files, function(x) {
    output_df <- read.table(file.path(exchange_with_tools_directory, "2) Output from DINAMelt", x),
                            sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL
                            )[, 2:5]

    for (i in 1:4) {
      output_df[[i]] <- as.numeric(gsub("[^-+.0-9]", "", output_df[[i]]))
    }
    output_mat <- as.matrix(output_df)
    return(output_mat)
  })
  results_mat <- do.call(rbind, output_mat_list)
  colnames(results_mat) <- c("Gibbs_free_energy", "enthalpy", "entropy", "Tm_Celcius")
  return(results_mat)
}



# Read in data ------------------------------------------------------------

only_fwd_melting_mat <- ReadTwoStateMeltingOutput("1) Only forward sequences")
fwd_and_rev_melting_mat <- ReadTwoStateMeltingOutput("2) Forward and reverse sequences")






# Process output from the DINAMelt two-state melting tool -----------------

only_fwd_melting_df <- data.frame(only_fwd_df,
                                  only_fwd_melting_mat,
                                  stringsAsFactors = FALSE,
                                  row.names = NULL
                                  )

fwd_and_rev_melting_df <- data.frame(fwd_and_rev_df,
                                     fwd_and_rev_melting_mat,
                                     stringsAsFactors = FALSE,
                                     row.names = NULL
                                     )


for (column_name in colnames(only_fwd_melting_mat)[1:3]) {
  assembly_df[[paste0("Only_fwd_", column_name)]] <- tapply(only_fwd_melting_df[[column_name]],
                                                            only_fwd_melting_df[["Well_number"]],
                                                            min
                                                            )
}
assembly_df[["Only_fwd_Tm_Celcius"]] <- tapply(only_fwd_melting_df[["Tm_Celcius"]],
                                               only_fwd_melting_df[["Well_number"]],
                                               max
                                               )


for (column_name in colnames(fwd_and_rev_melting_mat)[1:3]) {
  assembly_df[[paste0("Fwd_and_rev_", column_name)]] <- tapply(fwd_and_rev_melting_df[[column_name]],
                                                               fwd_and_rev_melting_df[["Well_number"]],
                                                               min
                                                               )
}
assembly_df[["Fwd_and_rev_Tm_Celcius"]] <- tapply(fwd_and_rev_melting_df[["Tm_Celcius"]],
                                                  fwd_and_rev_melting_df[["Well_number"]],
                                                  max
                                                  )





# Process minimum free energy (MFE) scores from ViennaRNA -----------------

ViennaRNA_mat <- matrix(minimum_free_energies_vec,
                        nrow = nrow(assembly_df),
                        dimnames = list(NULL, paste0("MFE_sg", 1:4))
                        )

ViennaRNA_mat <- cbind(ViennaRNA_mat,
                       "MinMFE_4sg" = apply(ViennaRNA_mat, 1, min),
                       "MinMFE_3sg" = apply(ViennaRNA_mat[, 2:4], 1, min)
                       )

assembly_df <- data.frame(assembly_df, ViennaRNA_mat, stringsAsFactors = FALSE, row.names = NULL)






# Build a data frame for individual guides --------------------------------

individual_guides_df <- data.frame(
  "Gene_symbol"      = rep(assembly_df[["Gene_symbol"]], 4),
  "Sequence"         = as.vector(as.matrix(assembly_df[, paste0("sg_", 1:4)])),
  "Num_reads"        = rep(assembly_df[["Num_reads"]], 4),
  "Fraction_correct" = as.vector(as.matrix(assembly_df[, paste0("Fraction_correct_sg", 1:4)])),
  "MFE"              = minimum_free_energies_vec,
  stringsAsFactors   = FALSE,
  row.names          = NULL
)





# Save data ---------------------------------------------------------------

save(list = c("assembly_df", "individual_guides_df",
              "only_fwd_melting_df", "fwd_and_rev_melting_df"
              ),
     file = file.path(intermediate_R_objects_directory,
                      "3) Import data from external tools.RData"
                      )
     )







