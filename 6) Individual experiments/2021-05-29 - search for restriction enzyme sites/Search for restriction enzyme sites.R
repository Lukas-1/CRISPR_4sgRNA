### 22nd December 2020 ###




# Import packages and source code -----------------------------------------

library("circRNAprofiler") # For getRegexPattern

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "14) Checking for identical subsequences.R")) # For CheckThatFactorIsInOrder
source(file.path(general_functions_directory, "17) Exporting CRISPR libraries as text files.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")

file_directory           <- file.path(CRISPR_root_directory, "6) Individual experiments/2021-05-29 - search for restriction enzyme sites")
file_input_directory     <- file.path(file_directory, "1) Input")
RData_directory          <- file.path(file_directory, "2) RData files")
file_output_directory    <- file.path(file_directory, "3) Output")




# Load data ---------------------------------------------------------------

load(file.path(CRISPRa_RData_directory, "28) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRa_by_gene_df <- full_4sg_by_gene_df
load(file.path(CRISPRko_RData_directory, "20) Distribute sgRNAs for the whole genome onto plates.RData"))
CRISPRko_by_gene_df <- full_4sg_by_gene_df
rm(list = c("full_4sg_by_gene_df", "full_4sg_by_well_df",
            "sg4_by_gene_df", "sg4_by_well_df", "shared_sgRNAs_df"
            )
   )




# Read in data ------------------------------------------------------------

unique_cutters <- read.table(file.path(file_input_directory, "unique_cutters.txt"),
                             header = TRUE,
                             check.names = FALSE, sep = "\t",
                             stringsAsFactors = FALSE
                             )
unique_cutters <- unlist(strsplit(unique_cutters[, 1], ", ?"))


rebase_file_df <- read.table(file.path(file_input_directory, "link_allenz"),
                             sep = "\n", header = FALSE,
                             stringsAsFactors = FALSE,
                             quote = "\"", comment.char = ""
                             )




# Define functions --------------------------------------------------------

WriteTable <- function(export_df, file_path) {
  are_character <- vapply(export_df, is.character, logical(1))
  for (i in which(are_character)) {
    export_df[[i]] <- ifelse(is.na(export_df[[i]]), "", export_df[[i]])
  }
  write.table(export_df,
              file      = file_path,
              sep       = "\t",
              row.names = FALSE,
              col.names = TRUE,
              quote     = FALSE
              )
}



FindSequenceHits <- function(CRISPR_df, regex_vec) {
  ID_vec <- CRISPR_df[["Combined_ID"]]
  if ("AltTSS_ID" %in% names(CRISPR_df)) {
    ID_vec <- ifelse(CRISPR_df[["Is_control"]] == "Yes", ID_vec, CRISPR_df[["AltTSS_ID"]])
  }
  IDs_fac <- factor(ID_vec, levels = unique(ID_vec))
  CheckThatFactorIsInOrder(IDs_fac)
  stopifnot(length(unique(split(CRISPR_df[["Rank"]], IDs_fac))) == 1)

  are_this_rank_list <- lapply(1:4, function(x) CRISPR_df[["Rank"]] %in% x)
  num_plasmids <- unique(vapply(are_this_rank_list, sum, integer(1)))
  stopifnot(identical(num_plasmids, nlevels(IDs_fac)))

  are_hits_list <- lapply(regex_vec, function(x) grepl(x, CRISPR_df[["Sequence_with_primers"]]))

  are_hits_mat_list <- lapply(are_hits_list, function(x) do.call(cbind, lapply(are_this_rank_list, function(y) x[y])))
  any_hits_list <- lapply(are_hits_mat_list, function(x) rowSums(x))
  num_hits_mat <- do.call(rbind, any_hits_list)
  mode(num_hits_mat) <- "integer"

  num_hits_vec <- as.integer(rowSums(num_hits_mat > 0))
  colnames(num_hits_mat) <- levels(IDs_fac)

  results_mat <- cbind("Num_cut_plasmids" = num_hits_vec,
                       "Total_num_plasmids" = num_plasmids,
                       num_hits_mat
                       )
  return(results_mat)
}





# Add flanking sequences (including the primer) ---------------------------

CRISPRa_by_gene_df[["Sequence_with_primers"]]  <- toupper(AddPrimers(CRISPRa_by_gene_df))
CRISPRko_by_gene_df[["Sequence_with_primers"]] <- toupper(AddPrimers(CRISPRko_by_gene_df))





# Parse the REBASE text file ----------------------------------------------

are_field_1 <- grepl("^<1>", rebase_file_df[, 1])
are_references_start <- rebase_file_df[, 1] == "References:"

rebase_indices <- seq(from = which(are_field_1)[[1]],
                      to = which(are_references_start)
                      )
enzyme_number <- 0L
enzyme_number_vec <- rep(NA, nrow(rebase_file_df))
for (i in rebase_indices) {
  if (are_field_1[[i]]) {
    enzyme_number <- enzyme_number + 1L
  }
  enzyme_number_vec[[i]] <- enzyme_number
}

enzyme_splits <- split(rebase_file_df[rebase_indices, 1], enzyme_number_vec[rebase_indices])

enzyme_names <- vapply(enzyme_splits, function(x) grep("^<1>", x, value = TRUE), "")
enzyme_names <- sub("<1>", "", enzyme_names, fixed = TRUE)

recognition_sites <- vapply(enzyme_splits, function(x) grep("^<5>", x, value = TRUE), "")
recognition_sites <- toupper(sub("<5>", "", recognition_sites, fixed = TRUE))

recognition_sites_only <- ifelse(recognition_sites == "?", NA, recognition_sites)
recognition_sites_only <- gsub("\\([^)]+\\)", "", recognition_sites_only)
recognition_sites_only <- sub("^", "", recognition_sites_only, fixed = TRUE)

methylation_sites <- vapply(enzyme_splits, function(x) grep("^<6>", x, value = TRUE), "")
methylation_sites <- sub("<6>", "", methylation_sites, fixed = TRUE)
methylation_sites <- ifelse(methylation_sites == "", NA, methylation_sites)

rebase_df <- data.frame(
  "Enzyme_name" = enzyme_names,
  "Recognition_site" = recognition_sites_only,
  "Recognition_and_cut_site" = recognition_sites,
  "Methylation_site" = methylation_sites,
  stringsAsFactors = FALSE
)

rebase_df <- rebase_df[!(is.na(rebase_df[["Recognition_site"]])), ]
row.names(rebase_df) <- NULL






# Search for matches in the CRISPR libraries ------------------------------

rebase_df[["Is_unique_cutter"]] <- rebase_df[["Enzyme_name"]] %in% unique_cutters

site_splits <- strsplit(rebase_df[["Recognition_site"]], ",", fixed = TRUE)
rebase_df[["Recognition_regex"]] <- vapply(site_splits,
                                           function(x) {
                                             if ((length(x) == 1) && (is.na(x))) {
                                               NA_character_
                                             } else {
                                               regex_vec <- vapply(x, function(y) circRNAprofiler::getRegexPattern(y, isDNA = TRUE), "", USE.NAMES = FALSE)
                                               paste0(regex_vec, collapse = "|")
                                             }},
                                           "", USE.NAMES = FALSE
                                           )

CRISPRa_hits_mat  <- FindSequenceHits(CRISPRa_by_gene_df,  rebase_df[["Recognition_regex"]])
CRISPRko_hits_mat <- FindSequenceHits(CRISPRko_by_gene_df, rebase_df[["Recognition_regex"]])

CRISPRa_rebase_df  <- data.frame(rebase_df["Enzyme_name"], CRISPRa_hits_mat,
                                 check.names = FALSE, stringsAsFactors = FALSE
                                 )
CRISPRko_rebase_df <- data.frame(rebase_df["Enzyme_name"], CRISPRko_hits_mat,
                                 check.names = FALSE, stringsAsFactors = FALSE
                                 )

colnames(CRISPRa_hits_mat)[1:2] <- paste0("CRISPRa_", tolower(colnames(CRISPRa_hits_mat)[1:2]))
colnames(CRISPRko_hits_mat)[1:2] <- paste0("CRISPRko_", tolower(colnames(CRISPRko_hits_mat)[1:2]))

rebase_df <- data.frame(rebase_df,
                        CRISPRa_hits_mat[, 1:2],
                        CRISPRko_hits_mat[, 1:2],
                        stringsAsFactors = FALSE
                        )





# Check for co-occurrence of cut sites ------------------------------------

CRISPRa_cut_by_PmeI <- CRISPRa_hits_mat[rebase_df[["Enzyme_name"]] == "PmeI", 3:ncol(CRISPRa_hits_mat)] > 0
CRISPRa_cut_by_FseI <- CRISPRa_hits_mat[rebase_df[["Enzyme_name"]] == "FseI", 3:ncol(CRISPRa_hits_mat)] > 0

CRISPRko_cut_by_PmeI <- CRISPRko_hits_mat[rebase_df[["Enzyme_name"]] == "PmeI", 3:ncol(CRISPRko_hits_mat)] > 0
CRISPRko_cut_by_FseI <- CRISPRko_hits_mat[rebase_df[["Enzyme_name"]] == "FseI", 3:ncol(CRISPRko_hits_mat)] > 0

sum(CRISPRa_cut_by_PmeI)
sum(CRISPRa_cut_by_FseI)

sum(CRISPRko_cut_by_PmeI)
sum(CRISPRko_cut_by_FseI)

table(CRISPRa_cut_by_PmeI, CRISPRa_cut_by_FseI)
table(CRISPRko_cut_by_PmeI, CRISPRko_cut_by_FseI)




# Re-order the REBASE data frame ------------------------------------------

new_order <- order(!(rebase_df[["Is_unique_cutter"]]),
                   rowSums(as.matrix(rebase_df[, c("CRISPRa_num_cut_plasmids", "CRISPRko_num_cut_plasmids")]))
                   )
rebase_reordered_df <- rebase_df[new_order, ]
row.names(rebase_reordered_df) <- NULL





# Export data -------------------------------------------------------------

rebase_export_df <- rebase_reordered_df[, !(names(rebase_reordered_df) %in% "Recognition_regex")]
rebase_export_df[["Is_unique_cutter"]] <- ifelse(rebase_export_df[["Is_unique_cutter"]], "Yes", "No")
WriteTable(rebase_export_df, file.path(file_output_directory, "REBASE_4sg_plasmids.tsv"))





# Save data ---------------------------------------------------------------

save(list = c("rebase_df", "rebase_reordered_df",
              "CRISPRa_rebase_df", "CRISPRko_rebase_df"
              ),
     file = file.path(RData_directory, "Search for restriction enzyme sites.RData")
     )





