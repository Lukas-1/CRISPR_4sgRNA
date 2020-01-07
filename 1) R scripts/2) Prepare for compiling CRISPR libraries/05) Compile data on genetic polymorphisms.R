### 7th August 2019 ###





# Import packages and source code -----------------------------------------





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "03) Read in data on genetic polymorphisms from dbSNP.RData"))
load(file.path(general_RData_directory, "04) Read in data on genetic polymorphisms from Kaviar.RData"))





# Define functions --------------------------------------------------------

CharacterVecToAFList <- function(AF_char_vec, prefix) {

  AF_vec <- gsub("(?<=,)\\.(?=(,|$))", ",", AF_char_vec, perl = TRUE)
  AF_vec <- gsub(",+$", "", AF_vec)
  AF_vec <- gsub(",,+", ",", AF_vec)
  AF_vec <- sub(paste0(prefix, ".,"), prefix, AF_vec, fixed = TRUE)

  AF_splits <- strsplit(AF_vec, ",", fixed = TRUE)
  AF_1_vec <- as.numeric(sub(prefix, "", sapply(AF_splits, "[[", 1), fixed = TRUE))
  AF_rest <- lapply(AF_splits, function(x) if (length(x) > 1) as.numeric(x[2:length(x)]) else numeric(0))

  AF_list <- mapply(c, as.list(AF_1_vec), AF_rest)
  return(AF_list)

}



# Process data ------------------------------------------------------------

AF1kGenomes_list <- CharacterVecToAFList(dbSNP_common_raw_df[, "AF_1kGenomes"], "CAF=")
TOPMED_list      <- CharacterVecToAFList(dbSNP_common_raw_df[, "AF_TOPMED"], "TOPMED=")
Kaviar_list      <- CharacterVecToAFList(Kaviar_common_raw_df[, "AF_Kaviar"], "AF=")

AF_1k_vec  <- sapply(AF1kGenomes_list, "[[", 1)
AF_TOPMED_vec <- sapply(TOPMED_list, "[[", 1)
AF_Kaviar_vec <- vapply(Kaviar_list, function(x) 1 - sum(x), numeric(1))




# Merge and build the data frame ------------------------------------------

shared_rsIDs <- intersect(dbSNP_common_raw_df[, "rsID"], Kaviar_common_raw_df[, "rsID"])

shared_rsIDs_indices_dbSNP  <- match(shared_rsIDs, dbSNP_common_raw_df[, "rsID"])
shared_rsIDs_indices_Kaviar <- match(shared_rsIDs, Kaviar_common_raw_df[, "rsID"])

dbSNP_indices <- seq_len(nrow(dbSNP_common_raw_df))
unique_indices_dbSNP <- setdiff(dbSNP_indices, shared_rsIDs_indices_dbSNP)

Kaviar_indices <- seq_len(nrow(Kaviar_common_raw_df))
unique_indices_Kaviar <- setdiff(Kaviar_indices, shared_rsIDs_indices_Kaviar)

common_polymorphisms_df <- rbind.data.frame(
  data.frame(
    dbSNP_common_raw_df[shared_rsIDs_indices_dbSNP, c("Chromosome", "Position", "rsID", "Reference")],
    "AF_1kGenomes" = AF_1k_vec[shared_rsIDs_indices_dbSNP],
    "AF_TOPMED"    = AF_TOPMED_vec[shared_rsIDs_indices_dbSNP],
    "AF_Kaviar"    = AF_Kaviar_vec[shared_rsIDs_indices_Kaviar],
    stringsAsFactors = FALSE,
    row.names = NULL
  ),
  data.frame(
    dbSNP_common_raw_df[unique_indices_dbSNP, c("Chromosome", "Position", "rsID", "Reference")],
    "AF_1kGenomes" = AF_1k_vec[unique_indices_dbSNP],
    "AF_TOPMED"    = AF_TOPMED_vec[unique_indices_dbSNP],
    "AF_Kaviar"    = NA_real_,
    stringsAsFactors = FALSE,
    row.names = NULL
  ),
  data.frame(
    Kaviar_common_raw_df[unique_indices_Kaviar, c("Chromosome", "Position", "rsID", "Reference")],
    "AF_1kGenomes" = NA_real_,
    "AF_TOPMED"    = NA_real_,
    "AF_Kaviar"    = AF_Kaviar_vec[unique_indices_Kaviar],
    stringsAsFactors = FALSE,
    row.names = NULL
  ),
  make.row.names = FALSE,
  stringsAsFactors = FALSE
)

common_polymorphisms_df <- common_polymorphisms_df[order(match(common_polymorphisms_df[, "Chromosome"], unique(dbSNP_common_raw_df[, "Chromosome"])),
                                                         common_polymorphisms_df[, "Position"],
                                                         common_polymorphisms_df[, "rsID"]
                                                         ), ]
row.names(common_polymorphisms_df) <- NULL





# Perform checks and tests ------------------------------------------------

gc()

dbSNP_chrom_pos_ref_vec  <- paste0(dbSNP_common_raw_df[shared_rsIDs_indices_dbSNP, "Chromosome"], "_",
                                   dbSNP_common_raw_df[shared_rsIDs_indices_dbSNP, "Position"],   "_",
                                   dbSNP_common_raw_df[shared_rsIDs_indices_dbSNP, "Reference"]
                                   )
Kaviar_chrom_pos_ref_vec <- paste0(Kaviar_common_raw_df[shared_rsIDs_indices_Kaviar, "Chromosome"], "_",
                                   Kaviar_common_raw_df[shared_rsIDs_indices_Kaviar, "Position"],   "_",
                                   Kaviar_common_raw_df[shared_rsIDs_indices_Kaviar, "Reference"]
                                   )
chrom_pos_ref_df <- data.frame("dbSNP"  = dbSNP_chrom_pos_ref_vec,
                               "Kaviar" = Kaviar_chrom_pos_ref_vec,
                               "rsID"   = dbSNP_common_raw_df[shared_rsIDs_indices_dbSNP, "rsID"],
                               stringsAsFactors = FALSE,
                               row.names = NULL
                               )
table(chrom_pos_ref_df[, 1] != chrom_pos_ref_df[, 2])
table(dbSNP_common_raw_df[, "rsID"] %in% Kaviar_common_raw_df[, "rsID"])





# Save data ---------------------------------------------------------------

save(list = "common_polymorphisms_df",
     file = file.path(general_RData_directory, "05) Compile data on genetic polymorphisms.RData")
     )









