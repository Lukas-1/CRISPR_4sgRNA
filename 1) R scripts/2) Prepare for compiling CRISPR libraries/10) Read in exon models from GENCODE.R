### 7th November 2019 ###



# Import packages and source code -----------------------------------------






# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
CRISPR_input_directory   <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gff3.gz
# on 6 September 2019

GENCODE_original_df <- read.table(file.path(CRISPR_input_directory, "Human genome", "GENCODE", "gencode.v32.annotation.gff3"),
                                  sep = "\t", quote = "", stringsAsFactors = FALSE, header = FALSE, row.names = NULL,
                                  fill = TRUE, check.names = FALSE
                                  )
colnames(GENCODE_original_df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")




# Split the .gff3 attributes column ---------------------------------------

GENCODE_splits <- strsplit(GENCODE_original_df[, "attributes"], ";", fixed = TRUE)

GENCODE_split_splits <- lapply(GENCODE_splits, function(x) strsplit(x, "=", fixed = TRUE))

fields_list  <- lapply(GENCODE_split_splits, function(x) vapply(x, function(y) y[[1]], ""))
entries_list <- lapply(GENCODE_split_splits, function(x) vapply(x, function(y) y[[2]], ""))

unique_fields <- unique(unlist(fields_list))

expanded_split_splits <- lapply(seq_along(GENCODE_splits), function(x) entries_list[[x]][match(unique_fields, fields_list[[x]])])

splits_df <- do.call(rbind.data.frame, c(expanded_split_splits, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
colnames(splits_df) <- unique_fields

GENCODE_df <- data.frame(GENCODE_original_df[, !(colnames(GENCODE_original_df) %in% c("score", "attributes"))],
                         splits_df,
                         stringsAsFactors = FALSE,
                         check.names = FALSE
                         )




# Save data ---------------------------------------------------------------

save(list = "GENCODE_df",
     file = file.path(general_RData_directory, "09) Read in exon models from GENCODE.RData")
     )












