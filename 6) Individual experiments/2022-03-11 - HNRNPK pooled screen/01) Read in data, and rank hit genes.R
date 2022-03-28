### 11th March 2022 ###



# Import packages and source code -----------------------------------------

library("readxl")
library("edgeR")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

file_directory        <- file.path(experiments_directory, "2022-03-11 - HNRNPK pooled screen")
file_input_directory  <- file.path(file_directory, "1) Input")
R_objects_directory   <- file.path(file_directory, "2) R objects")

screen_data_directory <- file.path(file_input_directory, "Screen data")
comparisons_folder    <- file.path(screen_data_directory, "Sample comparisons")




# Define functions --------------------------------------------------------

TidyComparisonDf <- function(input_df, library_df) {
  are_unique <- !(duplicated(as.list(input_df)))
  input_df <- input_df[, are_unique]
  are_all_same <- vapply(seq_along(input_df), function(x) {
    length(unique(input_df[[x]])) == 1
  }, logical(1))
  input_df <- input_df[, !(are_all_same)]
  library_vec <- paste0(library_df[, "Target Transcript"], "-", library_df[, "sgRNA Target Sequence"])
  matches_vec <- match(input_df[, "gene_id"], library_vec)
  stopifnot(!(anyNA(is.na(matches_vec))))
  input_df[, "Gene_symbol"] <- library_df[, "Target Gene Symbol"][matches_vec]
  input_df[, "Entrez_gene_ID"] <- library_df[, "Target Gene ID"][matches_vec]
  return(input_df)
}

ReadInComparison <- function(file_name) {
  data_df <- data.frame(read_excel(file.path(comparisons_folder, paste0(file_name, ".xlsx"))),
                        stringsAsFactors = FALSE, check.names = FALSE
                        )
  data_df <- TidyComparisonDf(data_df, Brunello_df)
  return(data_df)
}


# Read in the library -----------------------------------------------------

Brunello_2016_df <- data.frame(read_excel(file.path(file_input_directory, "Library",
                                                    "2016 - Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9 - Table S21 - Brunello.xlsx"
                                                    )),
                               stringsAsFactors = FALSE, check.names = FALSE
                               )

Brunello_2018_df <- data.frame(read_excel(file.path(file_input_directory, "Library",
                                                    "2018 - Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities - Data S1.xlsx"
                                                    ),
                                          sheet = "sgRNA annotations"
                                          ),
                               stringsAsFactors = FALSE, check.names = FALSE
                               )



# Read in gene selections -------------------------------------------------

Stefano_df <- read_excel(file.path(file_input_directory, "Screen data",
                                   "Stefano's Excel sheets",
                                   "Validation_Screen_first gene set.xlsx"
                                   ))




# Read in raw data --------------------------------------------------------

raw_dir <- path.expand(file.path(screen_data_directory, "Raw counts"))
raw_paths <- list.files(raw_dir, full.names = TRUE)

counts_mat <- do.call(cbind, lapply(raw_paths, function(x) read.delim(x)[, "Count"]))
file_names <- list.files(raw_dir)
file_names <- sub("o25644_1_", "", file_names, fixed = TRUE)
file_names <- sub("-result.txt", "", file_names, fixed = TRUE)
colnames(counts_mat) <- file_names

raw_df_list <- unique(lapply(raw_paths, function(x) {
  results_df <- read.delim(x, stringsAsFactors = FALSE)
  results_df[, names(results_df) != "Count"]
}))
stopifnot(length(raw_df_list) == 1)
raw_df <- data.frame(raw_df_list[[1]], counts_mat, check.names = FALSE)




# Include control sgRNAs in the Brunello library  -------------------------

are_controls <- !(Brunello_2018_df[, "sgRNA Sequence"] %in% Brunello_2016_df[, "sgRNA Target Sequence"])
Brunello_2018_df[are_controls, ]

controls_df <- data.frame(matrix(nrow = sum(are_controls),
                                 ncol = ncol(Brunello_2016_df),
                                 dimnames = list(NULL, names(Brunello_2016_df))
                                 ), check.names = FALSE)
controls_df[, "sgRNA Target Sequence"] <- Brunello_2018_df[, "sgRNA Sequence"][are_controls]
controls_df[, "Target Transcript"] <- "Control"
Brunello_df <- rbind.data.frame(Brunello_2016_df, controls_df,
                                stringsAsFactors = FALSE, make.row.names = FALSE
                                )



# Read in comparison data -------------------------------------------------

NT_v_HNRNPK_df   <- ReadInComparison("result--D14_HNRNPK--over--D14_NT")
base_v_NT_df     <- ReadInComparison("result--D14_NT--over--D1")
base_v_HNRNPK_df <- ReadInComparison("result--D14_HNRNPK--over--D1")




# Perform checks ----------------------------------------------------------

table(raw_df[, "GeneSymbol"] %in% NT_v_HNRNPK_df[, "Gene_symbol"])
table(Brunello_df[, "Target Gene Symbol"] %in% raw_df[, "GeneSymbol"])




# Compare read counts across samples --------------------------------------

are_NT <- grepl("Non-Targeting", raw_df[, "GeneSymbol"], fixed = TRUE)

cbind("NT"    = apply(counts_mat, 2, function(x) sum(x[are_NT])),
      "Total" = colSums(counts_mat)
      )

dgList <- edgeR::calcNormFactors(counts_mat, method = "TMM")
dgList
edgeR::cpm(dgList)





# Define intersection genes -----------------------------------------------

pass_criteria <- (NT_v_HNRNPK_df[, "log2 Ratio"] > 1) &
                 (NT_v_HNRNPK_df[, "fdr"] < 0.01)

NT_v_HNRNPK_df[["Is_hit"]] <- pass_criteria




# Summarize (intersected) upregulated gRNAs on a per-gene level -----------

genes_fac <- factor(NT_v_HNRNPK_df[pass_criteria, "Gene_symbol"])

up_genes_df <- data.frame(
  "Gene_symbol" = levels(genes_fac),
  "Number_of_guides" = tabulate(genes_fac),
  stringsAsFactors = FALSE
)

new_order <- order(up_genes_df[, "Number_of_guides"], decreasing = TRUE)

up_genes_df <- up_genes_df[new_order, ]
row.names(up_genes_df) <- NULL




# Define Stefano's gene selections ----------------------------------------

Stefano_top_genes           <- Stefano_df[[2]][71:89]
Stefano_non_essential_genes <- Stefano_df[[5]][71:92]
Stefano_top_genes           <- sapply(strsplit(Stefano_top_genes, " ", fixed = TRUE), "[[", 1)
Stefano_non_essential_genes <- sapply(strsplit(Stefano_non_essential_genes, " ", fixed = TRUE), "[[", 1)





# Save data ---------------------------------------------------------------

save(list = c("NT_v_HNRNPK_df", "base_v_NT_df", "base_v_HNRNPK_df",
              "raw_df", "up_genes_df",
              "Stefano_top_genes", "Stefano_non_essential_genes"
              ),
     file = file.path(R_objects_directory, "01) Read in data, and rank hit genes.RData")
     )





