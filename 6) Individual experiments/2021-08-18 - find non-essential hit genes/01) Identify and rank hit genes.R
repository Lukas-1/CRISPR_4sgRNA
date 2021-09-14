### 10 September 2021 ###



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

file_directory        <- file.path(experiments_directory, "2021-08-18 - find non-essential hit genes")
file_input_directory  <- file.path(file_directory, "1) Input")
screen_data_directory <- file.path(file_input_directory, "Pooled screen", "Data")
R_objects_directory   <- file.path(file_directory, "2) R objects")



# Read in data ------------------------------------------------------------

YM5_v_DMSO5_df <- read.csv(file.path(screen_data_directory, "20210817_SVS-YM5_DMSO5_Final.csv"),
                           stringsAsFactors = FALSE, check.names = FALSE
                           )

NT_v_KO_df <- read.delim(file.path(screen_data_directory, "SVS_KO_bp14-b4_1_0.01_Sidak_sgRNAList.txt"),
                         stringsAsFactors = FALSE, check.names = FALSE
                         )



# Define intersection genes -----------------------------------------------

are_up <- NT_v_KO_df[, "fold change"] > 2
up_in_KO_genes <- unique(NT_v_KO_df[["gene"]][are_up])

target_up_genes <- YM5_v_DMSO5_df[, "Gene_symbol"] %in% up_in_KO_genes

pass_criteria <- (YM5_v_DMSO5_df[, "log2 Ratio"] > 1) &
                 (YM5_v_DMSO5_df[, "fdr"] < 0.01)

are_selected <- pass_criteria & target_up_genes

YM5_v_DMSO5_df[["Is_hit"]] <- are_selected




# Summarize (intersected) upregulated gRNAs on a per-gene level -----------

genes_fac <- factor(YM5_v_DMSO5_df[are_selected, "Gene_symbol"])

up_genes_df <- data.frame(
  "Gene_symbol" = levels(genes_fac),
  "Number_of_guides" = tabulate(genes_fac),
  stringsAsFactors = FALSE
)

new_order <- order(up_genes_df[, "Number_of_guides"], decreasing = TRUE)

up_genes_df <- up_genes_df[new_order, ]
row.names(up_genes_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = c("YM5_v_DMSO5_df", "up_genes_df"),
     file = file.path(R_objects_directory, "01) Identify and rank hit genes.RData")
     )





