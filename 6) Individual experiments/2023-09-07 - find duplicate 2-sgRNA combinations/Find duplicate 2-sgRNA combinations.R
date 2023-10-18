### 7 September 2023



# Define folder paths -----------------------------------------------------

project_dir <- "~/CRISPR_4sgRNA/6) Individual experiments/2023-09-07 - find duplicate 2-sgRNA combinations"
input_dir <- file.path(project_dir, "1_input")
output_dir <- file.path(project_dir, "2_output")



# Read in data ------------------------------------------------------------

CRISPRa_df <- read.delim(file.path(input_dir, "CRISPRa_all_sublibraries_ordered_by_well.tsv"),
                         stringsAsFactors = FALSE
                         )



# Find duplicates ---------------------------------------------------------

sg_numbers <- rep(1:4, length.out = nrow(CRISPRa_df))
are_sg2 <- sg_numbers == 2
are_sg3 <- sg_numbers == 3

## Check the assumption that the 4 sgRNAs for each plasmid always follow each other
plate_splits <- strsplit(CRISPRa_df[, "Plate_string"], "_", fixed = TRUE)
plates_vec <- vapply(plate_splits, function(x) paste0(x[[1]], "-", x[[2]]), "")
plasmids_vec <- paste0(plates_vec, "_", CRISPRa_df[, "Well_number"])
stopifnot(all(table(plasmids_vec) == 4))
stopifnot(all(plasmids_vec[are_sg2] == plasmids_vec[are_sg3]))


sg2_3_vec <- paste0(toupper(CRISPRa_df[are_sg2, "sgRNA_sequence"]), "_",
                    toupper(CRISPRa_df[are_sg3, "sgRNA_sequence"])
                    )
num_occurrences_vec <- table(sg2_3_vec)[sg2_3_vec]
are_duplicated <- num_occurrences_vec > 1



# Collect results ---------------------------------------------------------

duplicates_df <- CRISPRa_df[are_sg2, ][are_duplicated, ]

duplicates_df <- data.frame(
  duplicates_df["Sublibrary_4sg"],
  "Plate" = plates_vec[are_sg2][are_duplicated],
  "Well" = duplicates_df["Well_number"],
  "Replaced_TF_gene" = NA,
  duplicates_df[, c("Is_obsolete", "Entrez_ID", "Gene_symbol", "TSS_ID", "Is_main_TSS")],
  "Number_of_target_genes" = NA,
  "sg2_sequence" = CRISPRa_df[, "sgRNA_sequence"][are_sg2][are_duplicated],
  "sg3_sequence" = CRISPRa_df[, "sgRNA_sequence"][are_sg3][are_duplicated],
  stringsAsFactors = FALSE
)
duplicates_sg2_3_vec <- paste0(toupper(duplicates_df[, "sg2_sequence"]), "_",
                               toupper(duplicates_df[, "sg3_sequence"])
                               )

new_order <- order(match(duplicates_sg2_3_vec, duplicates_sg2_3_vec))
duplicates_df <- duplicates_df[new_order, ]



# Figure out which duplicates are due to obsolete TFs ---------------------

duplicates_sg2_3_vec <- paste0(toupper(duplicates_df[, "sg2_sequence"]), "_",
                               toupper(duplicates_df[, "sg3_sequence"])
                               )
duplicates_sg2_3_fac <- factor(duplicates_sg2_3_vec, levels = unique(duplicates_sg2_3_vec))
obsolete_list <- split(duplicates_df[, "Is_obsolete"] %in% "Yes", duplicates_sg2_3_fac)
only_obsolete <- vapply(obsolete_list, function(x) (sum(x) == 1) && (length(x) == 2), logical(1))

only_obsolete_vec <- rep(only_obsolete, lengths(obsolete_list))
duplicates_df[, "Replaced_TF_gene"] <- ifelse(only_obsolete_vec, "Yes", "No")
duplicates_df[, "Number_of_target_genes"] <- rep(vapply(obsolete_list, function(x) sum(!(x)), integer(1)), lengths(obsolete_list))
duplicates_df <- duplicates_df[order(only_obsolete_vec), ]




# Export results ----------------------------------------------------------

write.csv(duplicates_df, file = file.path(output_dir, "Duplicates_sg2_sg3_Tgonfio.csv"),
          quote = FALSE, row.names = FALSE
          )



