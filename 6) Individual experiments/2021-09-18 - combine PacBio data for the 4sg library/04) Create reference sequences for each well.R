### 21st September 2021 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "03) Creating reference plasmid sequences.R"))



# Define folder paths -----------------------------------------------------

s2r1_directory           <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
s2rC_directory           <- file.path(experiments_directory, "2021-09-18 - combine PacBio data for the 4sg library")

p1_R_objects_directory   <- file.path(plate1_directory, "3) R objects")
s2r1_R_objects_directory <- file.path(s2r1_directory, "3) R objects")
s2rC_R_objects_directory <- file.path(s2rC_directory, "3) R objects")

file_output_directory    <- file.path(s2rC_directory, "5) Output")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p1_R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData"))
load(file.path(s2r1_R_objects_directory, "04) Create reference sequences for each well - annotated plasmid.RData"))
load(file.path(s2rC_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rC_R_objects_directory, "03) Import and process sgRNA sequences.RData"))




# Define the reference sequences ------------------------------------------

sg_sequences_df <- AddReferenceSequences(library_df,
                                         tracRNAs_vec,
                                         promoters_vec,
                                         plasmid_string
                                         )



# Insert barcodes and guides into the annotated plasmids ------------------

plate_matches <- match(sg_sequences_df[["Plate_name"]], plates_df[["Plate_name"]])
plate_barcodes <- plates_df[["Barcode_sequence"]][plate_matches]
rev_plate_barcodes <- as.character(reverseComplement(DNAStringSet(plate_barcodes)))

column_barcodes <- column_bc_vec[sg_sequences_df[["Well_number"]]]
row_barcodes <- row_bc_vec[sg_sequences_df[["Well_number"]]]

library_seq <- seq_len(nrow(library_df))
plasmid_lines_list <- lapply(library_seq, function(x) {
  use_lines <- plasmid_gbk[["sequence"]]
  sg_seqs <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  replace_seqs <- c(plate_barcodes[[x]],
                    column_barcodes[[x]],
                    sg_seqs,
                    row_barcodes[[x]],
                    rev_plate_barcodes[[x]]
                    )
  ReplaceNNNLines(use_lines, replace_seqs)
})



# Create the full plasmid sequences ---------------------------------------

barcoded_plasmids <- paste0(plate_barcodes,
                            "T",
                            column_barcodes,
                            sg_sequences_df[, "Whole_plasmid"],
                            row_barcodes,
                            "A",
                            rev_plate_barcodes
                            )
sg_sequences_df[["Barcoded_plasmid"]] <- barcoded_plasmids

plasmids_NNN <- paste0(plate_barcodes, "T", column_barcodes, plasmid_string,
                       row_barcodes, "A", rev_plate_barcodes
                       )
sg_sequences_df[["Barcoded_plasmid_NNN"]] <- plasmids_NNN




# Export the plain plasmid sequences --------------------------------------

plate_number_strings <- paste0("Plate", formatC(plates_df[["Plate_number"]], width = 3, flag = "0"))

fasta_titles <- paste0(">", sg_sequences_df[["Combined_ID"]])
fasta_list <- lapply(library_seq,
                     function(x) c(fasta_titles[[x]], barcoded_plasmids[[x]], "")
                     )

plain_path <- file.path(file_output_directory, "Plasmid sequences (plain)")

for (plate_number in unique(sg_sequences_df[["Plate_number"]])) {
  are_this_plate <- library_df[["Plate_number"]] == plate_number
  file_name <- paste0(plate_number_strings[[which(plates_df[["Plate_number"]] == plate_number)]], "_reference_sequences.fa")
  write.table(unlist(fasta_list[are_this_plate]),
              file = file.path(plain_path, file_name),
              quote = FALSE, row.names = FALSE, col.names = FALSE
              )
}



# Export the annotated plasmids -------------------------------------------

for (plate_number in unique(sg_sequences_df[["Plate_number"]])) {
  are_this_plate <- library_df[["Plate_number"]] == plate_number
  folder_path <- file.path(file_output_directory,
                           "Reference plasmid sequences (annotated)",
                           plate_number_strings[[which(plates_df[["Plate_number"]] == plate_number)]]
                           )
  dir.create(folder_path, showWarnings = FALSE)
  for (i in which(are_this_plate)) {
    file_name <- paste0(library_df[["Combined_ID"]][[i]], "_barcoded_cassette.gbk")
    export_lines <- c(plasmid_gbk[["preceding"]],
                      plasmid_lines_list[[i]],
                      plasmid_gbk[["following"]]
                      )
    write.table(x         = export_lines,
                file      = file.path(folder_path, file_name),
                col.names = FALSE,
                row.names = FALSE,
                quote     = FALSE
                )
  }
}



# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(s2rC_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData")
     )



