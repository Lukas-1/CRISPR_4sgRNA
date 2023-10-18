### 20th October 2020 ###




# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
plate1_directory      <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "03) Creating reference plasmid sequences.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory       <- "~/CRISPR_4sgRNA"
plate2_directory            <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-09-18 - PacBio - second 384-well plate")

p1_R_objects_directory      <- file.path(plate1_directory, "3) R objects")
p2_R_objects_directory      <- file.path(plate2_directory, "2) R objects")

annotated_plasmid_directory <- file.path(plate1_directory, "2) Input", "Annotated plasmids")

file_output_directory       <- file.path(plate2_directory, "3) Output")
reference_output_directory  <- file.path(file_output_directory, "Reference sequences")




# Load data ---------------------------------------------------------------

load(file.path(p1_R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(p1_R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData"))
load(file.path(p2_R_objects_directory, "01) Import and process sgRNA sequences.RData"))





# Define functions --------------------------------------------------------

PlasmidForRow <- function(row_index) {
  stopifnot("plasmid_list" %in% ls(envir = globalenv()))
  sg_seq <- sg_sequences_df[row_index, paste0("Sequence_sg", 1:4)]
  sg_lengths <- nchar(sg_seq)
  if (all(sg_lengths == 20)) {
    plasmid_name <- "pYJA5_with_4sg_20bp.gbk"
  } else if (all(sg_lengths == 19)) {
    plasmid_name <- "pYJA5_with_4sg_19bp.gbk"
  } else {
    stop("Unspecified empty plasmid input sequence!")
  }
  return(plasmid_name)
}




# Read in annotated plasmid sequences -------------------------------------

plasmid_files <- list.files(annotated_plasmid_directory)

plasmid_list <- sapply(plasmid_files, function(x) {
  ReadInPlasmid_gbk(file.path(annotated_plasmid_directory, x))
}, simplify = FALSE)





# Define the reference sequences ------------------------------------------

sg_sequences_df <- AddReferenceSequences(sg_sequences_df,
                                         tracRNAs_vec,
                                         promoters_vec,
                                         plasmid_string
                                         )



# Insert barcodes and guides into the annotated plasmids ------------------

plasmid_lines_list <- lapply(seq_len(nrow(sg_sequences_df)), function(x) {
  use_plasmid <- PlasmidForRow(x)
  use_lines <- plasmid_list[[use_plasmid]][["sequence"]]
  sg_seqs <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  well_number <- sg_sequences_df[["Well_number"]][[x]]
  replace_seqs <- c(column_bc_vec[[well_number]], sg_seqs, row_bc_vec[[well_number]])
  ReplaceNNNLines(use_lines, replace_seqs)
})




# Export the plain plasmid sequences --------------------------------------

well_names <- paste0("Well", formatC(sg_sequences_df[["Well_number"]], flag = "0", width = 3))
barcoded_plasmids <- paste0(column_bc_vec[sg_sequences_df[["Well_number"]]],
                            sg_sequences_df[["Whole_plasmid"]],
                            row_bc_vec[sg_sequences_df[["Well_number"]]]
                            )
sg_sequences_df[["Barcoded_plasmid"]] <- toupper(barcoded_plasmids)
fasta_titles <- paste0(">", well_names)
fasta_list <- lapply(seq_along(fasta_titles),
                     function(x) c(fasta_titles[[x]], barcoded_plasmids[[x]], "")
                     )

write.table(unlist(fasta_list),
            file = file.path(file_output_directory, "reference_sequences.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            )




# Export the annotated plasmids -------------------------------------------

for (i in seq_len(nrow(sg_sequences_df))) {
  well_number <- sg_sequences_df[["Well_number"]][[i]]
  file_name <- paste0(well_names[[i]], "_barcoded_cassette.gbk")
  use_plasmid <- PlasmidForRow(i)
  export_lines <- c(plasmid_list[[use_plasmid]][["preceding"]],
                    plasmid_lines_list[[i]],
                    plasmid_list[[use_plasmid]][["following"]]
                    )
  write.table(x         = export_lines,
              file      = file.path(reference_output_directory, file_name),
              col.names = FALSE,
              row.names = FALSE,
              quote     = FALSE
              )
}





# Save data ---------------------------------------------------------------

save(list = "sg_sequences_df",
     file = file.path(p2_R_objects_directory, "02) Create reference sequences for each well - sg_sequences_df.RData")
     )

save(list = "plasmid_list",
     file = file.path(p2_R_objects_directory, "02) Create reference sequences for each well - annotated plasmid.RData")
     )





