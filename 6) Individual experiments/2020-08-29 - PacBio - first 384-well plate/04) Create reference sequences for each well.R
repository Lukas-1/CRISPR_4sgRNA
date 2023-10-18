### 9th October 2020 ###





# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "03) Creating reference plasmid sequences.R"))





# Define folder paths -----------------------------------------------------

file_input_directory        <- file.path(file_directory, "2) Input")
R_objects_directory         <- file.path(file_directory, "3) R objects")

metadata_directory          <- file.path(file_input_directory, "Metadata")
tracRNAs_file               <- file.path(metadata_directory, "tracRNAs_4sg.txt")
promoters_file              <- file.path(metadata_directory, "promoters_4sg.txt")
full_plasmid_file           <- file.path(metadata_directory, "reference with 4sg and barcode example.fa")

annotated_plasmid_directory <- file.path(file_input_directory, "Annotated plasmids")

file_output_directory       <- file.path(file_directory, "5) Output")
reference_output_directory  <- file.path(file_output_directory, "Reference sequences")




# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Process and export barcodes.RData"))
load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))




# Read in data ------------------------------------------------------------

tracRNAs_vec <- read.table(tracRNAs_file, header = FALSE,
                           stringsAsFactors = FALSE
                           )[[1]]
promoters_vec <- read.table(promoters_file, header = FALSE,
                            stringsAsFactors = FALSE
                            )[[1]]
plain_plasmid_vec <- read.table(full_plasmid_file, header = FALSE, skip = 1,
                                stringsAsFactors = FALSE
                                )[[1]]




# Process the plain plasmid sequence --------------------------------------

plasmid_string <- toupper(paste0(plain_plasmid_vec, collapse = ""))
plasmid_string <- substr(plasmid_string, 11, nchar(plasmid_string) - 10) # Remove barcode




# Define functions --------------------------------------------------------

PlasmidForWell <- function(well_number) {
  stopifnot("plasmid_list" %in% ls(envir = globalenv()))
  search_string <- paste0("_well", well_number, ".gbk")
  plasmid_name <- grep(search_string, names(plasmid_list), fixed = TRUE)
  if (length(plasmid_name) == 0) {
    plasmid_name <- "pYJA5_with_4sg_20bp.gbk"
  }
  return(plasmid_name)
}




# Read in annotated plasmid sequences -------------------------------------

plasmid_files <- list.files(annotated_plasmid_directory)

plasmid_list <- sapply(plasmid_files, function(x) {
  ReadInPlasmid_gbk(file.path(annotated_plasmid_directory, x))
}, simplify = FALSE)





# Insert barcodes and guides into the annotated plasmids ------------------

plasmid_lines_list <- lapply(seq_len(384), function(x) {
  use_plasmid <- PlasmidForWell(x)
  use_lines <- plasmid_list[[use_plasmid]][["sequence"]]
  sg_seqs <- as.character(sg_sequences_df[x, paste0("Sequence_sg", 1:4)])
  replace_seqs <- c(column_bc_vec[[x]], sg_seqs, row_bc_vec[[x]])
  ReplaceNNNLines(use_lines, replace_seqs)
})





# Define the reference sequences ------------------------------------------

sg_sequences_df <- AddReferenceSequences(sg_sequences_df,
                                         tracRNAs_vec,
                                         promoters_vec,
                                         plasmid_string
                                         )




# Export the plain plasmid sequences --------------------------------------

well_names <- paste0("Well", formatC(seq_len(384), flag = "0", width = 3))
barcoded_plasmids <- paste0(column_bc_vec, sg_sequences_df[["Whole_plasmid"]], row_bc_vec)
sg_sequences_df[["Barcoded_plasmid"]] <- toupper(barcoded_plasmids)
fasta_titles <- paste0(">", well_names)
fasta_list <- lapply(seq_len(384),
                     function(x) c(fasta_titles[[x]], barcoded_plasmids[[x]], "")
                     )

write.table(unlist(fasta_list),
            file = file.path(file_output_directory, "reference_sequences.fa"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,
            )




# Export the annotated plasmids -------------------------------------------

for (i in seq_len(384)) {
  file_name <- paste0(well_names[[i]], "_barcoded_cassette.gbk")
  use_plasmid <- PlasmidForWell(i)
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

save(list = c("plasmid_string", "tracRNAs_vec", "promoters_vec"),
     file = file.path(R_objects_directory, "04) Create reference sequences for each well - constant sequences.RData")
     )

save(list = "sg_sequences_df",
     file = file.path(R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData")
     )

save(list = "plasmid_list",
     file = file.path(R_objects_directory, "04) Create reference sequences for each well - annotated plasmid.RData")
     )





