### 2022-09-13



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-09-02 - Illumina 4sg sequencing")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "02_reformat_CRISPRa_library.RData"))
load(file.path(rdata_dir, "06_extend_read_counts.RData"))
for (i in 1:8) {
  load(file.path(rdata_dir, paste0("03_look_up_sgRNAs_chunk", i, ".RData")))
}



# Create a combined matched_df --------------------------------------------

matched_df_chunks <- paste0("chunk", 1:8, "_matched_df")
matched_df <- do.call(rbind.data.frame,
                      c(lapply(matched_df_chunks, get),
                        stringsAsFactors = FALSE, make.row.names = FALSE
                        ))
rm(list = matched_df_chunks)




# Examine the total number of reads ---------------------------------------

num_reads <- nrow(matched_df)




# Explore sequences that feature mismatched bases -------------------------

### Examine reads with "N" base calls

have_N <- grepl("N", matched_df[, "Sequence_sg2"], fixed = TRUE) |
          grepl("N", matched_df[, "Sequence_sg3"], fixed = TRUE)
have_N_table <- table(have_N)
have_N_table / nrow(matched_df)

num_Ns_sg2 <- 20L - nchar(gsub("N", "", matched_df[, "Sequence_sg2"], fixed = TRUE))
num_Ns_sg3 <- 20L - nchar(gsub("N", "", matched_df[, "Sequence_sg3"], fixed = TRUE))
table(num_Ns_sg2)
table(num_Ns_sg3)



## Examine the proportion of reads featuring an sgRNA with a mismatched base

have_0MM <- (matched_df[, "Num_MM_sg2"] %in% 0) &
            (matched_df[, "Num_MM_sg3"] %in% 0)
have_match <- (matched_df[, "Num_MM_sg2"] %in% c(0, 1)) &
              (matched_df[, "Num_MM_sg3"] %in% c(0, 1))
have_1MM <- have_match & !(have_0MM)

have_0MM_table <- table(have_0MM)
have_1MM_table <- table(have_1MM)
have_0MM_table / nrow(matched_df)
have_1MM_table / nrow(matched_df)





# Check whether sgRNAs from the library are found in the data -------------

all_sg2_vec <- unique(ifelse(matched_df[, "Num_MM_sg2"] == 0,
                             matched_df[, "Sequence_sg2"],
                             matched_df[, "Correct_sgRNA_sg2"]
                             ))
all_sg3_vec <- unique(ifelse(matched_df[, "Num_MM_sg3"] == 0,
                             matched_df[, "Sequence_sg3"],
                             matched_df[, "Correct_sgRNA_sg3"]
                             ))

found_sg2 <- sg_sequences_df[, "Sequence_sg2"] %in% all_sg2_vec
found_sg3 <- sg_sequences_df[, "Sequence_sg3"] %in% all_sg3_vec

table(found_sg2, found_sg3)
table(found_sg2)
table(found_sg3)
table(found_sg2 | found_sg3)
table(found_sg2 & found_sg3)


combined_sg2orsg3 <- unique(c(all_sg2_vec, all_sg3_vec))
found_sg2_either <- sg_sequences_df[, "Sequence_sg2"] %in% combined_sg2orsg3
found_sg3_either <- sg_sequences_df[, "Sequence_sg3"] %in% combined_sg2orsg3

table(found_sg2_either, found_sg3_either)
table(found_sg2_either)
table(found_sg3_either)
table(found_sg2_either | found_sg3_either)
table(found_sg2_either & found_sg3_either)


table(!(found_sg2) & found_sg2_either)
table(!(found_sg3) & found_sg3_either)



# Examine sgRNAs that are not found in the sequencing data ----------------

table(matched_df[, "Num_MM_sg2"], useNA = "ifany")

are_found_sg2 <- matched_df[, "Num_MM_sg2"] %in% c(0, 1)
altered_sg2_reads <- matched_df[, "Sequence_sg2"][!(are_found_sg2)]
altered_sg2_reads <- altered_sg2_reads[!(altered_sg2_reads %in% sg_sequences_df[, "Sequence_sg2"])]
altered_sg2s_table <- table(altered_sg2_reads)

are_found_sg3 <- matched_df[, "Num_MM_sg3"] %in% c(0, 1)
altered_sg3_reads <- matched_df[, "Sequence_sg3"][!(are_found_sg3)]
altered_sg3_reads <- altered_sg3_reads[!(altered_sg3_reads %in% sg_sequences_df[, "Sequence_sg3"])]
altered_sg3s_table <- table(altered_sg3_reads)

table(altered_sg2s_table > 500)
head(sort(altered_sg2s_table, decreasing = TRUE), 10)

## Draw histograms

DrawHistogram(counts_df[, "Sum_MaySwitch_xMM"], num_breaks = 100,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )
DrawHistogram(altered_sg2s_table[altered_sg2s_table > 500], num_breaks = 150,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )
DrawHistogram(altered_sg2s_table, num_breaks = 150,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )

hist(counts_df[, "Sum_MaySwitch_xMM"], breaks = 600, col = "black", xlim = c(0, 20000))
hist(altered_sg2s_table[altered_sg2s_table > 500], breaks = 600, xlim = c(0, 20000), col = "black")
hist(altered_sg2s_table, breaks = 600, xlim = c(0, 20000), col = "black")




# Save data ---------------------------------------------------------------

save(list = "num_reads",
     file = file.path(rdata_dir, "07_examine_read_level_metrics_1_matched_df.RData")
     )



