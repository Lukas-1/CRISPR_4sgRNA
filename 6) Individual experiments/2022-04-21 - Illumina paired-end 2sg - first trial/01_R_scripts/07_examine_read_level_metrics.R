### 2022-04-09



# Load packages and source code -------------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_nanopore_dir <- file.path(experiments_directory, "2022-01-05 - first nanopore sequencing run")
source(file.path(first_nanopore_dir, "01_R_scripts", "1_R_functions", "02_creating_histograms.R"))



# Define paths ------------------------------------------------------------

project_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
rdata_dir <- file.path(project_dir, "03_R_objects")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "03_disambiguate_CRISPRoff_library.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk1.RData"))
load(file.path(rdata_dir, "04_look_up_sgRNAs_run2_chunk2.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__counts_df.RData"))
load(file.path(rdata_dir, "06_assign_sgRNAs_to_plasmids__lumi_df.RData"))



# Create a combined matched_df --------------------------------------------

matched_df <- rbind.data.frame(data.frame("Run" = 1L, run1_matched_df),
                               data.frame("Run" = 2L, run2_chunk1_matched_df),
                               data.frame("Run" = 2L, run2_chunk2_matched_df),
                               stringsAsFactors = FALSE,
                               make.row.names = FALSE
                               )
matched_df <- data.frame("Read_number" = seq_len(nrow(matched_df)),
                         matched_df, stringsAsFactors = FALSE
                         )
rm(list = c("run1_matched_df", "run2_chunk1_matched_df", "run2_chunk2_matched_df"))
gc()



# Explore sequences that feature mismatched bases -------------------------

### Examine reads with "N" base calls

have_N <- grepl("N", matched_df[, "Sequence_sg1"], fixed = TRUE) |
          grepl("N", matched_df[, "Sequence_sg2"], fixed = TRUE)
have_N_table <- table(have_N)
have_N_table / nrow(matched_df)

num_Ns_sg1 <- 19L - nchar(gsub("N", "", matched_df[, "Sequence_sg1"], fixed = TRUE))
num_Ns_sg2 <- 19L - nchar(gsub("N", "", matched_df[, "Sequence_sg2"], fixed = TRUE))
table(num_Ns_sg1)
table(num_Ns_sg2)



## Examine the proportion of reads featuring an sgRNA with a mismatched base

have_0MM <- (matched_df[, "Num_MM_sg1"] %in% 0) &
            (matched_df[, "Num_MM_sg2"] %in% 0)
have_match <- (matched_df[, "Num_MM_sg1"] %in% c(0, 1)) &
              (matched_df[, "Num_MM_sg2"] %in% c(0, 1))
have_1MM <- have_match & !(have_0MM)

have_0MM_table <- table(have_0MM)
have_1MM_table <- table(have_1MM)
have_0MM_table / nrow(matched_df)
have_1MM_table / nrow(matched_df)




# Create "sg_df" for convenience ------------------------------------------

sg_df <- data.frame(
  "sg1" = substr(toupper(CRISPRoff_df[, "protospacer_A"]), 2, 20),
  "sg2" = substr(toupper(CRISPRoff_df[, "protospacer_B"]), 2, 20),
  stringsAsFactors = FALSE
)



# Check whether sgRNAs from the library are found in the data -------------

all_sg1_vec <- unique(ifelse(matched_df[, "Num_MM_sg1"] == 0,
                             matched_df[, "Sequence_sg1"],
                             matched_df[, "Correct_sgRNA_sg1"]
                             ))
all_sg2_vec <- unique(ifelse(matched_df[, "Num_MM_sg2"] == 0,
                             matched_df[, "Sequence_sg2"],
                             matched_df[, "Correct_sgRNA_sg2"]
                             ))

found_sg1 <- sg_df[, "sg1"] %in% all_sg1_vec
found_sg2 <- sg_df[, "sg2"] %in% all_sg2_vec

table(found_sg1, found_sg2)
table(found_sg1)
table(found_sg2)
table(found_sg1 | found_sg2)
table(found_sg1 & found_sg2)


combined_sg1orsg2 <- unique(c(all_sg1_vec, all_sg2_vec))
found_sg1_either <- sg_df[, "sg1"] %in% combined_sg1orsg2
found_sg2_either <- sg_df[, "sg2"] %in% combined_sg1orsg2

table(found_sg1_either, found_sg2_either)
table(found_sg1_either)
table(found_sg2_either)
table(found_sg1_either | found_sg2_either)
table(found_sg1_either & found_sg2_either)


table(!(found_sg1) & found_sg1_either)
table(!(found_sg2) & found_sg2_either)



# Check the GC content of missing sgRNAs ----------------------------------

t.test(CRISPRoff_df[, "Num_GC_sg1"][found_sg1_either],
       CRISPRoff_df[, "Num_GC_sg1"][!(found_sg1_either)]
       )

boxplot(CRISPRoff_df[, "Num_GC_sg1"][found_sg1_either],
        CRISPRoff_df[, "Num_GC_sg1"][!(found_sg1_either)],
        ylim = c(0, 20), yaxs = "i", axes = FALSE, pch = 16,
        col = adjustcolor(brewer.pal(9, "Blues")[[4]]), cex = 0.8,
        boxwex = 0.5
        )
box(bty = "l")
axis(2, las = 1, tcl = -0.35, mgp = c(3, 0.5, 0))
mtext("Number of GC bases", side = 2, line = 2)
mtext(c("sgRNAs found\nin data", "sgRNAs absent\nfrom data"),
      side = 1, at = 1:2, line = 1.5
      )
title("Does GC content affect recovery of sgRNAs?", cex.main = 1)



# Examine the incidence of template switching -----------------------------

## There were 138'385'789 paired-end reads in total,
## combined from both sequencing runs (37.6% from the first run, the rest from
## the second.)

## 90.4% of paired-end reads contained at least one sgRNA sequence
## from the library (for either sg1 or sg2).
## One mismatch within the sgRNA sequence was allowed, as long as there was an
## unambiguous match to an sgRNA in the library.
nrow(lumi_df) / nrow(matched_df)


## 76.2% of paired-end reads contained two sgRNA sequences from the library.
## (In relative terms, 84.2% of just those reads that contained one sgRNA from
## the library also featured a second sgRNA from the library.)
sum(lumi_df[, "Num_matched_sgRNAs"] == 2) / nrow(matched_df)
sum(lumi_df[, "Num_matched_sgRNAs"] == 2) / nrow(lumi_df)


## 24.4% of reads (with two sgRNA sequences that are found in the library)
## had a template switch, i.e., the two sgRNAs belong to different plasmids.
## (For the purpose of this statistic, I only considered reads for which both
## sgRNAs were found in the library.)
have_2sg <- lumi_df[, "Num_matched_sgRNAs"] == 2
sum(lumi_df[, "Has_template_switch"][have_2sg]) / sum(have_2sg)



# Examine sgRNAs that are not found in the sequencing data ----------------

table(matched_df[, "Num_MM_sg1"], useNA = "ifany")

are_found_sg1 <- matched_df[, "Num_MM_sg1"] %in% c(0, 1)
altered_sg1_reads <- matched_df[, "Sequence_sg1"][!(are_found_sg1)]
altered_sg1_reads <- altered_sg1_reads[!(altered_sg1_reads %in% sg_df[, "sg1"])]
altered_sg1s_table <- table(altered_sg1_reads)

are_found_sg2 <- matched_df[, "Num_MM_sg2"] %in% c(0, 1)
altered_sg2_reads <- matched_df[, "Sequence_sg2"][!(are_found_sg2)]
altered_sg2_reads <- altered_sg2_reads[!(altered_sg2_reads %in% sg_df[, "sg2"])]
altered_sg2s_table <- table(altered_sg2_reads)

table(altered_sg1s_table > 500)
head(sort(altered_sg1s_table, decreasing = TRUE), 10)

## Draw histograms

DrawHistogram(counts_df[, "Sum_MaySwitch_xMM"], num_breaks = 100,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )
DrawHistogram(altered_sg1s_table[altered_sg1s_table > 500], num_breaks = 150,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )
DrawHistogram(altered_sg1s_table, num_breaks = 150,
              truncation_limit = 20000, x_axis_upper_limit = 20000
              )

hist(counts_df[, "Sum_MaySwitch_xMM"], breaks = 600, col = "black", xlim = c(0, 20000))
hist(altered_sg1s_table[altered_sg1s_table > 500], breaks = 600, xlim = c(0, 20000), col = "black")
hist(altered_sg1s_table, breaks = 600, xlim = c(0, 20000), col = "black")



