### 22nd September 2020 ###



# Define maps -------------------------------------------------------------

titles_list <- list(

  "Num_contaminated_reads" = "Percentage of reads that match sgRNAs from other wells (with 100% accuracy)",

  "Num_under_2kb"          = "Percentage of truncated reads (shorter than 2000 base pairs)",

  "Num_low_barcode_scores" = "Percentage of reads with low barcode scores",
  "Num_low_quality_scores" = "Percentage of reads with low average read quality scores",
  "Num_low_bc_or_qual"     = "Percentage of low-quality reads (barcode or read quality below the threshold)",
  "Mean_read_quality"      = "Mean read quality scores",

  "Longest_subsequence"    = "Length of the longest subsequence shared between any two gRNAs in that well",

  "Count_sg1_cr1"          = "Percentage of reads where sg1 (and its tracRNA) are 100% correct",
  "Count_sg2_cr2"          = "Percentage of reads where sg2 (and its tracRNA) are 100% correct",
  "Count_sg3_cr3"          = "Percentage of reads where sg3 (and its tracRNA) are 100% correct",
  "Count_sg4_cr4"          = "Percentage of reads where sg4 (and its tracRNA) are 100% correct",

  "Count_mean_sg1to4"      = "Percentage of individual gRNAs (+ tracRNA) that are 100% correct (mean of sg1-4)",

  "Count_at_least_1"       = "Percentage of reads with at least 1 sgRNA that is 100% correct",
  "Count_at_least_2"       = "Percentage of reads with at least 2 sgRNAs that are 100% correct",
  "Count_at_least_3"       = "Percentage of reads with at least 3 sgRNAs that are 100% correct",
  "Count_all_4"            = "Percentage of reads for which all 4 sgRNAs are 100% correct",

  "Count_all_4_promoters"  = "Percentage of reads for which all 4 sgRNAs (including the promoters) are 100% correct",
  "Count_whole_plasmid"    = "Percentage of reads where the entire plasmid sequence is 100% correct"
)




