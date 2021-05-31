### 22nd September 2020 ###



# Define maps -------------------------------------------------------------

titles_list <- list(

  "Num_contaminated_reads" = "Percentage of reads that match gRNAs from other wells (with 100% accuracy)",

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

  "Count_at_least_1"       = "Percentage of reads with at least 1 gRNA that is 100% correct",
  "Count_at_least_2"       = "Percentage of reads with at least 2 gRNAs that are 100% correct",
  "Count_at_least_3"       = "Percentage of reads with at least 3 gRNAs that are 100% correct",
  "Count_all_4"            = "Percentage of reads for which all 4 gRNAs are 100% correct",

  "Count_all_4_promoters"  = "Percentage of reads for which all 4 gRNAs (including the promoters) are 100% correct",
  "Count_whole_plasmid"    = "Percentage of reads where the entire plasmid sequence is 100% correct",

  "Binary_num_contaminated_reads" = "20% or more of the reads match gRNAs from other wells (with 100% accuracy)",

  "Binary_all_four_guides"        = "For all 4 gRNAs, the majority of reads contain the correct sequence",
  "Binary_count_mean_sg1to4"      = "The mean percentage of correct gRNAs (for sg1-4) is 50% or higher",

  "Binary_count_at_least_1"       = "In the majority of reads, at least 1 gRNA is 100% correct",
  "Binary_count_at_least_2"       = "In the majority of reads, at least 2 gRNAs are 100% correct",
  "Binary_count_at_least_3"       = "In the majority of reads, at least 3 gRNAs are 100% correct",
  "Binary_count_all_4"            = "In the majority of reads, all 4 gRNAs are 100% correct",

  "Binary_count_all_4_promoters"  = "In the majority of reads, all 4 gRNAs (including the promoters) are 100% correct",
  "Binary_count_whole_plasmid"    = "In the majority of reads, the entire plasmid sequence is 100% correct",

  "Binary_count_sg1_cr1"          = "In the majority of reads, sg1 (and its tracRNA) are 100% correct",
  "Binary_count_sg2_cr2"          = "In the majority of reads, sg2 (and its tracRNA) are 100% correct",
  "Binary_count_sg3_cr3"          = "In the majority of reads, sg3 (and its tracRNA) are 100% correct",
  "Binary_count_sg4_cr4"          = "In the majority of reads, sg4 (and its tracRNA) are 100% correct"


)




