### 7th August 2019 ###



# Import packages and source code -----------------------------------------

library("vcfR")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
polymorphisms_directory <- file.path(CRISPR_input_directory, "Human genome", "Polymorphisms")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")

NCBI_vcf_path <- file.path(polymorphisms_directory, "common_all_20180418.vcf")




# Read in data ------------------------------------------------------------

# Downloaded from ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
# on 5 August 2019
vcfR_read <- read.vcfR(file = NCBI_vcf_path)




# Process data ------------------------------------------------------------

info_splits <- strsplit(vcfR_read@fix[, "INFO"], ";", fixed = TRUE)

CAF_raw_vec <- vapply(info_splits, function(x) grep("CAF=", x, value = TRUE, fixed = TRUE), "")
TOPMED_raw_vec <- vapply(info_splits, function(x) {
  result_string <- grep("TOPMED=", x, value = TRUE, fixed = TRUE)
  if (length(result_string) == 0) {
    return(NA_character_)
  } else {
    return(result_string)
  }
}, "")


dbSNP_common_raw_df <- data.frame(
  "Chromosome"   = vcfR_read@fix[, "CHROM"],
  "Position"     = as.integer(vcfR_read@fix[, "POS"]),
  "rsID"         = vcfR_read@fix[, "ID"],
  "Reference"    = vcfR_read@fix[, "REF"],
  "Alternative"  = vcfR_read@fix[, "ALT"],
  "AF_1kGenomes" = CAF_raw_vec,
  "AF_TOPMED"    = TOPMED_raw_vec,
  stringsAsFactors = FALSE
)




# Save data ---------------------------------------------------------------

save(list = "dbSNP_common_raw_df",
     file = file.path(general_RData_directory, "03) Read in data on genetic polymorphisms from dbSNP.RData")
     )



