### 7th August 2019 ###





# Import packages and source code -----------------------------------------

library("vcfR")





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR_4sgRNA"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
polymorphisms_directory <- file.path(CRISPR_input_directory, "Human genome", "Polymorphisms")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")

Kaviar_vcf_path <- file.path(polymorphisms_directory, "Kaviar-160204-Public-hg38-trimACgt3.vcf")




# Read in data ------------------------------------------------------------

# Downloaded from http://db.systemsbiology.net/kaviar/Kaviar.downloads.html
# (The file labelled "Only variants seen > 3 times (1.1 GB each file)", GRCh38)
# Direct link: http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg38-trimACgt3.vcf.tar
# accessed on 7 August 2019

vcfR_read <- read.vcfR(file = Kaviar_vcf_path)





# Process data ------------------------------------------------------------

info_splits <- strsplit(vcfR_read@fix[, "INFO"], ";", fixed = TRUE)

AF_raw_vec <- vapply(info_splits, function(x) grep("AF=", x, value = TRUE, fixed = TRUE), "")

Kaviar_common_raw_df <- data.frame(
  "Chromosome"   = vcfR_read@fix[, "CHROM"],
  "Position"     = as.integer(vcfR_read@fix[, "POS"]),
  "rsID"         = vcfR_read@fix[, "ID"],
  "Reference"    = vcfR_read@fix[, "REF"],
  "Alternative"  = vcfR_read@fix[, "ALT"],
  "AF_Kaviar"    = AF_raw_vec,
  stringsAsFactors = FALSE
)






# Save data ---------------------------------------------------------------

save(list = "Kaviar_common_raw_df",
     file = file.path(general_RData_directory, "04) Read in data on genetic polymorphisms from Kaviar.RData")
     )






