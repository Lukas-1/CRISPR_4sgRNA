### 18th November 2019 ###



# Set global variables ----------------------------------------------------

preferred_rsID_column   <- "all22_SNP_IDs_vcf"
preferred_AF_max_column <- "all22_SNP_AF_max_Kaviar"
SNP_frequency_cutoff    <- 0.001 # ==> 0.1%


SNP_column_stems <- c(
    "SNP_IDs_vcf",     "SNP_AFs_1kGenomes", "SNP_AF_max_1kGenomes", "SNP_AF_sum_1kGenomes",
                       "SNP_AFs_TOPMED",    "SNP_AF_max_TOPMED",    "SNP_AF_sum_TOPMED",
                       "SNP_AFs_Kaviar",    "SNP_AF_max_Kaviar",    "SNP_AF_sum_Kaviar",

    "SNP_IDs_1kG_ph1", "SNP_AFs_1kG_ph1",   "SNP_AF_max_1kG_ph1",   "SNP_AF_sum_1kG_ph1",
    "SNP_IDs_1kG_ph3", "SNP_AFs_1kG_ph3",   "SNP_AF_max_1kG_ph3",   "SNP_AF_sum_1kG_ph3",
    "SNP_IDs_gnomAD",  "SNP_AFs_gnomAD",    "SNP_AF_max_gnomAD",    "SNP_AF_sum_gnomAD"
)

SNP_column_names_list <- sapply(c("sgRNA", "PAM", "all23", "all22"), function(x) paste0(x, "_", SNP_column_stems), simplify = FALSE)
SNP_column_names_list[["all22"]] <- grep("_SNP_AF(s|_sum)_", SNP_column_names_list[["all22"]], value = TRUE, invert = TRUE)
SNP_column_names <- unlist(SNP_column_names_list, use.names = FALSE)
