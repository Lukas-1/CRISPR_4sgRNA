### 18th November 2019 ###



# Set global variables ----------------------------------------------------

preferred_rsID_column   <- "all22_SNP_IDs_vcf"
# preferred_AF_sum_column <- "all23_SNP_AF_sum_Kaviar"
preferred_AF_max_column <- "all22_SNP_AF_max_Kaviar"
SNP_frequency_cutoff    <- 0.001 # 0.1%


SNP_column_names <- c("sgRNA_SNP_IDs_vcf",     "sgRNA_SNP_AFs_1kGenomes", "sgRNA_SNP_AF_max_1kGenomes", "sgRNA_SNP_AF_sum_1kGenomes",
                      "sgRNA_SNP_AFs_TOPMED",  "sgRNA_SNP_AF_max_TOPMED", "sgRNA_SNP_AF_sum_TOPMED",
                      "sgRNA_SNP_AFs_Kaviar",  "sgRNA_SNP_AF_max_Kaviar", "sgRNA_SNP_AF_sum_Kaviar",
                      "sgRNA_SNP_IDs_1kG_ph1", "sgRNA_SNP_AFs_1kG_ph1",   "sgRNA_SNP_AF_max_1kG_ph1", "sgRNA_SNP_AF_sum_1kG_ph1",
                      "sgRNA_SNP_IDs_1kG_ph3", "sgRNA_SNP_AFs_1kG_ph3",   "sgRNA_SNP_AF_max_1kG_ph3", "sgRNA_SNP_AF_sum_1kG_ph3",
                      "sgRNA_SNP_IDs_gnomAD",  "sgRNA_SNP_AFs_gnomAD",    "sgRNA_SNP_AF_max_gnomAD",  "sgRNA_SNP_AF_sum_gnomAD",

                      "PAM_SNP_IDs_vcf",     "PAM_SNP_AFs_1kGenomes", "PAM_SNP_AF_max_1kGenomes", "PAM_SNP_AF_sum_1kGenomes",
                      "PAM_SNP_AFs_TOPMED",  "PAM_SNP_AF_max_TOPMED", "PAM_SNP_AF_sum_TOPMED",
                      "PAM_SNP_AFs_Kaviar",  "PAM_SNP_AF_max_Kaviar", "PAM_SNP_AF_sum_Kaviar",
                      "PAM_SNP_IDs_1kG_ph1", "PAM_SNP_AFs_1kG_ph1",   "PAM_SNP_AF_max_1kG_ph1", "PAM_SNP_AF_sum_1kG_ph1",
                      "PAM_SNP_IDs_1kG_ph3", "PAM_SNP_AFs_1kG_ph3",   "PAM_SNP_AF_max_1kG_ph3", "PAM_SNP_AF_sum_1kG_ph3",
                      "PAM_SNP_IDs_gnomAD",  "PAM_SNP_AFs_gnomAD",    "PAM_SNP_AF_max_gnomAD",  "PAM_SNP_AF_sum_gnomAD",

                      "all23_SNP_IDs_vcf",     "all23_SNP_AFs_1kGenomes", "all23_SNP_AF_max_1kGenomes", "all23_SNP_AF_sum_1kGenomes",
                      "all23_SNP_AFs_TOPMED",  "all23_SNP_AF_max_TOPMED", "all23_SNP_AF_sum_TOPMED",
                      "all23_SNP_AFs_Kaviar",  "all23_SNP_AF_max_Kaviar", "all23_SNP_AF_sum_Kaviar",
                      "all23_SNP_IDs_1kG_ph1", "all23_SNP_AFs_1kG_ph1",   "all23_SNP_AF_max_1kG_ph1", "all23_SNP_AF_sum_1kG_ph1",
                      "all23_SNP_IDs_1kG_ph3", "all23_SNP_AFs_1kG_ph3",   "all23_SNP_AF_max_1kG_ph3", "all23_SNP_AF_sum_1kG_ph3",
                      "all23_SNP_IDs_gnomAD",  "all23_SNP_AFs_gnomAD",    "all23_SNP_AF_max_gnomAD",  "all23_SNP_AF_sum_gnomAD",

                      "all22_SNP_IDs_vcf",     "all22_SNP_AF_max_1kGenomes", "all22_SNP_AF_max_TOPMED", "all22_SNP_AF_max_Kaviar",
                      "all22_SNP_IDs_1kG_ph1", "all22_SNP_AF_max_1kG_ph1",
                      "all22_SNP_IDs_1kG_ph3", "all22_SNP_AF_max_1kG_ph3",
                      "all22_SNP_IDs_gnomAD",  "all22_SNP_AF_max_gnomAD"
                      )


