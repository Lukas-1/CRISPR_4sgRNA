### 14th September 2021 ###


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")

file_directory        <- file.path(experiments_directory, "2021-09-13 - visualize flow cytometry experiments")
file_input_directory  <- file.path(file_directory, "2) Input")
FACS_data_directory   <- file.path(file_input_directory, "Data for BFP+ cells/Gating 2")
R_objects_directory   <- file.path(file_directory, "3) R objects")



# Read in data ------------------------------------------------------------

BFP_df_list <- lapply(list.files(FACS_data_directory),
                      function(x) data.frame("Sample_name" = x,
                                             read.csv(file.path(FACS_data_directory, x),
                                                      check.names = FALSE,
                                                      stringsAsFactors = FALSE
                                                      ),
                                             stringsAsFactors = FALSE,
                                             check.names = FALSE
                                             )
                      )

BFP_df <- do.call(rbind.data.frame,
                  c(BFP_df_list, stringsAsFactors = FALSE)
                  )



# Tidy data ---------------------------------------------------------------

BFP_df[["Sample_name"]] <- sub("HEK293_", "", BFP_df[["Sample_name"]], fixed = TRUE)
BFP_df[["Sample_name"]] <- sub(".exported.FCS3.csv", "", BFP_df[["Sample_name"]], fixed = TRUE)

full_splits <- strsplit(BFP_df[["Sample_name"]], "_", fixed = TRUE)
first_half_vec <- sapply(full_splits, "[[", 1)
first_half_splits <- strsplit(first_half_vec, "-stained-", fixed = TRUE)

BFP_df <- data.frame(
  BFP_df["Sample_name"],
  "Sample_number"  = sapply(full_splits, "[[", 2),
  "Target_protein" = sapply(first_half_splits, "[[", 1),
  "gRNA_used"      = sapply(first_half_splits, "[[", 2),
  BFP_df[, 2:ncol(BFP_df)],
  stringsAsFactors = FALSE,
  check.names      = FALSE
)



# Re-order data -----------------------------------------------------------

sg_levels <- c("WT", "NT", paste0("sg", 1:4), "4sg")

new_order <- order(match(BFP_df[["Target_protein"]], BFP_df[["Target_protein"]]),
                   match(BFP_df[["gRNA_used"]], sg_levels)
                   )
BFP_df <- BFP_df[new_order, ]
row.names(BFP_df) <- NULL




# Create a data frame of just the sample info -----------------------------

samples_df <- BFP_df[!(duplicated(BFP_df[["Sample_name"]])), 1:4]
row.names(samples_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = c("BFP_df", "samples_df"),
     file = file.path(R_objects_directory, "01) Read in and tidy event-level FACS data.RData")
     )





