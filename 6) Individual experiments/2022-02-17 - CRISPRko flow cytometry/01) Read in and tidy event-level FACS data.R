## 2022-02-17


# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-02-17 - CRISPRko flow cytometry")
file_input_directory  <- file.path(file_directory, "2) Input")
FACS_data_directory   <- file.path(file_input_directory, "csv files")
R_objects_directory   <- file.path(file_directory, "3) R objects")



# Read in data ------------------------------------------------------------

flow_df_list <- lapply(list.files(FACS_data_directory),
                       function(x) data.frame("Sample_name" = x,
                                              read.csv(file.path(FACS_data_directory, x),
                                                       check.names = FALSE,
                                                       stringsAsFactors = FALSE
                                                       ),
                                              stringsAsFactors = FALSE,
                                              check.names = FALSE
                                              )
                       )

flow_df <- do.call(rbind.data.frame,
                   c(flow_df_list, stringsAsFactors = FALSE)
                   )



# Tidy data ---------------------------------------------------------------

flow_df[["Sample_name"]] <- sub("export_HEK293_", "", flow_df[["Sample_name"]], fixed = TRUE)
flow_df[["Sample_name"]] <- sub("_Singlets.csv", "", flow_df[["Sample_name"]], fixed = TRUE)

flow_df[["Sample_name"]] <- sub("CD47 CD109", "none", flow_df[["Sample_name"]], fixed = TRUE)

full_splits <- strsplit(flow_df[["Sample_name"]], "_", fixed = TRUE)
name_vec <- sapply(full_splits, "[[", 2)
name_splits <- strsplit(name_vec, " ", fixed = TRUE)

flow_df <- data.frame(
  flow_df["Sample_name"],
  "Sample_number"  = as.integer(sapply(full_splits, "[[", 3)),
  "Target_protein" = sapply(name_splits, "[[", 1),
  "gRNA_used"      = sapply(name_splits, "[[", 2),
  "Unstained"      = sapply(name_splits, "[[", 3) == "unstained",
  "Target_signal"  = NA,
  flow_df[, 2:ncol(flow_df)],
  stringsAsFactors = FALSE,
  check.names      = FALSE
)

flow_df[, "Target_signal"] <- ifelse(flow_df[, "Target_protein"] == "none",
                                     NA,
                                     ifelse(flow_df[, "Target_protein"] == "CD47",
                                            flow_df[, "APC-A"],
                                            flow_df[, "DsRed-A"]
                                            )
                                     )



# Re-order data -----------------------------------------------------------

sg_levels <- c("WT", "NT", paste0("sg", 1:4), "4sg")

new_order <- order(match(flow_df[["Target_protein"]], flow_df[["Target_protein"]]),
                   match(flow_df[["gRNA_used"]], sg_levels)
                   )
flow_df <- flow_df[new_order, ]
row.names(flow_df) <- NULL




# Create a data frame of just the sample info -----------------------------

samples_df <- flow_df[!(duplicated(flow_df[["Sample_name"]])), 1:5]
row.names(samples_df) <- NULL




# Save data ---------------------------------------------------------------

save(list = c("flow_df", "samples_df"),
     file = file.path(R_objects_directory, "01) Read in and tidy event-level FACS data.RData")
     )





