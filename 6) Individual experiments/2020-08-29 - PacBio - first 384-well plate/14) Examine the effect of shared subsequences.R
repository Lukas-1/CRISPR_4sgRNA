### 22nd September 2020 ###



# Import packages and source code -----------------------------------------


CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "10) Examining the effect of shared subsequences.R"))




# Define folder paths -----------------------------------------------------

file_output_directory  <- file.path(file_directory, "5) Output")
R_objects_directory    <- file.path(file_directory, "3) R objects")
plots_output_directory <- file.path(file_output_directory, "Figures")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "10) Process demultiplexed PacBio reads.RData"))






# Draw plots --------------------------------------------------------------

old_mar <- par(mar = use_mar)

for (var in use_vars) {
  PlotBySharedSubsequence(sl7_ccs3_df_list[["original_summary_df"]], var)
}




# Start loop --------------------------------------------------------------

for (use_filtered in c(FALSE, TRUE)) {
  for (smrtlink_version in c(7, 9)) {

    file_prefix <- paste0("Shared subsequences - SmrtLink ", smrtlink_version, " - ")

    if (use_filtered) {
      file_postfix <- " - filtered"
      df_name <- "filtered_summary_df"
    } else {
      file_postfix <- " - original"
      df_name <- "original_summary_df"
    }
    file_postfix <- paste0(file_postfix, ".pdf")

    if (smrtlink_version == 7) {
      ccs3_df_list <- sl7_ccs3_df_list
      ccs5_df_list <- sl7_ccs5_df_list
    } else if (smrtlink_version == 9) {
      ccs3_df_list <- sl9_ccs3_df_list
      ccs5_df_list <- sl9_ccs5_df_list
    }



    # Export plots ------------------------------------------------------------

    pdf(file = file.path(plots_output_directory,
                         paste0("SmrtLink ", smrtlink_version),
                         "Shared subsequences",
                         paste0(file_prefix, "CCS3", file_postfix)
                         ),
        width = pdf_width,
        height = pdf_height
        )
    par("mar" = use_mar)
    for (var in use_vars) {
      PlotBySharedSubsequence(ccs3_df_list[["original_summary_df"]], var)
    }
    dev.off()


    pdf(file = file.path(plots_output_directory,
                         paste0("SmrtLink ", smrtlink_version),
                         "Shared subsequences",
                         paste0(file_prefix, "CCS5", file_postfix)
                         ),
        width = pdf_width,
        height = pdf_height
        )
    par("mar" = use_mar)
    for (var in use_vars) {
      PlotBySharedSubsequence(ccs5_df_list[["original_summary_df"]], var)
    }
    dev.off()



# End loop ----------------------------------------------------------------

  }
}


PlotBySharedSubsequence(sl7_ccs5_df_list[["original_summary_df"]], "Count_mean_sg1to4")







