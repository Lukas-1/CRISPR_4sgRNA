### 22nd September 2020 ###



# Import packages and source code -----------------------------------------


CRISPR_root_directory <- "~/CRISPR"
file_directory        <- file.path(CRISPR_root_directory, "6) Individual experiments/2020-08-29 - PacBio - first 384-well plate")
R_functions_directory <- file.path(file_directory, "1) R functions")

source(file.path(R_functions_directory, "01) Define titles and labels.R"))
source(file.path(R_functions_directory, "10) Examining the effect of shared subsequences.R"))
source(file.path(R_functions_directory, "09) Producing heatmaps.R"))




# Define folder paths -----------------------------------------------------

file_output_directory  <- file.path(file_directory, "5) Output")
R_objects_directory    <- file.path(file_directory, "3) R objects")
plots_output_directory <- file.path(file_output_directory, "Figures")





# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(R_objects_directory, "11) Process demultiplexed PacBio reads.RData"))





# Draw plots --------------------------------------------------------------

old_mar <- par(mar = use_mar)

for (var in use_vars) {
  PlotBySharedSubsequence(sl7_ccs3_df_list[["original_summary_df"]], var)
}




# Start loop --------------------------------------------------------------

for (filter_category in c("original", "filtered reads", "filtered gRNAs")) {
  for (smrtlink_version in c(7, 9)) {

    file_prefix <- paste0("Shared subsequences - SmrtLink ", smrtlink_version, " - ")

    if (filter_category == "original") {
      file_postfix <- " - original"
      df_name <- "original_summary_df"
    } else if (filter_category == "filtered reads") {
      file_postfix <- " - filtered"
      df_name <- "filtered_summary_df"
    } else if (filter_category == "filtered gRNAs") {
      file_postfix <- " - filtered gRNAs"
      df_name <- "filtered_gRNAs_df"
    }
    file_postfix <- paste0(file_postfix, ".pdf")

    if (smrtlink_version == 7) {
      ccs3_df_list <- sl7_ccs3_df_list
      ccs5_df_list <- sl7_ccs5_df_list
      ccs7_df_list <- sl7_ccs7_df_list
    } else if (smrtlink_version == 9) {
      ccs3_df_list <- sl9_ccs3_df_list
      ccs5_df_list <- sl9_ccs5_df_list
      ccs7_df_list <- sl9_ccs7_df_list
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
      PlotBySharedSubsequence(ccs3_df_list[[df_name]], var)
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
      PlotBySharedSubsequence(ccs5_df_list[[df_name]], var)
    }
    dev.off()


    pdf(file = file.path(plots_output_directory,
                         paste0("SmrtLink ", smrtlink_version),
                         "Shared subsequences",
                         paste0(file_prefix, "CCS7", file_postfix)
                         ),
        width = pdf_width,
        height = pdf_height
        )
    par("mar" = use_mar)
    for (var in use_vars) {
      PlotBySharedSubsequence(ccs7_df_list[[df_name]], var)
    }
    dev.off()



# End loop ----------------------------------------------------------------

  }
}


manuscript_vars <- c("Num_under_2kb", "Count_all_4")
manuscript_directory <- file.path(file_output_directory, "Figures", "For the manuscript")

axis_labels <- c(
  "Num_under_2kb" = "Truncated reads (< 2000 bp)",
  "Count_all_4"   = "All 4 sgRNAs + tracRNAs correct"
)


for (var in manuscript_vars) {
  pdf(file.path(manuscript_directory, paste0("Shared sub-sequences - ", var, ".pdf")),
      width = 3.9, height = 2.6
      )
  par(cex = 0.7, lwd = 0.8, mai = rep(0.5, 4))
  PlotBySharedSubsequence(sl7_ccs7_df_list[["filtered_summary_df"]],
                          var,
                          grid_light_gray = "gray85",
                          grid_dark_gray  = "gray95",
                          use_spacing     = 0.4,
                          use_boxwex      = 0.75,
                          grid_lwd        = 0.75,
                          y_axis_label    = axis_labels[[var]],
                          use_title       = "",
                          corr_line       = 1.3,
                          bold_corr       = FALSE,
                          x_axis_label    = "Longest shared subsequence (base pairs)"
                          )
  dev.off()
}




PlotBySharedSubsequence(sl7_ccs7_df_list[["filtered_summary_df"]], "Count_all_4")














