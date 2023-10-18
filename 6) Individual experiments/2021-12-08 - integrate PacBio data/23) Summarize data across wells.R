### 4th January 2022 ###



# Import packages and source code -----------------------------------------

CRISPR_root_directory <- "~/CRISPR_4sgRNA"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory      <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
s2r1_directory        <- file.path(experiments_directory, "2021-04-03 - PacBio - first Sequel-II run")
R_functions_directory <- file.path(plate1_directory, "1) R functions")

source(file.path(R_functions_directory, "02) Analyzing reads.R"))
source(file.path(R_functions_directory, "08) Processing demultiplexed PacBio reads.R"))
source(file.path(R_functions_directory, "09) Producing heatmaps.R")) # For VerticalAdjust and related functions
source(file.path(R_functions_directory, "20) Summarizing data across wells.R"))



# Define folder paths -----------------------------------------------------

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
file_output_directory    <- file.path(s2rI_directory, "5) Output", "Figures", "Summaries across wells")
manuscript_directory     <- file.path(s2rI_directory, "5) Output", "Figures", "Manuscript")
PNGs_output_directory    <- file.path(s2rI_directory, "5) Output", "PNGs", "Summaries across wells")



# Load data ---------------------------------------------------------------

load(file.path(s2rI_R_objects_directory, "01) Process and export plate barcodes.RData"))
load(file.path(s2rI_R_objects_directory, "04) Create reference sequences for each well - sg_sequences_df.RData"))
load(file.path(s2rI_R_objects_directory, "09) Characterize contaminations (using aligned reads).RData"))
load(file.path(s2rI_R_objects_directory, "10) Identify and characterize deletions.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - ccs_df_lists.RData"))
load(file.path(s2rI_R_objects_directory, "11) Process demultiplexed PacBio reads - plates_analysis_list.RData"))




# Consolidate colony-picked control wells ---------------------------------

manhattan_dist_list <- MakeDistanceList(manhattan_distance = TRUE)

use_reads_df <- SampleControlReads(plates_analysis_list[["individual_reads_df"]],
                                   plates_df
                                   )
are_to_keep <- sg_sequences_df[, "Plate_number"] %in% use_reads_df[, "Plate_number"]
sg_sequences_df <- sg_sequences_df[are_to_keep, ]

use_analysis_list <- plates_analysis_list
use_analysis_list[["individual_reads_df"]] <- use_reads_df
use_analysis_list[["individual_reads_df"]][["Passed_filters"]] <- TRUE

ccs3_zmws <- GetCCS3_ZMWs(use_analysis_list[["individual_reads_df"]])
ccs7_zmws <- GetCCS7_ZMWs(use_analysis_list[["individual_reads_df"]])

use_IDs <- intersect(sg_sequences_df[["Combined_ID"]], use_reads_df[["Combined_ID"]])

ccs3_df_list <- SummarizeWells(use_analysis_list,
                               use_zmws           = ccs3_zmws,
                               ID_column          = "Combined_ID",
                               unique_IDs         = sg_sequences_df[, "Combined_ID"],
                               deletions_df       = deletions_df,
                               aligned_contam_df  = contam_df,
                               filter_cross_plate = TRUE
                               )

ccs7_df_list <- SummarizeWells(plates_analysis_list,
                               use_zmws           = ccs7_zmws,
                               ID_column          = "Combined_ID",
                               unique_IDs         = sg_sequences_df[, "Combined_ID"],
                               deletions_df       = deletions_df,
                               aligned_contam_df  = contam_df,
                               filter_cross_plate = TRUE
                               )




# Define plate selections -------------------------------------------------

are_CRISPRa  <- grepl("HA", plates_df[["Plate_name"]], fixed = TRUE)
are_CRISPRko <- grepl("HO", plates_df[["Plate_name"]], fixed = TRUE)

are_library <- are_CRISPRa | are_CRISPRko

plate_selection_list <- list(
  "All plates"    = plates_df[["Plate_name"]][are_library],
  "Colony-picked" = plates_df[["Plate_name"]][plates_df[["Colony_picked"]]],
  "CRISPRa"       = plates_df[["Plate_name"]][are_CRISPRa],
  "CRISPRko"      = plates_df[["Plate_name"]][are_CRISPRko]
)

plate_selection_titles_list <- list(
  "All plates" = "Both libraries",
  "CRISPRa"       = "CRISPRa library",
  "CRISPRko"      = "CRISPRko library"
)




# Create labels for individual plates -------------------------------------

plate_labels <- paste0("Plate #", plates_df[["Plate_number"]], " \u2013 ", plates_df[["Plate_name"]])
names(plate_labels) <- plates_df[["Plate_name"]]
plate_labels <- plate_labels[!(plates_df[["Colony_picked"]])]
are_single_plates <- c(rep(FALSE, length(plate_selection_titles_list)),
                       rep(TRUE, length(plate_labels))
                       )
plate_selection_titles_list <- c(plate_selection_titles_list, as.list(plate_labels))




# Prepare file name strings for plate selections --------------------------

plate_matches <- match(names(plate_selection_titles_list)[are_single_plates],
                       plates_df[["Plate_name"]]
                       )
plate_numbers <- plates_df[["Plate_number"]][plate_matches]
plate_selection_prefixes <- c(paste0("0", letters[seq_len(sum(!(are_single_plates)))]),
                              formatC(plate_numbers, width = 2, flag = "0")
                              )



# Produce example plots ---------------------------------------------------

use_df <- ccs3_df_list[["original_summary_df"]]

SummaryBoxPlot(use_df, "All plates")
SummaryBoxPlot(use_df, "CRISPRa")

SummaryBoxPlot(use_df, "CRISPRa", sina_plot = TRUE, draw_whiskers = TRUE,
               box_lwd = 2
               )



LollipopPlot(use_df,
             "All plates",
             c("Num_contaminated_reads",
               "Num_contaminated_reads_aligned",
               "Num_cross_plate_contaminated"
               ),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "All plates",
             c("Num_contaminated_reads"),
             use_y_limits = c(0, 0.05)
             )

LollipopPlot(use_df,
             "All plates",
             c("Num_reads_with_sgRNA_deletion",
               "Num_reads_with_deletions_exceeding_20bp",
               "Num_reads_with_deletions_spanning_tracrRNAs",
               "Num_reads_with_deletions_spanning_promoters"
               )
             )


ReadCountsBoxPlot(use_df, c("CRISPRa", "CRISPRko"))

SummaryStackedBars(use_df, "CRISPRa", consider_tracrRNAs = FALSE)
SummaryStackedBars(use_df, "CRISPRko", consider_tracrRNAs = FALSE)

SummaryStackedBars(use_df, "CRISPRa", consider_tracrRNAs = TRUE)
SummaryStackedBars(use_df, "CRISPRko", consider_tracrRNAs = TRUE)



# Create box plots for the thesis -----------------------------------------

manuscript_side_labels <- list(
  c("Colony-", "picked", "controls"),
  c("Wells in", "CRISPR", "library")
)

for (show_library in c("CRISPRa", "CRISPRo")) {

  for (include_promoters in c(FALSE, TRUE)) {

    use_side_labels <- manuscript_side_labels
    if (show_library == "CRISPRa") {
      use_side_labels[[2]][[2]] <- "T.gonfio"
    } else {
      use_side_labels[[2]][[2]] <- "T.spiezzo"
    }

    for (draw_figure in c("A", "C", "E", "F")) {

      figure_prefix <- draw_figure
      if (include_promoters) {
        figure_prefix <- file.path("Not used", draw_figure)
      }

      if (draw_figure == "A") {
        custom_labels <- list(
          "Count_at_least_1" = expression("" >= "1 gRNA",  "correct"),
          "Count_at_least_2" = expression("" >= "2", "correct"),
          "Count_at_least_3" = expression("" >= "3", "correct"),
          "Count_all_4"      = expression("All 4", "correct")
        )
        if (include_promoters) {
          selection_name <- "Promoter_num_guides"
          names(custom_labels) <- sub("Count_at", "Count_pr_at", names(custom_labels), fixed = TRUE)
        } else {
          selection_name <- "At_least_num_guides"
        }
      } else if (draw_figure == "C") {
        custom_labels <- lapply(1:4, function(x) c(paste0("sg", x), "correct"))
        if (include_promoters) {
          selection_name <- "Count_pr_sg_cr"
          names(custom_labels) <- paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
        } else {
          selection_name <- "Count_sg_cr"
          names(custom_labels) <- paste0("Count_sg", 1:4, "_cr", 1:4)
        }
      } else if (draw_figure == "E") {
        selection_name <- "Deletions"
        custom_labels <- list(
          "Num_reads_with_deletions_exceeding_20bp"     = expression("Any deletion", "(" >= "20 bp)"),
          "Num_reads_with_deletions_spanning_tracrRNAs" = c("spanning", "tracrRNAs"),
          "Num_reads_with_deletions_spanning_promoters" = c("spanning", "promoters")
        )
      } else if (draw_figure == "F") {
        selection_name <- "Contaminations"
        custom_labels <- list(
          "Num_contaminated_reads" = c("Well-to-well", "contamination")
        )
      }

      file_name <- paste0(figure_prefix, ") Box plot - ", selection_name, " - ", show_library, " library.emf")
      emf(file  = file.path(manuscript_directory,
                            if (include_promoters) "Fig. S4" else "Fig. 4",
                            "Individual plots",
                            file_name
                            ),
          width = if (draw_figure == "F") 2.7 else 3.4, height = 1.75,
          emfPlus = FALSE, coordDPI = 7000
          )
      old_par <- par(lwd = 0.8, cex = 0.7, mai = c(0.35, 0.5, 0.2, 0.75))
      SummaryBoxPlot(use_df,
                     plate_names        = plate_selection_list[[show_library]],
                     use_columns        = names(custom_labels),
                     custom_title       = if (show_library == "CRISPRa") "T.gonfio library" else "T.spiezzo library",
                     title_line         = 0.5,
                     title_cex          = 1,
                     bottom_labels_list = custom_labels,
                     y_label_line       = 2.65,
                     bottom_start_y     = 0.1,
                     draw_whiskers      = TRUE,
                     embed_PNG          = TRUE,
                     embed_PNG_res      = 900,
                     # PNG_adjust         = 0.00325,
                     set_mar            = FALSE,
                     points_centered    = TRUE,
                     side_labels_list   = use_side_labels,
                     side_legend_x      = 0.8,
                     side_legend_x_gap  = -0.2,
                     side_y_gap         = 2.25,
                     bottom_y_gap       = 1.8,
                     bottom_large_y     = 0.8,
                     median_lwd         = 2,
                     violin_wex         = if (draw_figure %in% c("E", "F")) 0.4 else 0.45, #0.4, 0.52
                     box_wex            = if (draw_figure == "E") 0.55 else if (draw_figure == "F") 1 else 0.4, #0.55
                     use_side_gap       = if (draw_figure == "E") 0.5 else if (draw_figure == "F") 0.8 else 0.35,
                     controls_x_gap     = if (draw_figure == "F") 0.45 else 0.4, # 0.35
                     sina_plot          = TRUE,
                     shift_left         = 0.03,
                     box_lwd            = 1.75,
                     grid_lwd           = 0.7
                     )
      par(old_par)
      dev.off()
    }
  }
}




# Create box plots for the manuscript -------------------------------------

manuscript_side_labels <- list(
  c("Colony-", "picked", "controls"),
  c("Wells in", "CRISPR", "library")
)

for (show_library in c("CRISPRa", "CRISPRo")) {

  for (include_promoters in c(FALSE, TRUE)) {

    use_side_labels <- manuscript_side_labels
    if (show_library == "CRISPRa") {
      use_side_labels[[2]][[2]] <- "T.gonfio"
    } else {
      use_side_labels[[2]][[2]] <- "T.spiezzo"
    }

    for (draw_figure in c("A", "C", "E", "F")) {

      figure_prefix <- draw_figure
      if (include_promoters) {
        figure_prefix <- file.path("Not used", draw_figure)
      }

      if (draw_figure == "A") {
        custom_labels <- list(
          "Count_at_least_1" = expression("" >= "1 gRNA",  "correct"),
          "Count_at_least_2" = expression("" >= "2", "correct"),
          "Count_at_least_3" = expression("" >= "3", "correct"),
          "Count_all_4"      = expression("All 4", "correct")
        )
        if (include_promoters) {
          selection_name <- "Promoter_num_guides"
          names(custom_labels) <- sub("Count_at", "Count_pr_at", names(custom_labels), fixed = TRUE)
        } else {
          selection_name <- "At_least_num_guides"
        }
      } else if (draw_figure == "C") {
        custom_labels <- lapply(1:4, function(x) c(paste0("sg", x), "correct"))
        if (include_promoters) {
          selection_name <- "Count_pr_sg_cr"
          names(custom_labels) <- paste0("Count_pr", 1:4, "_sg", 1:4, "_cr", 1:4)
        } else {
          selection_name <- "Count_sg_cr"
          names(custom_labels) <- paste0("Count_sg", 1:4, "_cr", 1:4)
        }
      } else if (draw_figure == "E") {
        selection_name <- "Deletions"
        custom_labels <- list(
          "Num_reads_with_deletions_exceeding_20bp"     = expression("Any deletion", "(" >= "20 bp)"),
          "Num_reads_with_deletions_spanning_tracrRNAs" = c("spanning", "tracrRNAs"),
          "Num_reads_with_deletions_spanning_promoters" = c("spanning", "promoters")
        )
      } else if (draw_figure == "F") {
        selection_name <- "Contaminations"
        custom_labels <- list(
          "Num_contaminated_reads" = c("Well-to-well", "contamination")
        )
      }

      file_name <- paste0(figure_prefix, ") Box plot - ", selection_name, " - ", show_library, " library.pdf")
      pdf(file  = file.path(manuscript_directory,
                            if (include_promoters) "Fig. S4" else "Fig. 4",
                            "Individual plots",
                            file_name
                            ),
          width = if (draw_figure == "F") 2.7 else 3.4, height = 1.75
          )
      old_par <- par(lwd = 0.8, cex = 0.7, mai = c(0.35, 0.5, 0.2, 0.75))
      SummaryBoxPlot(use_df,
                     plate_names        = plate_selection_list[[show_library]],
                     use_columns        = names(custom_labels),
                     custom_title       = if (show_library == "CRISPRa") "T.gonfio library" else "T.spiezzo library",
                     title_line         = 0.5,
                     title_cex          = 1,
                     bottom_labels_list = custom_labels,
                     y_label_line       = 2.65,
                     bottom_start_y     = 0.1,
                     draw_whiskers      = TRUE,
                     embed_PNG          = TRUE,
                     embed_PNG_res      = 900,
                     set_mar            = FALSE,
                     points_centered    = TRUE,
                     side_labels_list   = use_side_labels,
                     side_legend_x      = 0.8,
                     side_legend_x_gap  = -0.2,
                     side_y_gap         = 2.25,
                     bottom_y_gap       = 1.8,
                     bottom_large_y     = 0.8,
                     median_lwd         = 2,
                     violin_wex         = if (draw_figure %in% c("E", "F")) 0.4 else 0.52,
                     box_wex            = if (draw_figure == "E") 0.55 else if (draw_figure == "F") 1 else 0.6,
                     use_side_gap       = if (draw_figure == "E") 0.5 else if (draw_figure == "F") 0.8 else 0.3,
                     controls_x_gap     = if (draw_figure == "F") 0.425 else 0.35
                     )
      par(old_par)
      dev.off()
    }
  }
}




# Create stacked bar plots for the thesis ---------------------------------

for (show_library in c("CRISPRa", "CRISPRko")) {
  file_name <- paste0("E) Stacked bar plot - ", show_library, " library.emf")
  emf(file  = file.path(manuscript_directory, "Fig. 4", "Individual plots", file_name),
      width = 3.4, height = 1.75, emfPlus = FALSE, coordDPI = 7000
      )
  old_par <- par(lwd = 0.8, cex = 0.7, mai = c(0.35, 0.5, 0.2, 0.95))
  if (show_library == "CRISPRa") {
    library_name <- "T.gonfio"
  } else {
    library_name <- "T.spiezzo"
  }
  SummaryStackedBars(use_df,
                     show_library,
                     consider_tracrRNAs = TRUE,
                     top_title          = paste0(library_name, " library"),
                     top_title_line     = 0.3,
                     top_title_cex      = 1,
                     top_title_font     = 1,
                     set_mar            = FALSE,
                     x_labels_line      = 0.3,
                     y_label_line       = 2.65,
                     small_y_gap        = 1,
                     lines_x_start      = 0.7,
                     large_gap_ratio    = 1.45,
                     use_side_gap       = 0.6,
                     four_colors        = c("#F9F4EC", "#DB678B", "#a3c7eb", "#601A4A", NA),
                     grid_lwd           = 0.7
                     )
  par(old_par)
  dev.off()
}




# Create stacked bar plots for the manuscript -----------------------------

for (show_library in c("CRISPRa", "CRISPRko")) {
  file_name <- paste0("E) Stacked bar plot - ", show_library, " library.emf")
  pdf(file  = file.path(manuscript_directory, "Fig. 4", "Individual plots", file_name),
      width = 3.4, height = 1.75
      )
  old_par <- par(lwd = 0.8, cex = 0.7, mai = c(0.35, 0.5, 0.2, 0.95))
  if (show_library == "CRISPRa") {
    library_name <- "T.gonfio"
  } else {
    library_name <- "T.spiezzo"
  }
  SummaryStackedBars(use_df,
                     show_library,
                     consider_tracrRNAs = TRUE,
                     top_title          = paste0(library_name, " library"),
                     top_title_line     = 0.3,
                     top_title_cex      = 1,
                     top_title_font     = 1,
                     set_mar            = FALSE,
                     x_labels_line      = 0.3,
                     y_label_line       = 2.65,
                     small_y_gap        = 1,
                     lines_x_start      = 0.7,
                     large_gap_ratio    = 1.45,
                     use_side_gap       = 0.6
                     )
  par(old_par)
  dev.off()
}



# Export violin/box plots of total read counts ----------------------------

pdf(file.path(manuscript_directory,  "Fig. S4", "Individual plots",
              paste0("B) Read counts.pdf")
              ),
    width = 2.2, height = 1.7
    )
par(cex = 0.7, lwd = 0.8, mai = c(0.35, 0.47, 0.15, 0.5))
ReadCountsBoxPlot(use_df, c("CRISPRa", "CRISPRko"),
                  selection_labels = c("T.gonfio library", "T.spiezzo library"),
                  x_labels_line = 0.3,
                  y_label_line  = 2.4,
                  side_gap      = 0.5,
                  embed_PNG     = TRUE
                  )
dev.off()



# Export summary stacked bar plots ----------------------------------------

DrawAllSummaryBarPlots()



# Export lollipop plots and violin/box plots ------------------------------

DrawAllLollipopsAndViolins(export_PNGs = TRUE)




# Save data ---------------------------------------------------------------

save(list = c("plate_selection_list",
              "plate_selection_titles_list",
              "plate_selection_prefixes"
              ),
     file = file.path(s2rI_R_objects_directory, "23) Summarize data across wells - plate selections.RData")
     )



