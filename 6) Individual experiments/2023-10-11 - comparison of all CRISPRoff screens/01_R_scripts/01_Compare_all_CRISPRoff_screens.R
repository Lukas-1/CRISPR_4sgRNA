### 2023-10-11


# Load packages and source code -------------------------------------------

library("beeswarm")

CRISPR_root_directory    <- "~/CRISPR_4sgRNA"
experiments_directory    <- file.path(CRISPR_root_directory, "6) Individual experiments")
first_illumina_trial_dir <- file.path(experiments_directory, "2022-04-21 - Illumina paired-end 2sg - first trial")
R_functions_dir          <- file.path(first_illumina_trial_dir, "01_R_scripts", "R_functions")

source(file.path(R_functions_dir, "05_creating_figures_from_count_data.R")) # For VerticalAdjust and MakeEmptyPlot




# Define paths ------------------------------------------------------------

off_2sg_rdata_dir     <- file.path(experiments_directory,
                                  "2022-06-21 - Illumina paired-end 2sg - correct reference",
                                  "03_R_objects"
                                   )
off_4sg_rdata_dir     <- file.path(experiments_directory,
                                  "2022-09-02 - Illumina 4sg sequencing",
                                  "03_R_objects"
                                   )
pool_rdata_dir        <- file.path(experiments_directory,
                                  "2023-09-28 - prepooled vs postpooled - Illumina",
                                  "03_R_objects"
                                  )
project_dir           <- file.path(experiments_directory, "2023-10-11 - comparison of all CRISPRoff screens")
output_dir            <- file.path(project_dir, "02_output")



# Load data ---------------------------------------------------------------

load(file.path(off_2sg_rdata_dir, "09_create_figures_from_count_data.RData"))
load(file.path(off_2sg_rdata_dir, "10_recreate_figures_of_Nunez_et_al.RData"))
load(file.path(off_4sg_rdata_dir, "09_create_figures_from_count_data.RData"))
load(file.path(pool_rdata_dir, "10_create_figures_from_count_data.RData"))





# Define functions --------------------------------------------------------

BeePoints <- function(points_list,
                      side_gap          = 0.5,
                      left_gap          = side_gap,
                      right_gap         = side_gap,
                      y_lower_limit     = NULL,
                      y_upper_limit     = NULL,
                      y_axis_label      = NULL,
                      use_tcl           = 0.375,
                      y_axis_label_line = 1.8,
                      y_axis_mgp        = 0.55,
                      y_axis_n          = 5,
                      group_labels_y    = 1.25,
                      points_color      = brewer.pal(9, "Blues")[[9]],
                      point_cex         = 1,
                      mean_lwd          = 1
                      ) {


  ## Determine point positions
  num_groups <- length(points_list)
  group_positions <- seq_len(num_groups)
  groups_vec <- rep(group_positions, each = 2)
  group_limits <- c((min(group_positions) - left_gap) - (num_groups * 0.04),
                    (max(group_positions) + right_gap) + (num_groups * 0.04)
                    )

  ## Prepare the data axis
  if (is.null(y_lower_limit)) {
    y_lower_limit <- min(pretty(c(0, min(unlist(points_list)))))
  }
  if (is.null(y_upper_limit)) {
    y_upper_limit <- max(pretty(c(0, max(unlist(points_list)))))
  }
  numeric_limits <- c(y_lower_limit, y_upper_limit)

  ## Draw lines
  MakeEmptyPlot(x_limits = group_limits, y_limits = numeric_limits)
  use_width <- 0.25
  final_width <- use_width * ((max(group_positions) - min(group_positions)) / (num_groups - 1))
  SEMs <- vapply(points_list, function(x) sd(x) / length(x), numeric(1))
  means_vec <- vapply(points_list, mean, numeric(1))
  segments(x0  = group_positions - final_width,
           x1  = group_positions + final_width,
           y0  = means_vec,
           lwd = mean_lwd
           )


  ## Draw points
  beeswarm_df <- beeswarm(
    points_list,
    cex = point_cex,
    do.plot = FALSE,
    priority = "density",
    spacing = 1.01
  )
  points(x   = beeswarm_df[, "x"],
         y   = beeswarm_df[, "y"],
         cex = point_cex,
         pch = 16,
         col = points_color,
         xpd = NA
         )
  text(x      = beeswarm_df[, "x"],
       y      = beeswarm_df[, "y"],
       labels = unlist(lapply(points_list, function(x) seq_along(x))),
       cex    = 0.25,
       col    = brewer.pal(9, "Blues")[[3]],
       xpd    = NA
       )

  ## Draw the y axis
  tick_locations <- pretty(numeric_limits, n = y_axis_n)
  axis(2,
       at     = tick_locations,
       labels = format(tick_locations),
       las    = 2,
       mgp    = c(3, y_axis_mgp, 0),
       tcl    = -(use_tcl),
       lwd    = par("lwd")
       )
  if (!(is.null(y_axis_label))) {
    mtext(VerticalAdjust(y_axis_label),
          side = 2,
          line = y_axis_label_line,
          cex  = par("cex")
          )
  }

  text(x      = group_positions + diff(grconvertX(c(0, 0.25), from = "lines", to = "user")),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 0.5), from = "lines", to = "user")),
       labels = names(points_list),
       adj    = c(1, 0.5),
       srt    = 45,
       xpd    = NA
       )

  box(bty = "l")

}




# Prepare data ------------------------------------------------------------

SSMD_vec_list <- list(
  "Re-analysis"  = separation_original_mat["Robust SSMD", ],
  "CRISPRoff"    = separation_CRISPRoff_mat["Robust SSMD", ],
  "T.gonfio"     = separation_4sg_mat["Robust SSMD", ],
  "Pre-pooled"   = separation_prepooled_mat["Robust SSMD", ],
  "Post-pooled"  = separation_postpooled_mat["Robust SSMD", ]
)
SSMD_vec_list <- lapply(SSMD_vec_list, function(x) -(x))

AUC_vec_list <- list(
  "Re-analysis"  = AUC_original_vec[c("rep1", "rep2")],
  "CRISPRoff"    = AUC_CRISPRoff_vec[c("rep1", "rep2")],
  "T.gonfio"     = AUC_4sg_vec[c("rep1", "rep2")],
  "Pre-pooled"   = AUC_prepooled_vec[c("rep1", "rep2")],
  "Post-pooled"  = AUC_postpooled_vec[c("rep1", "rep2")]
)


for (create_PDF in c(FALSE, TRUE)) {
  if (create_PDF) {
    pdf(file.path(output_dir, "AUC and SSMD.pdf"),
        width = 4.5, height = 4.5
        )
  }

  BeePoints(AUC_vec_list, y_lower_limit = 0.8, y_axis_label = "AUC",
            y_axis_label_line = 2.4, point_cex = 0.9,
            mean_lwd = 1.25
            )
  title("Area under the ROC curve", cex.main = 1, font.main = 1)

  BeePoints(SSMD_vec_list, y_axis_label = "SSMD*",
            y_axis_label_line = 2.4, point_cex = 0.9,
            mean_lwd = 1.25
            )
  title("Robust estimate of the SSMD", cex.main = 1, font.main = 1)


  if (create_PDF) {
    dev.off()
  }
}












