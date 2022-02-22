## 2022-02-17


# Import packages and source code -----------------------------------------

library("vioplot")
library("RColorBrewer")
# library("flowWorkspace")



# Define folder paths -----------------------------------------------------

CRISPR_root_directory <- "~/CRISPR"
experiments_directory <- file.path(CRISPR_root_directory, "6) Individual experiments")
file_directory        <- file.path(experiments_directory, "2022-02-17 - CRISPRko flow cytometry")
file_input_directory  <- file.path(file_directory, "2) Input")
R_objects_directory   <- file.path(file_directory, "3) R objects")
file_input_directory  <- file.path(file_directory, "4) Output")



# Load data ---------------------------------------------------------------

load(file.path(R_objects_directory, "01) Read in and tidy event-level FACS data.RData"))



# Define functions --------------------------------------------------------

PrepareFACSData <- function(input_df,
                            use_target,
                            target_column = "Target_protein",
                            group_column  = "gRNA_used",
                            data_column   = "Target_signal",
                            exclude_WT    = TRUE,
                            equal_numbers = TRUE
                            ) {

  are_this_target <- input_df[, target_column] %in% use_target
  if (exclude_WT) {
    are_WT <- input_df[, group_column] %in% "WT"
    are_selected <- are_this_target & !(are_WT)
  } else {
    are_selected <- are_this_target
  }

  groups_vec <- input_df[are_selected, group_column]
  groups_fac <- factor(groups_vec, levels = unique(groups_vec))

  # BiexpFunction <- flowWorkspace::flowjo_biexp()

  results_df <- data.frame(
    input_df[are_selected, target_column, drop = FALSE],
    "Group"          = groups_fac,
    "Original_data"  = input_df[[data_column]][are_selected],
    # "Biexp_data"     = BiexpFunction(input_df[[data_column]][are_selected]),
    stringsAsFactors = FALSE
  )

  if (equal_numbers) {
    num_to_select <- min(tabulate(groups_fac))
    are_chosen <- rep(FALSE, length(groups_fac))
    set.seed(1)
    for (i in seq_len(nlevels(groups_fac))) {
      are_this_level <- as.integer(groups_fac) == i
      chosen_indices <- sample(seq_len(sum(are_this_level)), num_to_select)
      are_chosen[are_this_level][chosen_indices] <- TRUE
    }
    results_df <- results_df[are_chosen, ]
    row.names(results_df) <- NULL
  }
  return(results_df)
}



FACSViolinBox <- function(input_df,
                          use_target,
                          use_brewer = "Blues",
                          exclude_WT = TRUE
                          ) {

  sub_df <- PrepareFACSData(input_df, use_target, exclude_WT = exclude_WT)

  DrawViolinBox(sub_df[, "Original_data"],
                sub_df[, "Group"],
                use_brewer = use_brewer,
                use_title = use_target
                )

  return(invisible(NULL))
}



DrawViolinBox <- function(numeric_vec,
                          groups_fac,
                          use_brewer = "Blues",
                          use_title = NULL
                          ) {


  ## Determine group positions
  num_groups <- nlevels(groups_fac)
  group_positions <- seq_len(num_groups)
  group_limits <- c((min(group_positions) - 0.3) - (num_groups * 0.04),
                     max(group_positions) + 0.3  + (num_groups * 0.04)
                    )


  ## Prepare the data axis
  y_range <- range(numeric_vec)
  y_span <- y_range[[2]] - y_range[[1]]
  y_space_fraction <- 0.02
  final_y_limits <- c(y_range[[1]] - (y_span * y_space_fraction),
                      y_range[[2]] + (y_span * y_space_fraction)
                      )


  ## Set up the plot canvas
  plot(1,
       xlim = group_limits,
       ylim = final_y_limits,
       xaxs = "i",
       yaxs = "i",
       type = "n",
       axes = FALSE,
       ann  = FALSE
       )

  abline(h = median(numeric_vec[groups_fac == "NT"]),
         col = "gray85"
         )


  ## Draw the title
  if (!(is.null(use_title))) {
    title(use_title, cex.main = 1.1, font.main = 1, line = 2)
  }


  ## Draw the violin
  use_wex <- 0.8
  vioplot(numeric_vec ~ groups_fac,
          at       = group_positions,
          pchMed   = NA,
          drawRect = FALSE,
          col      = colorRampPalette(brewer.pal(9, use_brewer)[c(3, 4)])(3)[[2]],
          border   = NA,
          wex      = use_wex,
          add      = TRUE,
          axes     = FALSE
          )


  ## Draw the jittered points
  jittered_vec  <- group_positions[as.integer(groups_fac)] +
                   rnorm(n = length(groups_fac), mean = 0, sd = 0.05)
  points_alpha <- 0.3
  alpha_hex <- substr(rgb(1, 1, 1, points_alpha), 8, 9)
  points(x   = jittered_vec,
         y   = numeric_vec,
         cex = 0.3,
         col = paste0(brewer.pal(9, use_brewer)[[8]], alpha_hex),
         pch = 16
         )


  ## Draw the superimposed box plots
  box_alpha <- 0.5
  alpha_hex <- substr(rgb(1, 1, 1, box_alpha), 8, 9)
  boxplot(numeric_vec ~ groups_fac,
          at         = group_positions,
          boxwex     = use_wex * 0.4,
          outline    = FALSE,
          names      = rep.int("", length(group_positions)),
          whisklty   = "blank",
          staplewex  = 0,
          whisklwd   = 0,
          staplelty  = 0,
          medlwd     = par("lwd") * 3,
          col        = paste0(brewer.pal(9, use_brewer)[[1]], alpha_hex),
          border     = brewer.pal(9, use_brewer)[[8]],
          add        = TRUE,
          axes       = FALSE,
          lwd        = 1
          )


  ## Draw the axis and axis labels
  axis(2,
       mgp      = c(3, 0.45, 0),
       gap.axis = 0,
       tcl      = -0.3,
       las      = 1
       )
  box(bty = "l")

  y_axis_label <- "Fluorescence (biexponential transform)"
  mtext(text = y_axis_label, side = 2, line = 3)

  mtext(text = levels(groups_fac),
        at   = group_positions,
        side = 1,
        line = 0.6
        )

  return(invisible(NULL))
}




# Export plots ------------------------------------------------------------

FACSViolinBox(flow_df, "CD47", "Blues", exclude_WT = FALSE)
FACSViolinBox(flow_df, "CD109", "Purples", exclude_WT = FALSE)


pdf(file.path(file_output_directory, "FACS data - BFP+ population.pdf"),
    width = 6, height = 5.8
    )
par(mar = c(3, 4.2, 4, 2))
FACSViolinBox(flow_df, "CD2",   "Blues")
FACSViolinBox(flow_df, "CD4",   "Purples")
FACSViolinBox(flow_df, "CD43",  "Greens")
FACSViolinBox(flow_df, "CD200", "Greys")
dev.off()





