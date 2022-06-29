# 2021-12-27


# Load packages and source code -------------------------------------------

library("RColorBrewer")



# Define labels -----------------------------------------------------------

controls_labels <- list(
  "NT"   = c("Non-targeting", "controls"),
  "Pos"  = c("Positive", "controls", expression("(" * italic("GBA") * " gene)")),
  "Gene" = c("Genes in ", "CRISPRa", "library")
)

pos_ctrl_color <- brewer.pal(5, "Reds")[[4]]
NT_ctrl_color  <- brewer.pal(5, "Blues")[[4]]
custom_color   <- brewer.pal(5, "Purples")[[4]]
controls_colors <- c(NT_ctrl_color, pos_ctrl_color, "black")



# Define pairs of metrics -------------------------------------------------

pairs_list <- list(

   "deltaNT" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_deltaNT",
      "title" = "Volcano plot (p values from untransformed data)"
   ),
   "activation" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_act",
      "title" = "Volcano plot (p values from % activation)"
   ),
   "log2" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_log2",
      "title" = "Volcano plot (p values from log2-transformed data)"
   ),
   "act_log2" = c(
      "x_var" = "Log2FC_rep1",
      "y_var" = "p_value_act_log2",
      "title" = "Volcano plot (p values from % activation, log2 data)"
   ),


   "deltaNT_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_deltaNT_Glo",
      "title" = "Volcano plot (untransformed data, CellTiter-Glo-norm.)"
   ),
   "activation_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_act_Glo",
      "title" = "Volcano plot (% activation, CellTiter-Glo-norm.)"
   ),
   "log2_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_log2_Glo",
      "title" = "Volcano plot (log2, CellTiter-Glo-normalized)"
   ),
   "act_log2_Glo" = c(
      "x_var" = "Log2FC_Glo_rep1",
      "y_var" = "p_value_act_log2_Glo",
      "title" = "Volcano plot (% activation, log2, CellTiter-Glo-norm.)"
   )
)



# Define functions --------------------------------------------------------

VolcanoFlashPlot <- function(input_df,
                             log_fc_column,
                             y_column,
                             show_only_genes   = FALSE,
                             show_title        = "",
                             point_size        = 0.8,
                             indicate_p_values = 0.05,              # for the grey background
                             indicate_log2FCs  = 0.5,               # for the grey background
                             label_p_values    = indicate_p_values, # for the gene names
                             label_log2FCs     = indicate_log2FCs,  # for the gene names
                             indicate_areas    = FALSE,
                             indicate_lines    = FALSE,
                             label_points      = FALSE,
                             tiny_labels       = FALSE,
                             use_mai           = NULL,
                             use_mgp           = c(2.8, 0.55, 0),
                             x_axis_mgp        = use_mgp,
                             ...
                             ) {


  ## Prepare data for plotting
  if (grepl("_rep", log_fc_column, fixed = TRUE)) {
    rep2_column <- sub("_rep1", "_rep2", log_fc_column, fixed = TRUE)
    log_fc_vec <- rowMeans(input_df[, c(log_fc_column, rep2_column)])
  } else {
    log_fc_vec <- input_df[, log_fc_column]
  }
  y_value_vec <- input_df[, y_column]

  if (grepl("PercActivation", log_fc_column, fixed = TRUE)) {
    x_label <- "% activation"
    log_fc_vec <- log_fc_vec * 100
  } else {
    x_label <- expression("log"["2"] ~ "fold change")
  }
  if (grepl("SSMD", y_column, fixed = TRUE)) {
    y_label <- "SSMD"
  } else {
    y_value_vec <- -log10(y_value_vec)
    y_label <- expression("" - "log"["10"] ~ plain("p") ~ "value")
  }

  ## Prepare data subsets
  are_NT      <- input_df[, "Is_NT_ctrl"]
  are_posctrl <- input_df[, "Is_pos_ctrl"]
  are_gene    <- !(is.na(input_df[, "Entrez_ID"]))
  if ("Custom_color" %in% names(input_df)) {
    are_custom_color <- input_df[, "Custom_color"]
  } else {
    are_custom_color <- rep(FALSE, nrow(input_df))
  }
  if (show_only_genes) {
    are_valid <- are_gene
  } else {
    are_valid <- are_NT | are_posctrl | are_gene | are_custom_color
  }

  ## Prepare graphical parameters
  if (is.null(use_mai)) {
    use_mai <- c(0.85, 0.8, 0.7, 1.5)
    if (show_only_genes) {
      use_mai[[4]] <- 0.7
    }
  }

  old_mar <- par(mai = use_mai)

  ## Set up the plot region
  plot(1,
       xlim = range(log_fc_vec[are_valid]),
       ylim = range(c(0, y_value_vec[are_valid])),
       type = "n",
       bty  = "n",
       axes = FALSE,
       ann  = FALSE
       )
  use_tcl <- -0.35
  axis(1, mgp = x_axis_mgp, lwd = par("lwd"), tcl = use_tcl)
  mtext(x_label, side = 1, line = x_axis_mgp[[1]], cex = par("cex"))
  axis(2, mgp = use_mgp, lwd = par("lwd"), las = 1, tcl = use_tcl)
  mtext(y_label, side = 2, line = use_mgp[[1]], cex = par("cex"))

  ## Draw indicator lines
  abline(h = 0, lty = "dotted", col = "grey50")
  abline(v = 0, lty = "dotted", col = "grey50")
  if (indicate_lines) {
    if (!(is.null(indicate_p_values))) {
      abline(h = -log10(indicate_p_values), col = "grey90")
    }
    if (!(is.null(indicate_log2FCs))) {
      abline(v = c(-(indicate_log2FCs), indicate_log2FCs),
             col = "grey90"
             )
    }
  }

  ## Indicate the plot regions that lie above the cutoffs
  if (indicate_areas) {
    stopifnot(length(indicate_p_values) == 1)
    stopifnot(length(indicate_log2FCs) == 1)
    rect(xleft   = par("usr")[[1]],
         xright  = -(indicate_log2FCs),
         ybottom = -log10(indicate_p_values),
         ytop    = par("usr")[[4]],
         col     = "gray92",
         border  = NA
         )
    rect(xleft   = indicate_log2FCs,
         xright  =  par("usr")[[2]],
         ybottom = -log10(indicate_p_values),
         ytop    = par("usr")[[4]],
         col     = "gray92",
         border  = NA
         )
  }
  box()


  ## Draw the points of the volcano/flashlight plot
  points(log_fc_vec[are_gene],
         y_value_vec[are_gene],
         pch = 16,
         col = adjustcolor("black", alpha.f = 0.3),
         cex = point_size
         )

  if (!(show_only_genes)) {

    if (any(are_custom_color)) {
      points(log_fc_vec[are_custom_color],
             y_value_vec[are_custom_color],
             pch = 16,
             col = adjustcolor(custom_color, alpha.f = 0.5),
             cex = point_size
             )
    }

    points(log_fc_vec[are_posctrl],
           y_value_vec[are_posctrl],
           pch = 16,
           col = adjustcolor(pos_ctrl_color, alpha.f = 0.5),
           cex = point_size
           )

    points(log_fc_vec[are_NT],
           y_value_vec[are_NT],
           pch = 16,
           col = adjustcolor(NT_ctrl_color, alpha.f = 0.5),
           cex = point_size
           )

    ## Draw a legend for the points
    DrawSideLegend(labels_list = controls_labels, use_colors = controls_colors, ...)
  } else {
    are_custom_gene <- are_gene & are_custom_color
    if (any(are_custom_color)) {
      points(log_fc_vec[are_custom_gene],
             y_value_vec[are_custom_gene],
             pch = 16,
             col = adjustcolor(custom_color, alpha.f = 0.5),
             cex = point_size
             )
    }
  }

  ## Highlight genes that pass the cutoffs
  are_to_highlight <- are_gene &
                     (y_value_vec >= (-(log10(indicate_p_values)))) &
                     (abs(log_fc_vec) >= indicate_log2FCs)
  if (label_points) {
    are_to_label <- are_to_highlight & (abs(log_fc_vec) >= label_log2FCs)
    points(log_fc_vec[are_to_label],
           y_value_vec[are_to_label],
           pch = 16,
           cex = point_size
           )
    text(x      = log_fc_vec[are_to_label],
         y      = y_value_vec[are_to_label],
         labels = input_df[are_to_label, "Gene_symbol"],
         pos    = 3,
         offset = if (tiny_labels) 0.15 else 0.2,
         cex    = if (tiny_labels) 0.12 else 0.6,
         col    = if (tiny_labels) "orange" else "black",
         font   = 4,
         xpd    = NA
         )

  } else if (indicate_areas) {
    points(log_fc_vec[are_to_highlight],
           y_value_vec[are_to_highlight],
           pch = 16,
           col = "black",
           cex = point_size
           )
  }

  title(show_title, cex.main = 1.1)
  par(old_mar)
  return(invisible(NULL))
}



ExportAllVolcanoAndFlashlightPlots <- function(input_df) {

  plot_types <- c("Volcano", "Dual flashlight (logFC)", "Dual flashlight (% activation)")

  for (use_device in c("none", "pdf", "png")) {

    if (use_device == "pdf") {
      message("Exporting PDF files...")
    } else if (use_device == "png") {
      message("Exporting PNG files...")
    } else if (use_device == "none") {
      message("Diplaying plots in the editor...")
    }

    for (plot_type in plot_types) {

      if (plot_type == "Volcano") {
        message("... creating volcano plots...")
        folder_name <- "Volcano plots"
        PNG_prefix <- "Volcano plot"
        x_sub <- "Log2FC"
        y_sub <- "p_value"
        title_prefix <- "Volcano plot"
      } else if (plot_type == "Dual flashlight (logFC)") {
        message("... creating dual-flashlight plots (using log2FC)...")
        folder_name <- "Dual flashlight plots (logFC)"
        PNG_prefix <- "Dual flashlight (logFC)"
        x_sub <- "Log2FC"
        y_sub <- "SSMD"
        title_prefix <- "Dual flashlight plot"
      } else if (plot_type == "Dual flashlight (% activation)") {
        message("... creating dual-flashlight plots (using % activation)...")
        folder_name <- "Dual flashlight plots (pct activation)"
        PNG_prefix <- "Dual flashlight (% activation)"
        x_sub <- "PercActivation"
        y_sub <- "SSMD"
        title_prefix <- "Dual flashlight plot"
      }

      for (only_genes in c(FALSE, TRUE)) {

        file_name <- folder_name

        if (only_genes) {
          PDF_name <- paste0(file_name, " - genes only.pdf")
          use_width <- base_width
        } else {
          PDF_name <- paste0(file_name, " - with controls.pdf")
          use_width <- base_width + 0.8
        }
        if (use_device == "pdf") {
          pdf(file = file.path(output_dir, "Figures", folder_name, PDF_name),
              width = use_width, height = base_height
              )
        }
        for (i in seq_along(pairs_list)) {
          plot_title <- sub("Volcano plot", title_prefix, pairs_list[[i]][["title"]])
          if (plot_type != "Volcano") {
            plot_title <- sub("p values", "SSMD", plot_title, fixed = TRUE)
          }
          if (use_device == "png") {
            PNG_name <- paste0(i, ") ", gsub("%", "percent", plot_title, fixed = TRUE))
            if (only_genes) {
              PNG_name <- paste0(PNG_name, " - genes only.png")
            } else {
              PNG_name <- paste0(PNG_name, " - with controls.png")
            }
            png(filename = file.path(output_dir, "Figures", folder_name, "PNGs", PNG_name),
                width = use_width, height = base_height, units = "in", res = 600
                )
          }
          VolcanoFlashPlot(input_df,
                           sub("Log2FC", x_sub, pairs_list[[i]][["x_var"]], fixed = TRUE),
                           sub("p_value", y_sub, pairs_list[[i]][["y_var"]], fixed = TRUE),
                           show_title = FormatPlotMath(plot_title),
                           show_only_genes = only_genes
                           )
          if (use_device == "png") {
            dev.off()
          }
        }

        if (use_device == "pdf") {
          dev.off()
        }
      }
    }
  }
  return(invisible(NULL))
}


