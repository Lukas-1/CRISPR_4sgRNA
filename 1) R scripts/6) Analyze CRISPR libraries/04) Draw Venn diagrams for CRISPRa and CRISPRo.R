### 26th February 2021 ##





# Import packages and source code -----------------------------------------

library("venneuler")
library("sp")
library("rgeos")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))




# Define folder paths -----------------------------------------------------

CRISPR_root_directory    <- "~/CRISPR"
RData_directory          <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(RData_directory, "1) General")
CRISPRa_RData_directory  <- file.path(RData_directory, "2) CRISPRa")
CRISPRko_RData_directory <- file.path(RData_directory, "3) CRISPRko")
file_output_directory    <- file.path(CRISPR_root_directory, "5) Output", "Analysis", "Venn diagrams")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "12) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData"))
load(file.path(CRISPRa_RData_directory, "19) For problematic genes, pick 4 guides without reference to the TSS.RData"))
load(file.path(CRISPRko_RData_directory, "11) Pick 4 guides per gene.RData"))




# Define functions --------------------------------------------------------

MixColors <- function (rcol_in) {
  # This was taken entirely from eulerr:::mix_colors
  rgb_in <- t(grDevices::col2rgb(rcol_in))
  lab_in <- grDevices::convertColor(rgb_in, "sRGB", "Lab",
                                    scale.in = 255
                                    )
  mean_col <- colMeans(lab_in)
  rgb_out <- grDevices::convertColor(mean_col, "Lab",
                                     "sRGB", scale.out = 1
                                     )
  grDevices::rgb(rgb_out)
}



Rotate <- function(x, y = NULL, mx = NULL, my = NULL, theta = pi/3, asp = 1) {
  ## Adapted from DescTools::Rotate()
  xy <- xy.coords(x, y)
  if (is.null(mx))
    mx <- mean(xy$x)
  if (is.null(my))
    my <- mean(xy$y)
  dx <- xy$x - mx
  dy <- xy$y - my
  ptx <- mx + cos(theta) * dx - sin(theta) * dy/asp
  pty <- my + sin(theta) * dx * asp + cos(theta) * dy
  return(cbind("x" = ptx, "y" = pty))
}


ThreeSetsVenn <- function(use_factor,
                          x_limits        = c(0.2, 0.8),
                          y_limits        = c(0.17, 0.83),
                          use_colors      = c("#C2DCF5", "#C0EABD", "#FFBCAD"),
                          flip_horizontal = FALSE,
                          flip_vertical   = FALSE,
                          rotate_by       = 0,
                          scale_by        = 1,
                          positions_list  = NULL,
                          line_positions  = NULL,
                          shift_x         = 0,
                          shift_y         = 0,
                          embed_PNG       = FALSE
                          ) {

  x_range <- x_limits[[2]] - x_limits[[1]]
  y_range <- y_limits[[2]] - y_limits[[1]]

  x_mid <- x_limits[[1]] + (x_range / 2)
  y_mid <- y_limits[[1]] + (y_range / 2)

  new_half_x_range <- (x_range / scale_by) / 2
  new_half_y_range <- (y_range / scale_by) / 2
  x_limits <- c(x_mid - new_half_x_range, x_mid + new_half_x_range)
  y_limits <- c(x_mid - new_half_y_range, x_mid + new_half_y_range)

  numbers_vec <- tabulate(use_factor)
  names(numbers_vec) <- gsub(", ", "&", levels(use_factor), fixed = TRUE)
  venn_object <- venneuler(numbers_vec)

  if (flip_horizontal) {
    venn_object[["centers"]][, "x"] <- -(venn_object[["centers"]][, "x"] - x_mid) + x_mid
  }
  if (flip_vertical) {
    venn_object[["centers"]][, "y"] <- -(venn_object[["centers"]][, "y"] - y_mid) + y_mid
  }
  venn_object[["centers"]] <- Rotate(venn_object[["centers"]], theta = (rotate_by * pi) / 180)


  line_positions <- lapply(line_positions, function(x) {
    x[1:2] <- x[1:2] + shift_x
    return(x)
  })
  positions_list <- lapply(positions_list, function(x) {
    x[[1]] <- x[[1]] + shift_x
    return(x)
  })
  venn_object[["centers"]][, "x"] <- venn_object[["centers"]][, "x"] + shift_x

  line_positions <- lapply(line_positions, function(x) {
    x[3:4] <- x[3:4] + shift_y
    return(x)
  })
  positions_list <- lapply(positions_list, function(x) {
    x[[2]] <- x[[2]] + shift_y
    return(x)
  })
  venn_object[["centers"]][, "y"] <- venn_object[["centers"]][, "y"] + shift_y


  num_circles <- length(venn_object[["diameters"]])

  edges <- 200
  s <- seq.int(edges) / edges * 2 * pi
  sx <- cos(s) / 2 # VD uses diameter, not radius
  sy <- sin(s) / 2


  circles_mat_list <- lapply(seq_len(num_circles), function(i) {
    cbind("x" = venn_object[["centers"]][i, 1] + venn_object[["diameters"]][[i]] * sx,
          "y" = venn_object[["centers"]][i, 2] + venn_object[["diameters"]][[i]] * sy
          )
  })

  if (embed_PNG) {
    PDF_device <- dev.cur()
    temp_path <- file.path(file_output_directory, "temp.png")
    par(mar = rep(0, 4), omi = rep(0.1, 4))
    temp_width  <- par("pin")[[1]]
    temp_height <- par("pin")[[2]]
    current_par <- par(no.readonly = TRUE)
    png(filename = temp_path,
        width    = temp_width,
        height   = temp_height,
        units    = "in",
        res      = 900,
        bg       = "white"
        )
    par(lwd = current_par[["lwd"]])
    par(cex = current_par[["cex"]])
    plot.new()
    par(mar = rep(0, 4))
  } else {
    plot.new()
    par(mar = rep(0, 4), omi = rep(0.1, 4))
  }
  plot.window(x_limits, y_limits, "", asp = 1)


  # Plot all circles
  for (i in seq_along(circles_mat_list)) {
    polygon(circles_mat_list[[i]][, "x"],
            circles_mat_list[[i]][, "y"],
            col = use_colors[[i]],
            border = NA
            )
  }

  MakePolygon <- function(x, y, use_ID = 1) {
    x <- c(x, x[1])
    y <- c(y, y[1])
    result_p <- SpatialPolygons(list(Polygons(list(Polygon(cbind(rev(x), rev(y)))), ID = use_ID)))
    return(result_p)
  }


  comb_mat <- combn(1:3, 2)
  for (i in seq_len(ncol(comb_mat))) {
    index_1 <- comb_mat[1, i]
    index_2 <- comb_mat[2, i]
    x1 <- circles_mat_list[[index_1]][, "x"]
    x2 <- circles_mat_list[[index_2]][, "x"]
    y1 <- circles_mat_list[[index_1]][, "y"]
    y2 <- circles_mat_list[[index_2]][, "y"]
    p1 <- MakePolygon(x1, y1, use_ID = 1)
    p2 <- MakePolygon(x2, y2, use_ID = 3)
    tmp <- gIntersection(p1, p2)
    mixed_color <- MixColors(c(use_colors[[index_1]], use_colors[[index_2]]))
    plot(tmp, add = TRUE, col = mixed_color, border = NA)
  }
  x3 <- circles_mat_list[[1]][, "x"]
  y3 <- circles_mat_list[[1]][, "y"]
  p3 <- MakePolygon(x3, y3, use_ID = 3)
  new_tmp <- gIntersection(tmp, p3)
  mixed_color <- MixColors(use_colors)
  plot(new_tmp, add = TRUE, col = mixed_color, border = NA)


  if (embed_PNG) {
    dev.off()
    raster_array <- png::readPNG(temp_path)
    file.remove(temp_path)
    dev.set(PDF_device)
    plot.new()
    par(mar = rep(0, 4), omi = rep(0.1, 4))
    plot.window(x_limits, y_limits, "", asp = 1)
    rasterImage(raster_array,
                xleft = par("usr")[[1]], xright = par("usr")[[2]],
                ybottom = par("usr")[[3]], ytop = par("usr")[[4]]
                )
  }


  if (!(is.null(positions_list))) {
    number_table <- table(use_factor)
    for (combo in names(positions_list)) {
      x_pos <- positions_list[[combo]][["x"]]
      y_pos <- positions_list[[combo]][["y"]]
      is_combo <- combo %in% names(number_table)
      if (is_combo && !(grepl(",", combo, fixed = TRUE))) {
        use_label <- combo
        if (use_label %in% c("GPPa", "GPPo")) {
          use_label <- "CRISPick"
        }
        y_distance <- (par("usr")[[4]] - par("usr")[[3]]) * 0.048
        text(x = x_pos, y = y_pos + y_distance, labels = use_label, font = 2)
        y_pos <- y_pos - y_distance
      }
      if (combo %in% names(number_table)) {
        draw_text <- number_table[[combo]]
        if (draw_text >= 10000) {
          draw_text <- FormatThousands(draw_text)
        }
      } else {
        draw_text <- combo
      }
      text(x = x_pos, y = y_pos, labels = draw_text)
    }

    if (!(is.null(line_positions))) {
      for (line_position in line_positions) {
        do.call(segments, c(as.list(line_position), lwd = 1, col = "gray80"))
      }
    }
  }
  return(invisible(NULL))
}



Get4sgFactor <- function(CRISPR_df) {
  CRISPR_df <- FilterCRISPRDf(CRISPR_df)
  are_chosen <- Are4sg(CRISPR_df, sublibraries_all_entrezs_list, show_messages = FALSE)
  all_sources_fac <- ReformatSourceToFactor(CRISPR_df[["Source"]][are_chosen])
  return(all_sources_fac)
}



# Define label positions --------------------------------------------------

CRISPRko_number_labels_list <- list(

  "T.spiezzo" = c(
    "x" = 0.77, "y" = 0.72
  ),

  "Brunello" = c(
    "x" = 0.29, "y" = 0.337
  ),

  "TKOv3" = c(
    "x" = 0.7614468, "y" = 0.4131915
  ),

  "GPPo" = c(
    "x" = 0.42, "y" = 0.73
  ),

  "Brunello, TKOv3" = c(
    "x" = 0.56, "y" = 0.26
  ),

  "Brunello, GPPo" = c(
    "x" = 0.355, "y" = 0.53
  ),

  "TKOv3, GPPo" = c(
    "x" = 0.632, "y" = 0.5949787
  ),

  "Brunello, TKOv3, GPPo" = c(
    "x" = 0.5551489, "y" = 0.475
  )
)


CRISPRko_line_positions <- list(
  c("x0" = CRISPRko_number_labels_list[["Brunello, TKOv3"]][["x"]],
    "x1" = 0.55 ,
    "y0" = CRISPRko_number_labels_list[["Brunello, TKOv3"]][["y"]] + 0.027,
    "y1" = 0.325
    )
)



CRISPRa_number_labels_list <- list(

  "T.gonfio"= c(
    "x" = 0.69, "y" = 0.80
  ),

  "Calabrese" = c(
    "x" = 0.355, "y" = 0.36
  ),

  "hCRISPRa-v2" = c(
    "x" = 0.77, "y" = 0.43
  ),

  "GPPa" = c(
    "x" = 0.385, "y" = 0.6985547
  ),

  "Calabrese, hCRISPRa-v2" = c(
    "x" = 0.47, "y" = 0.27
  ),

  "Calabrese, GPPa" = c(
    "x" = 0.3259246, "y" = 0.4718941
  ),

  "hCRISPRa-v2, GPPa" = c(
    "x" = 0.5797845, "y" = 0.6187702
  ),

  "Calabrese, hCRISPRa-v2, GPPa" = c(
    "x" = 0.4746140, "y" = 0.4682675
  )
)


CRISPRa_line_positions <- list(
  c("x0" = CRISPRa_number_labels_list[["Calabrese, hCRISPRa-v2"]][["x"]],
    "x1" = 0.475 ,
    "y0" = CRISPRa_number_labels_list[["Calabrese, hCRISPRa-v2"]][["y"]] + 0.027,
    "y1" = 0.325
    )
)




# Prepare Venn diagram ----------------------------------------------------

CRISPRa_fac <- Get4sgFactor(merged_replaced_CRISPRa_df)
levels(CRISPRa_fac) <- sub("GPP", "GPPa", levels(CRISPRa_fac), fixed = TRUE)

CRISPRko_fac <- Get4sgFactor(merged_CRISPRko_df)
levels(CRISPRko_fac) <- sub("GPP", "GPPo", levels(CRISPRko_fac), fixed = TRUE)

combined_vec <- c(as.character(CRISPRa_fac), as.character(CRISPRko_fac))
combined_fac <- factor(combined_vec,
                       levels = c(levels(CRISPRa_fac), levels(CRISPRko_fac))
                       )




# Create Euler plots using the eulerr package -----------------------------

three_colors <- c(brewer.pal(9, "Blues")[[4]], brewer.pal(9, "Greens")[[2]], brewer.pal(9, "Greys")[[3]])
venn_colors <- rep(three_colors, times = 2)

three_colors <- c(brewer.pal(9, "Blues")[[5]], brewer.pal(9, "Greens")[[5]], brewer.pal(9, "Reds")[[5]])
three_colors <- vapply(three_colors, Palify, fraction_pale = 0.6, "")

three_colors <- c("#C4DFEF", "#C7E7C8", "#FDC3B7")

venn_colors <- rep(three_colors, times = 2)


combined_venn <- PlotVennDiagram(combined_fac, draw_plot = TRUE,
                                 quantities_font = 1, quantities_cex = 1,
                                 use_padding = 0.5, use_colors = venn_colors
                                 )




# Create Euler plots using the venneuler package --------------------------

ThreeSetsVenn(CRISPRa_fac, flip_vertical = TRUE, rotate_by = 188,
              scale_by = length(CRISPRa_fac) / length(CRISPRko_fac)
              )
ThreeSetsVenn(CRISPRko_fac, rotate_by = 170)





# Export Venn diagrams ----------------------------------------------------

pdf(file.path(file_output_directory, paste0("Venn diagram - CRISPRa & CRISPRo", ".pdf")),
    height = 3, width = 6
    )

print(combined_venn)
dev.off()



pdf(file.path(file_output_directory, paste0("CRISPRa Venn diagram", ".pdf")),
    width = horizontal_width, height = horizontal_height
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ThreeSetsVenn(CRISPRa_fac, flip_vertical = TRUE, rotate_by = 188,
              scale_by = length(CRISPRa_fac) / length(CRISPRko_fac),
              positions_list = CRISPRa_number_labels_list,
              line_positions = CRISPRa_line_positions,
              shift_x = -0.08, shift_y = -0.015
              )
dev.off()


pdf(file.path(file_output_directory, paste0("CRISPRko Venn diagram", ".pdf")),
    width = horizontal_width, height = horizontal_height
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ThreeSetsVenn(CRISPRko_fac, rotate_by = 170,
              positions_list = CRISPRko_number_labels_list,
              line_positions = CRISPRko_line_positions,
              shift_x = -0.05
              )
dev.off()





library("devEMF")

emf(file.path(file_output_directory, paste0("CRISPRa Venn diagram", ".emf")),
    width = horizontal_width, height = horizontal_height, emfPlus = FALSE
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ThreeSetsVenn(CRISPRa_fac, flip_vertical = TRUE, rotate_by = 188,
              scale_by = length(CRISPRa_fac) / length(CRISPRko_fac),
              positions_list = CRISPRa_number_labels_list,
              line_positions = CRISPRa_line_positions,
              shift_x = -0.08, shift_y = -0.015,
              embed_PNG = TRUE
              )
dev.off()


emf(file.path(file_output_directory, paste0("CRISPRko Venn diagram", ".emf")),
    width = horizontal_width, height = horizontal_height, emfPlus = FALSE
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ThreeSetsVenn(CRISPRko_fac, rotate_by = 170,
              positions_list = CRISPRko_number_labels_list,
              line_positions = CRISPRko_line_positions,
              shift_x = -0.05,
              embed_PNG = TRUE
              )
dev.off()
