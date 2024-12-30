## 2024-11-18



# Define labels -----------------------------------------------------------

feature_labels <- c(
  "promoter1_hU6"  = "promoter 1", # (hU6)
  "sg1"            = "sgRNA 1",
  "tracrRNA1"      = "tracrRNA 1",
  "polyT_1"        = "poly(T)",

  "EM7_promoter"   = "EM7 promoter",
  "pre_TpR"        = "unannotated",
  "TpR_DHFR"       = "DHFR (TpR)",
  "polyT_TpR"      = "poly(T)",

  "promoter2_mU6"  = "promoter 2", # (mU6)
  "sg2"            = "sgRNA 2",
  "tracrRNA2"      = "tracrRNA 2",
  "polyT_2"        = "poly(T)",

  "promoter3_hH1"  = "promoter 3", # (hH1)
  "sg3"            = "sgRNA 3",
  "tracrRNA3"      = "tracrRNA 3",
  "polyT_3"        = "poly(T)",

  "promoter4_h7SK" = "promoter 4", # (h7SK)
  "sg4"            = "sgRNA 4",
  "tracrRNA4"      = "tracrRNA 4",
  "polyT_4"        = "poly(T)"
)



# Define functions --------------------------------------------------------

PValueToStars <- function(p_value_vec) {
  ifelse(p_value_vec < 0.001,
         "***",
         ifelse(p_value_vec < 0.01,
                "**",
                ifelse(p_value_vec < 0.05, "*", NA_real_)
                ))
}



ErrorDumbBells <- function(display_df,
                           sg_A = NULL,
                           sg_B = NULL,
                           x_upper_limit = 25,
                           x_axis_label = "Error rate"
                           ) {

  ## Check parameters
  stopifnot(max(display_df[, c("Fraction_incorrect_nonswitched", "Fraction_incorrect_switched")]) < (x_upper_limit / 100))
  sg_NULL <- is.null(sg_A) + is.null(sg_B)
  if (sg_NULL == 1L) {
    stop("The parameters 'sg_A' and 'sg_B' must either both be specified, or both be NULL!")
  }

  ## Check which sgRNAs are involved (i.e., the sgRNAs between which a
  ## template switch is either seen or not seen to have occurred.)
  are_error_free <- rowSums(display_df[, c("Num_incorrect_nonswitched", "Num_incorrect_switched")] == 0) == 2
  error_free_sgs <- display_df[are_error_free, "Feature"]
  stopifnot(all(error_free_sgs %in% paste0("sg", 1:4)))
  if (sum(are_error_free) == 4) {
    num_features_per_sg <- 3L
    if (sg_NULL != 0) {
      stop("Please specify between which pair of sgRNAs a template switch was looked for!")
    }
  } else if (sum(are_error_free) == 2) {
    if (sg_NULL != 2) {
      stop("The pair of sgRNAs involved in the potential template switch can",
           " be inferred from the data and should not be specified!"
           )
    }
    num_features_per_sg <- 4L
    sg_A <- as.integer(substr(error_free_sgs[[1]], 3, 3))
    sg_B <- as.integer(substr(error_free_sgs[[2]], 3, 3))
  } else {
    stop("Unexpected number of sgRNAs with an error rate of 0!")
  }

  ## Identify the features that should be highlighted (e.g. labelled in bold)
  ## because they lie in between the pair of sgRNAs between which a template switch may have occurred.
  are_marked <- display_df[, "Feature"] %in% paste0("sg", c(sg_A, sg_B))
  highlighted_indices <- which(are_marked)[[1]]:which(are_marked)[[2]]
  are_highlighted <- seq_along(are_marked) %in% highlighted_indices


  ## Count the number of reads with or without a template switch.
  num_switched <- unique(round(display_df[!(are_error_free), "Num_incorrect_switched"] /
                               display_df[!(are_error_free), "Fraction_incorrect_switched"]
                               ))
  num_non_switched <- unique(round(display_df[!(are_error_free), "Num_incorrect_nonswitched"] /
                                   display_df[!(are_error_free), "Fraction_incorrect_nonswitched"]
                                   ))

  ## If all sgRNAs are fully mapped (i.e., error-free), all 4 sgRNA sequences can be removed from the plot.
  if (sum(are_error_free) == 4) {
    display_df <- display_df[!(are_error_free), ]
    are_marked <- rep(FALSE, nrow(display_df))
    are_highlighted <- are_highlighted[!(are_error_free)]
  }

  ## Prepare data for plotting
  x_vec_1 <- display_df[, "Fraction_incorrect_nonswitched"] * 100
  x_vec_2 <- display_df[, "Fraction_incorrect_switched"] * 100


  ## Visually group related features by positioning them closer together on the y axis
  num_features <- nrow(display_df)
  feature_groups <- rep(NA, num_features)
  for (i in 1:4) {
    are_this_sg <- grepl(paste0("^(promoter|sg|tracrRNA|polyT_)", i), display_df[, "Feature"])
    feature_groups[are_this_sg] <- i
    stopifnot(sum(are_this_sg) == num_features_per_sg)
  }
  feature_groups[is.na(feature_groups)] <- 5L
  y_pos <- num_features + 1 - RepositionByGroups(feature_groups, gap_ratio = 1.75)


  ## Define colors
  switched_dark_color     <- adjustcolor(hcl.colors(9, "Reds")[[3]], alpha.f = 0.6)
  switched_light_color    <- adjustcolor(hcl.colors(9, "Reds")[[5]], alpha.f = 0.5)
  nonswitched_dark_color  <- adjustcolor(hcl.colors(9, "Blues")[[3]], alpha.f = 0.6)
  nonswitched_light_color <- adjustcolor(hcl.colors(9, "Blues")[[5]], alpha.f = 0.5)

  ## Set up plot region
  par(mar = c(5, 8, 3.8, 3))
  plot.new()
  plot.window(xlim = c(0, x_upper_limit),
              ylim = c(1 - (1 / ((1 + sqrt(5)) / 2)), num_features),
              xaxs = "i",
              yaxs = "i"
              )

  ## Add a title/legend to the top of the plot
  text(x      = mean(par("usr")[1:2]),
       y      = par("usr")[[3]] - diff(grconvertY(c(0, 4), from = "lines", to = "user")),
       labels = paste0(num_switched, " reads with a template switch and ",
                       num_non_switched, " reads without a template switch ",
                       "were included."
                       ),
       cex    = 0.7,
       col    = "gray70",
       xpd    = NA,
       adj    = c(0.5, 0.5)
       )

  switch_text <- paste0("Template switch between sg", sg_A, " & sg", sg_B, "?")
  switch_y <- max(y_pos) + diff(grconvertY(c(0, 2), from = "lines", to = "user"))
  text(x      = 0,
       y      = switch_y,
       labels = switch_text,
       xpd    = NA,
       adj    = c(0, 0.5)
       )
  point_x <- 0 + strwidth(switch_text, units = "user") +
             diff(grconvertX(c(0, 0.9), from = "chars", to = "user"))
  points(x   = point_x,
         y   = switch_y,
         pch = 21,
         col = switched_dark_color,
         bg  = switched_light_color,
         xpd = NA
         )
  text(x      = point_x + diff(grconvertX(c(0, 0.4), from = "chars", to = "user")),
       y      = switch_y,
       labels = "Yes",
       xpd    = NA,
       adj    = c(0, 0.5)
       )

  points(x   = point_x + diff(grconvertX(c(0, 2.5), from = "chars", to = "user")),
         y   = switch_y,
         pch = 21,
         col = nonswitched_dark_color,
         bg  = nonswitched_light_color,
         xpd = NA
         )
  text(x      = point_x + diff(grconvertX(c(0, 2.9), from = "chars", to = "user")),
       y      = switch_y,
       labels = "No",
       xpd    = NA,
       adj    = c(0, 0.5)
       )


  ## Add an x axis and gridlines
  tick_pos <- axTicks(1)
  segments(x0  = tick_pos,
           y0  = par("usr")[[3]],
           y1  = max(y_pos),
           col = "gray90",
           xpd = NA
           )
  axis(1, mgp = c(3, 0.5, 0), tcl = -0.4, at = tick_pos, labels = paste0(format(tick_pos), "%"))
  title(xlab = x_axis_label, mgp = c(1.8, 1, 0))


  ## Label features on the y axis
  text(x      = par("usr")[[1]] - diff(grconvertX(c(0, 0.5), from = "lines", to = "user")),
       y      = y_pos,
       labels = feature_labels[display_df[, "Feature"]],
       font   = ifelse(are_highlighted, 2, 1),
       col    = ifelse(are_highlighted, "black", "gray30"),
       adj    = c(1, 0.5),
       xpd    = NA
       )


  ## Draw stars to indicate p values
  text(x      = par("usr")[[2]] + diff(grconvertX(c(0, 0.7), from = "lines", to = "user")),
       y      = y_pos[!(are_marked)] - diff(grconvertY(c(0, 0.1), from = "lines", to = "user")),
       labels = PValueToStars(p.adjust(display_df[!(are_marked), "P_two_sided"])),
       adj    = c(0.5, 0.5),
       col    = "gray50",
       xpd    = NA
       )


  ## Add gridlines for the y axis
  segments(x0  = 0,
           x1  = x_upper_limit,
           y0  = y_pos,
           col = "gray60",
           lty = "dotted",
           xpd = NA
           )


  ## Add line segments between the points of the dumbbells
  point_radius <- (par("cxy")[[1]] / pi) * par("cex")
  lower_vec  <- mapply(min, x_vec_1[!(are_marked)], x_vec_2[!(are_marked)]) + point_radius
  higher_vec <- mapply(max, x_vec_1[!(are_marked)], x_vec_2[!(are_marked)]) - point_radius
  are_to_show <- lower_vec < higher_vec
  segments(x0  = lower_vec[are_to_show],
           x1  = higher_vec[are_to_show],
           y0  = y_pos[!(are_marked)][are_to_show],
           col = "gray75",
           lwd = 2,
           xpd = NA
           )

  ## Plot points
  point_cex <- 1.2
  point_lwd <- 1.5
  points(x   = x_vec_1[!(are_marked)],
         y   = y_pos[!(are_marked)],
         col = nonswitched_dark_color,
         bg  = nonswitched_light_color,
         cex = point_cex,
         lwd = point_lwd,
         pch = 21,
         xpd = NA
         )
  points(x   = x_vec_2[!(are_marked)],
         y   = y_pos[!(are_marked)],
         col = switched_dark_color,
         bg  = switched_light_color,
         cex = point_cex,
         lwd = point_lwd,
         pch = 21,
         xpd = NA
         )

  return(invisible(NULL))
}

