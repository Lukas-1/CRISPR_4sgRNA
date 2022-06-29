# 2021-11-13


# Import packages and source code -----------------------------------------

library("readxl")

project_dir <- "~/R_projects/CRISPRa_TF"
misc_functions_dir <- file.path(project_dir, "1_R_scripts", "1_R_functions", "4_Misc_functions")
source(file.path(misc_functions_dir, "01_Converting_plate_layouts.R"))



# Define folder paths -----------------------------------------------------

input_dir         <- file.path(project_dir, "2_input")
general_rdata_dir <- file.path(project_dir, "3_R_objects", "1_General")
output_dir        <- file.path(project_dir,"4_output")

plate_layout_96_path <- file.path(input_dir, "2021-09-29_96formatplates.xlsx")
long_layout_96_path  <- file.path(input_dir, "ACTIVATION LIBRARY- Plate Layout.xlsx")



# Load data ---------------------------------------------------------------

load(file.path(input_dir, "CRISPRa_4sg_df.RData"))



# Read in data ------------------------------------------------------------

sheets_vec <- readxl::excel_sheets(plate_layout_96_path)
layout_96_mat_list <- lapply(sheets_vec,
                             function(x) as.matrix(suppressMessages(read_excel(plate_layout_96_path, sheet = x)),
                                                   stringsAsFactors = FALSE,
                                                   check.names = FALSE
                                                   )
                             )

plate_numbers <- as.integer(vapply(layout_96_mat_list, function(x) x[11, 1], ""))
stopifnot(identical(plate_numbers, seq_along(layout_96_mat_list)))

layout_96_mat_list <- lapply(layout_96_mat_list, function(x) {
  result_mat <- x[1:8, 2:ncol(x)]
  dimnames(result_mat) <- NULL
  return(result_mat)
})


sheets_vec <- readxl::excel_sheets(long_layout_96_path)
long_96_mat_list <- lapply(sheets_vec,
                           function(x) as.matrix(suppressMessages(read_excel(long_layout_96_path, sheet = x, col_names = FALSE)),
                                                 stringsAsFactors = FALSE,
                                                 check.names = FALSE
                                                 )
                           )


converted_96_mat_list <- lapply(long_96_mat_list, function(x) {
  result_mat <- matrix(x[, 2], nrow = 8, ncol = 12)
  return(result_mat)
})

stopifnot(identical(layout_96_mat_list, converted_96_mat_list))




# Define functions --------------------------------------------------------

PlateMappings <- function() {

  main_plate_96_column_indices <- list(
    "yellow" = rep(2:12, times = 7),
    "red"    = rep(2:12, times = 7),
    "green"  = rep(1:11, times = 7)
  )
  main_plate_384_column_indices <- list(
    "yellow" = rep(seq(from = 3, to = 23, by = 2), times = 7)
  )
  main_plate_384_column_indices <- list(
    "yellow" = main_plate_384_column_indices[["yellow"]],
    "red"    = main_plate_384_column_indices[["yellow"]],
    "green"  = main_plate_384_column_indices[["yellow"]] - 1L
  )

  main_plate_96_row_indices <- list(
    "yellow" = rep(1:7, each = 11),
    "red"    = rep(2:8, each = 11),
    "green"  = rep(1:7, each = 11)
  )
  main_plate_384_row_indices <- list(
    "yellow" = rep(seq(from = 2, to = 14, by = 2), each = 11)
  )
  main_plate_384_row_indices <- list(
    "yellow" = main_plate_384_row_indices[["yellow"]],
    "red"    = main_plate_384_row_indices[["yellow"]] + 1L,
    "green"  = main_plate_384_row_indices[["yellow"]]
  )

  coords_96wp_mat <- matrix(seq_len(96), nrow = 8, ncol = 12, byrow = TRUE)
  coords_384wp_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)

  vertical_384_columns   <- seq(from = 3, to = 17, by = 2)
  horizontal_384_columns <- seq(from = 3, to = 23, by = 2)

  edge_plate_96_vertical <- list(
    "yellow" = coords_96wp_mat[, 1],
    "red"    = coords_96wp_mat[, 1],
    "green"  = coords_96wp_mat[, 12]
  )

  edge_plate_384_vertical <- list(
    "yellow" = coords_384wp_mat[3, vertical_384_columns],
    "red"    = coords_384wp_mat[5, vertical_384_columns],
    "green"  = coords_384wp_mat[7, vertical_384_columns]
  )

  edge_plate_96_horizontal <- list(
    "yellow" = coords_96wp_mat[8, 2:12],
    "red"    = coords_96wp_mat[1, 2:12],
    "green"  = coords_96wp_mat[8, 1:11]
  )

  edge_plate_384_horizontal <- list(
    "yellow" = coords_384wp_mat[2, horizontal_384_columns],
    "red"    = coords_384wp_mat[4, horizontal_384_columns],
    "green"  = coords_384wp_mat[6, horizontal_384_columns]
  )

  colors_df_list <- lapply(c("yellow", "red", "green"), function(x) {

    wells_96_vec <- vapply(seq_along(main_plate_96_column_indices[[x]]), function(y) {
      coords_96wp_mat[main_plate_96_row_indices[[x]][[y]], main_plate_96_column_indices[[x]][[y]]]
    }, integer(1))

    wells_384_vec <- vapply(seq_along(main_plate_384_column_indices[[x]]), function(y) {
      coords_384wp_mat[main_plate_384_row_indices[[x]][[y]], main_plate_384_column_indices[[x]][[y]]]
    }, integer(1))

    yellow_df <- data.frame(
      "Color" = x,
      rbind.data.frame(
        data.frame("Plate_type_384"  = "main",
                   "Well_number_96"  = wells_96_vec,
                   "Well_number_384" = wells_384_vec,
                   stringsAsFactors  = FALSE
                   ),
        data.frame("Plate_type_384"  = "edge",
                   "Well_number_96"  = c(edge_plate_96_horizontal[[x]],
                                         edge_plate_96_vertical[[x]]
                                         ),
                   "Well_number_384" = c(edge_plate_384_horizontal[[x]],
                                         edge_plate_384_vertical[[x]]
                                         ),
                   stringsAsFactors  = FALSE
                   ),
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind.data.frame, c(colors_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  return(results_df)
}




IntegratePlateLayouts <- function(layouts_96wp_mat_list, mappings_df) {

  ## Specify well coordinates
  coords_384wp_mat <- matrix(seq_len(384), nrow = 16, ncol = 24, byrow = TRUE)
  alphanumeric_384wp_vec <- unlist(lapply(LETTERS[1:16], function(x) {
    lapply(1:24, function(y) {
      paste0(x, y)
    })
  }))
  alphanumeric_96wp_vec <- unlist(lapply(LETTERS[1:8], function(x) {
    lapply(1:12, function(y) {
      paste0(x, y)
    })
  }))

  ## Define control well positions (empty wells = uninfected)
  control_rows <- seq(3, 15, by = 2)
  pos_controls_384_vec <- c(coords_384wp_mat[control_rows, 4],
                            coords_384wp_mat[control_rows, 20]
                            )
  NT_controls_384_vec <- c(coords_384wp_mat[control_rows, 2],
                           coords_384wp_mat[control_rows, 22]
                           )


  ## Initialize the results data frame
  results_df <- data.frame(
    "Plate_number_384" = rep(1:3, each = 384),
    "Well_number_384"  = rep(seq_len(384), times = 3),
    "Well_coords_384"  = NA,
    "Well_ID_384"      = NA,
    "Target_ID"        = NA,
    "Plate_number_96"  = NA,
    "Well_number_96"   = NA,
    "Well_coords_96"   = NA,
    "Well_ID_96"       = NA,
    stringsAsFactors   = FALSE
  )

  ## Define variables for the loop
  are_used_96wp <- rep(FALSE, 6)
  are_used_96wp[seq_along(layouts_96wp_mat_list)] <- TRUE
  use_colors <- rep(c("yellow", "red", "green"), 2)
  main_plates <- c(rep(1L, 3), rep(2L, 3))


  ## Populate the data frame with the target genes
  for (i in seq_along(layouts_96wp_mat_list)) {

    if (!(are_used_96wp[[i]])) {
      next
    }
    this_96wp_targets <- as.vector(t(layouts_96wp_mat_list[[i]]))

    are_this_main_384wp <- results_df[["Plate_number_384"]] == main_plates[[i]]
    are_this_edge_384wp <- results_df[["Plate_number_384"]] == 3
    are_this_96wp <- mappings_df[, "Color"] == use_colors[[i]]

    are_main <- are_this_96wp & (mappings_df[, "Plate_type_384"] == "main")
    main_indices_384 <- mappings_df[are_main, "Well_number_384"]
    main_indices_96  <- mappings_df[are_main, "Well_number_96"]
    results_df[are_this_main_384wp, "Well_number_96"][main_indices_384] <- main_indices_96
    results_df[are_this_main_384wp, "Target_ID"][main_indices_384] <- this_96wp_targets[main_indices_96]
    results_df[are_this_main_384wp, "Plate_number_96"][main_indices_384] <- i

    are_edge <- are_this_96wp & (mappings_df[, "Plate_type_384"] == "edge")
    edge_indices_384 <- mappings_df[are_edge, "Well_number_384"]
    edge_indices_96  <- mappings_df[are_edge, "Well_number_96"]
    if (i %in% 4:6) {
      edge_indices_384 <- edge_indices_384 + as.integer((384 / 2)) # these 96-well plates go in the bottom half of the edge 384-well plate
    }
    results_df[are_this_edge_384wp, "Well_number_96"][edge_indices_384] <- edge_indices_96
    results_df[are_this_edge_384wp, "Target_ID"][edge_indices_384] <- this_96wp_targets[edge_indices_96]
    results_df[are_this_edge_384wp, "Plate_number_96"][edge_indices_384] <- i
  }


  ## Create well coordinates and combined well IDs
  results_df[["Well_coords_96"]] <- alphanumeric_96wp_vec[results_df[["Well_number_96"]]]
  results_df[["Well_coords_384"]] <- alphanumeric_384wp_vec[results_df[["Well_number_384"]]]

  results_df[["Well_ID_96"]] <- ifelse(is.na(results_df[["Plate_number_96"]]),
                                       NA,
                                       paste0("Plate",
                                              formatC(results_df[["Plate_number_96"]],
                                                      width = 2, flag = "0"
                                                      ),
                                              "_",
                                              results_df[["Well_coords_96"]]
                                              )
                                       )

  results_df[["Well_ID_384"]] <- paste0("Plate",
                                        formatC(results_df[["Plate_number_384"]],
                                                width = 2, flag = "0"
                                                ),
                                        "_",
                                        results_df[["Well_coords_384"]]
                                        )

  ## Add control wells
  for (i in 1:3) {
    are_this_384wp <- results_df[["Plate_number_384"]] == i
    results_df[are_this_384wp, "Target_ID"][pos_controls_384_vec] <- "Pos. control"
    results_df[are_this_384wp, "Target_ID"][NT_controls_384_vec]  <- "Own NT control"

    well_IDs_vec <- results_df[are_this_384wp, "Well_ID_96"]
    stopifnot(!(any(duplicated(well_IDs_vec[!(is.na(well_IDs_vec))]))))
  }

  return(results_df)
}




# Integrate plate mappings ------------------------------------------------

plate_mappings_df <- PlateMappings()

plate_groups <- c(rep(1, 4), rep(2:4, each = 6))
stopifnot(length(plate_groups) == 22)

plate_group_splits <- split(seq_along(layout_96_mat_list), plate_groups)

integrated_df_list <- lapply(plate_group_splits, function(x) {
  IntegratePlateLayouts(layout_96_mat_list[x], plate_mappings_df)
})




# Correct the plate numbers -----------------------------------------------

for (i in (seq_len(length(integrated_df_list) - 1L) + 1L)) {
  integrated_df_list[[i]][["Plate_number_96"]] <- integrated_df_list[[i]][["Plate_number_96"]] +
                                                  max(integrated_df_list[[i - 1]][["Plate_number_96"]], na.rm = TRUE)
  integrated_df_list[[i]][["Plate_number_384"]] <- integrated_df_list[[i]][["Plate_number_384"]] +
                                                   max(integrated_df_list[[i - 1]][["Plate_number_384"]])
}

integrated_df <- do.call(rbind.data.frame, c(integrated_df_list, stringsAsFactors = FALSE, make.row.names = FALSE))

integrated_df[["Well_ID_96"]] <- ifelse(is.na(integrated_df[["Plate_number_96"]]),
                                        NA,
                                        paste0("Plate",
                                               formatC(integrated_df[["Plate_number_96"]],
                                                       width = 2, flag = "0"
                                                       ),
                                               "_",
                                               integrated_df[["Well_coords_96"]]
                                               )
                                        )


integrated_df[["Plate_number_384"]] <- as.character(as.roman(1:12)[integrated_df[["Plate_number_384"]]])
integrated_df[["Well_ID_384"]] <- paste0("Plate",
                                         integrated_df[["Plate_number_384"]],
                                         "_",
                                         integrated_df[["Well_coords_384"]]
                                         )



# Interpret flags ---------------------------------------------------------

target_splits <- strsplit(integrated_df[, "Target_ID"], "_", fixed = TRUE)

have_PC <- grepl("_PC_", integrated_df[, "Target_ID"], fixed = TRUE)
have_empty <- grepl("empty", integrated_df[, "Target_ID"], ignore.case = TRUE)


target_splits <- lapply(target_splits, function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    if (!(is.na(suppressWarnings(as.integer(x[[2]]))))) {
      results_vec <- paste0(x[[1]], "_", x[[2]])
      if (length(x) > 2) {
        results_vec <- c(results_vec, x[3:length(x)])
      }
      return(results_vec)
    } else {
      return(x)
    }
  }
})

target_splits <- lapply(target_splits, function(x) {
  if (length(x) >= 3) {
    return(c(x[[1]], paste0(x[2:length(x)], collapse = "_")))
  } else {
    return(x)
  }
})

have_flag <- lengths(target_splits) > 1
flags_vec <- rep(NA_character_, nrow(integrated_df))
flags_vec[have_flag] <- sapply(target_splits[have_flag], "[[", 2)

integrated_df[have_flag & (grepl("rep", flags_vec, ignore.case = TRUE)), ]

cleaned_flags_vec <- vapply(flags_vec, function(x) {
  if (is.na(x)) {
    NA_character_
  } else if (grepl("PC_", x, fixed = TRUE)) {
    "Tuebingen pos. control"
  } else if (grepl("Rep", x, fixed = TRUE)) {
    "Rep"
  } else if (identical(x, "Empty")) {
    "Empty"
  } else {
    NA_character_
  }
}, "")

integrated_df[, "Target_ID_stripped"] <- sapply(target_splits, "[[", 1)
integrated_df[, "Target_flag"] <- cleaned_flags_vec

not_from_96wp <- is.na(integrated_df[, "Plate_number_96"])
integrated_df[not_from_96wp, "Target_flag"] <- integrated_df[not_from_96wp, "Target_ID"]




# Perform manual corrections ----------------------------------------------

integrated_df[grepl("Sox11", integrated_df[, "Target_ID"], ignore.case = TRUE), ]
integrated_df[integrated_df[, "Target_ID_stripped"] %in% "Sox11", "Target_ID_stripped"] <- "SOX11"

integrated_df[grepl("REST", integrated_df[, "Target_ID"], ignore.case = TRUE), ]
integrated_df[integrated_df[, "Target_ID_stripped"] %in% "REST", "Target_ID_stripped"] <- "REST_2"

integrated_df[grepl("NFE2L2", integrated_df[, "Target_ID"], ignore.case = TRUE), ]
integrated_df[integrated_df[, "Target_ID_stripped"] %in% "NFE2L2", "Target_ID_stripped"] <- "NFE2L2_1"




# Create analogous TSS IDs for the 4sg library ----------------------------

are_to_use <- (!(CRISPRa_4sg_df[, "Is_obsolete"] %in% "Yes")) &
                (CRISPRa_4sg_df[, "Sublibrary_4sg"] %in% c("Misc / changed TFs", "Transcription factors"))

use_4sg_df <- CRISPRa_4sg_df[are_to_use, ]

stopifnot(identical(unique(as.integer(table(use_4sg_df[, "AltTSS_ID"]))), 4L))

TSS_splits <- lapply(split(use_4sg_df[, "AltTSS_ID"], use_4sg_df[, "Entrez_ID"]), unique)
num_TSSs_vec <- lengths(TSS_splits)

TSS_vec <- rep(NA, nrow(use_4sg_df))
for (entrez_ID in names(TSS_splits)[num_TSSs_vec >= 2]) {
  num_TSSs <- num_TSSs_vec[[entrez_ID]]
  are_this_gene <- use_4sg_df[["Entrez_ID"]] == entrez_ID
  TSS_vec[are_this_gene] <- rep(seq_len(num_TSSs), each = 4)
}

IDs_4sg_vec <- paste0(use_4sg_df[["Gene_symbol"]],
                      ifelse(is.na(TSS_vec), "", paste0("_", TSS_vec))
                      )



# Add plate and well numbers ----------------------------------------------

plate_string_splits <- strsplit(use_4sg_df[, "Plate_string"], "_", fixed = TRUE)
plates_vec <- sapply(plate_string_splits, "[[", 2)
plates_vec <- paste0("HA_", sub("tf", "", plates_vec, fixed = TRUE))
use_4sg_df[, "Library_plate"] <- plates_vec
use_4sg_df[, "Library_coords"] <- ConvertWellNumbers(use_4sg_df[, "Well_number"])



# Examine unexpected target IDs -------------------------------------------

stopifnot(all(IDs_4sg_vec %in% integrated_df[, "Target_ID_stripped"]))

known_own_IDs <- setdiff(integrated_df[not_from_96wp, "Target_ID"], NA)

other_IDs <- setdiff(integrated_df[, "Target_ID_stripped"], IDs_4sg_vec)
other_IDs <- other_IDs[!(is.na(other_IDs))]
other_IDs <- setdiff(other_IDs, known_own_IDs)

are_empty <- grepl("empty", integrated_df[, "Target_ID_stripped"], ignore.case = TRUE)
stopifnot(all(is.na(integrated_df[are_empty, "Target_flag"])))
integrated_df[are_empty, "Target_flag"] <- "Empty"

are_scrambled <- grepl("scrambled", integrated_df[, "Target_ID_stripped"], ignore.case = TRUE)
stopifnot(all(is.na(integrated_df[are_scrambled, "Target_flag"])))
integrated_df[are_scrambled, "Target_flag"] <- "Scrambled"

are_mCherry <- grepl("mCherry", integrated_df[, "Target_ID_stripped"], ignore.case = TRUE)
stopifnot(all(is.na(integrated_df[are_mCherry, "Target_flag"])))
integrated_df[are_mCherry, "Target_flag"] <- "mCherry"


integrated_df[grepl("NFE2L2", integrated_df[, "Target_ID"], ignore.case = TRUE), ]
integrated_df[grepl("Sox11", integrated_df[, "Target_ID"], ignore.case = TRUE), ]

integrated_df[grepl("HA[0-9]{1,3}_", integrated_df[, "Target_ID"], ignore.case = TRUE), ]
integrated_df[grepl("scrambled", integrated_df[, "Target_ID"], ignore.case = TRUE), ]

## Questions:
# - Do the "NFE2L2" positive controls target TSS1 or TSS2? Is my assumption correct that they target TSS1?
# - Is my assumption correct that the repetition well labelled "REST" targets TSS2 (since the well for TSS2 was empty?)
# - What are the following wells: "HA5_B11_Rep", "HA5_P12_Rep", "HA4_K8_Rep", "HA4_F9_Rep" and "HA5_N18_Rep"?
# - Do "scrambled" wells correspond to non-targeting controls?



# Integrate 4sg library data ----------------------------------------------

matches_vec <- match(integrated_df[, "Target_ID_stripped"], IDs_4sg_vec)

use_columns <- c("Gene_symbol", "Entrez_ID", "TSS_ID", "Is_main_TSS",
                 "Library_plate", "Library_coords"
                 )

add_df <- use_4sg_df[matches_vec, use_columns]

columns_96wp <- grep("_96", names(integrated_df), value = TRUE, fixed = TRUE)
first_columns <- setdiff(names(integrated_df), c("Target_ID_stripped", columns_96wp))

layout_df <- data.frame(integrated_df[, first_columns],
                        use_4sg_df[matches_vec, use_columns],
                        integrated_df[, columns_96wp],
                        row.names = NULL,
                        stringsAsFactors = FALSE
                        )



# Define positive and negative controls -----------------------------------

layout_df[, "Is_NT_ctrl"]  <- layout_df[, "Target_flag"] %in% c("Own NT control", "Scrambled")
layout_df[, "Is_pos_ctrl"] <- layout_df[, "Target_flag"] %in% "Pos. control"




# Arrange into 384-well plate format --------------------------------------

plates_fac <- factor(integrated_df[, "Plate_number_384"],
                     levels = unique(integrated_df[, "Plate_number_384"])
                     )
layouts_384_mat_list <- lapply(split(integrated_df, plates_fac), function(x) {
  matrix(x[["Target_ID"]], nrow = 16, ncol = 24, byrow = TRUE)
})



# Export data -------------------------------------------------------------

write.csv(layout_df,
          file      = file.path(output_dir, "Plate_layout", "CRISPRa TF - all plates.csv"),
          row.names = FALSE,
          na        = ""
          )

for (i in seq_along(layouts_384_mat_list)) {
  plate_name <- paste0(formatC(i, width = 2, flag = "0"),
                       ") CRISPRa TF - Plate_", as.roman(i), ".csv"
                       )
  export_mat <- layouts_384_mat_list[[i]]
  colnames(export_mat) <- as.character(seq_len(ncol(export_mat)))
  write.csv(export_mat,
            file      = file.path(output_dir, "Plate_layout", "12_plates", plate_name),
            row.names = LETTERS[seq_len(nrow(export_mat))],
            na        = ""
            )
}



# Save data ---------------------------------------------------------------

save(layout_df, file = file.path(general_rdata_dir, "01_convert_plate_layouts.RData"))



