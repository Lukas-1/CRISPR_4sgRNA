### 3rd August 2021 ###



# Define functions --------------------------------------------------------

AddFilterPassed <- function(contamin_df) {

  ## Integrate data on thresholds passed

  threshold_vec <- ifelse(contamin_df[["ZMW"]] %in% ccs7_df_list[["individual_reads_df"]][["ZMW"]],
                          "CCS7",
                          ifelse(contamin_df[["ZMW"]] %in% ccs5_df_list[["individual_reads_df"]][["ZMW"]],
                                 "CCS5",
                                 ifelse(contamin_df[["ZMW"]] %in% ccs3_df_list[["individual_reads_df"]][["ZMW"]],
                                        "CCS3",
                                        "None"
                                        )
                                 )
                          )
  contamin_df[["Filter_passed"]] <- threshold_vec


  ## Re-order the columns

  move_after_column <- "Are_cross_plate"
  insert_column <- "Filter_passed"
  move_after_index <- which(names(contamin_df) == move_after_column)
  are_before <- seq_along(contamin_df) < move_after_index
  are_after <- seq_along(contamin_df) > move_after_index
  are_current <- names(contamin_df) == insert_column
  columns_reordered <- c(names(contamin_df)[are_before],
                         move_after_column, insert_column,
                         names(contamin_df)[are_after & !(are_current)]
                         )
  contamin_df <- contamin_df[, columns_reordered]

  return(contamin_df)
}




MakeExportContamDf <- function(contamin_df) {

  results_df <- contamin_df[order(contamin_df[["Are_cross_plate"]], decreasing = TRUE), ]
  results_df <- results_df[results_df[["Filter_passed"]] != "None", ]
  row.names(results_df) <- NULL

  results_df[["Are_cross_plate"]] <- ifelse(results_df[["Are_cross_plate"]], "Yes", "No")
  names(results_df)[names(results_df) == "Are_cross_plate"] <- "Cross_plate"
  return(results_df)
}




