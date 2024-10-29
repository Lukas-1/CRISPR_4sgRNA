### 2024-09-22



# Load packages -----------------------------------------------------------

library("shiny")



# Define folder paths -----------------------------------------------------

project_dir   <- file.path("~", "CRISPR_4sgRNA", "shiny")
libraries_dir <- file.path(project_dir, "2_input", "our_CRISPR_libraries")
rdata_dir     <- file.path(project_dir, "3_RData")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "2_tidy_CRISPR_libraries.RData"))




# Define labels -----------------------------------------------------------

library_list <- list(
  "T.gonfio (CRISPR activation)" = "CRISPRa",
  "T.spiezzo (CRISPR knockout)" = "CRISPRko"
)

wells_vec <- paste0(rep(LETTERS[1:16], each = 24), rep(1:24, times = 16))




# Define functions --------------------------------------------------------

with_tooltip <- function(value, tooltip) {
  # from https://glin.github.io/reactable/articles/cookbook/cookbook.html#tooltips-1
  tags$abbr(style = "text-decoration: underline; text-decoration-style: dotted; cursor: help",
            title = tooltip,
            value
            )
}

small_width <- 70
small_header_font <- "0.7em"

GenerateColumnsList <- function(input_df) {

  columns_list <- list(

    "Sublibrary_4sg" = reactable::colDef(
      "Sublibrary",
      minWidth = 200
    ),
    "Plate_ID" = reactable::colDef(
      "Plate",
      minWidth = small_width
    ),
    "Well_coords" = reactable::colDef(
      "Well",
      minWidth = small_width
    ),
    "Is_obsolete" = reactable::colDef(
      header = with_tooltip("Obsolete plasmid?", "Some transcription factor plasmids were re-designed, and plasmids for these genes appear twice.\nThe obsolete version is on plates 1 to 5, and the final version on plate 5+.")
    ),
    "Entrez_ID" = reactable::colDef(
      "NCBI gene ID",
      align = "left"
    ),
    "Other_target_Entrez_IDs" = reactable::colDef(
      "Other targets (NCBI IDs)"
    ),
    "Other_Entrez_IDs_4sg" = reactable::colDef(
      header = with_tooltip("Other targets of the plasmid", "The 4 sgRNAs acting in concert will disrupt genes that lie between the first and last cut site.")
    ),
    "Gene_symbol" = reactable::colDef(
      "Gene symbol",
      style = list(fontStyle = "italic", fontWeight = "bold")
    ),
    "Other_target_symbols" = reactable::colDef(
      "Other targets (gene symbols)",
      style = list(fontStyle = "italic")
    ),
    "Other_symbols_4sg" = reactable::colDef(
      header = with_tooltip("Other targets of the plasmid", "The 4 sgRNAs acting in concert will disrupt genes that lie between the first and last cut site.")
    ),
    "Original_symbol" = reactable::colDef(
      show = FALSE,
      searchable = TRUE
    ),
    "TSS_ID" = reactable::colDef(
      header = with_tooltip("TSS ID", "TSS = Transcription Start Site\nSome genes have multiple TSSs, which can result in different protein products.")
    ),
    "Is_main_TSS" = reactable::colDef(
      header = with_tooltip("Is main TSS?", "Single: There is only one")
    ),
    "Rank" = reactable::colDef(
      "sgRNA number",
      align = "center"
    ),
    "Num_overlaps" = reactable::colDef(
      header = with_tooltip("Overlaps", '4 sgRNAs spaced at least 50 bp apart were chosen whenever possible.\n"1|50 bp" indicates that one pair of sgRNAs had cut locations within 50 bp of each other.')
    ),
    "Source" = reactable::colDef(
      "Source library"
    ),
    "sgRNA_sequence" = reactable::colDef(
      "Sequence",
      style = list(fontFamily = "monospace", fontSize = "0.9em")
    ),
    "PAM" = reactable::colDef(
      header = with_tooltip("PAM", "Protospacer Adjacent Motif"),
      style = list(fontFamily = "monospace")
    ),
    "Calabrese_rank" = reactable::colDef(
      "Calabrese rank"
    ),
    "GPP_rank" = reactable::colDef(
      "CRISPick rank"
    ),
    "hCRISPRa_v2_rank" = reactable::colDef(
      "hCRISPRa-v2 rank"
    ),
    "Predicted_score"= reactable::colDef(
      "hCRISPRa-v2 predicted score"
    ),
    "Empirical_score"= reactable::colDef(
      "hCRISPRa-v2 empirical score"
    ),
    "Chromosome"= reactable::colDef(
      "Chromo-some",
      minWidth = small_width
    ),
    "Strand"= reactable::colDef(
      "Strand",
      minWidth = small_width
    ),
    "Cut_location" = reactable::colDef(
      "Cut location"
    ),
    "Pseudo_cut_location" = reactable::colDef(
      header = with_tooltip('"Cut location"', "The position where the Cas9 enzyme would cut, if it were active.\nThe hg38 reference genome is used.")
    ),
    "Distance_from_TSS" = reactable::colDef(
      "Distance from TSS"
    ),
    "GuideScan_efficiency" = reactable::colDef(
      header = with_tooltip("GuideScan efficiency score", "Based on the Rule Set 2 scores by Doench et al., 2016")
    ),
    "CRISPOR_Doench_efficacy" = reactable::colDef(
      header = with_tooltip("CRISPOR efficacy score", "Also based on the Rule Set 2 scores by Doench et al., 2016")
    ),
    "CRISPOR_Graf_status" = reactable::colDef(
      header = with_tooltip("Graf criteria", "Does the sgRNA fulfil the two criteria of Graf et al., 2019?")
    ),
    "GuideScan_specificity" = reactable::colDef(
      header = with_tooltip("GuideScan specificity score", "Includes off-target sites with up to 3 base mismatches.\nSee Perez et al., 2017.")
    ),
    "CRISPOR_3MM_specificity" = reactable::colDef(
      header = with_tooltip("3MM CRISPOR specificity score", "Includes off-target sites with up to 3 base mismatches.\nSee Concordet et al., 2018.")
    ),
    "CRISPOR_4MM_specificity" = reactable::colDef(
      header = with_tooltip("4MM CRISPOR specificity score", "Includes off-target sites with up to 4 base mismatches.\nSee Concordet et al., 2018.")
    ),
    "CRISPOR_CFD_specificity"= reactable::colDef(
      header = with_tooltip("CRISPOR CFD specificity", "Original CRISPOR specificity score, as included in the output from the CRISPOR webtool.")
    ),
    "Num_0MM" = reactable::colDef(
      header = with_tooltip("Number of 0MM locations", "Number of potential sgRNA binding sites with zero base mismatches.")
    ),
    "Num_1MM" = reactable::colDef(
      "Number of 1MM locations"
    ),
    "GuideScan_Num_2MM" = reactable::colDef(
      header = with_tooltip("GuideScan 2MM", "Number of potential binding sites with 2 mismatched bases, according to the GuideScan tool."),
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "GuideScan_Num_3MM" = reactable::colDef(
      "GuideScan 3MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "CRISPOR_Num_0MM" = reactable::colDef(
      "CRISPOR 0MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "CRISPOR_Num_1MM" = reactable::colDef(
      "CRISPOR 1MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "CRISPOR_Num_2MM" = reactable::colDef(
      "CRISPOR 2MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "CRISPOR_Num_3MM" = reactable::colDef(
      "CRISPOR 3MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "CRISPOR_Num_4MM" = reactable::colDef(
      "CRISPOR 4MM",
      minWidth = small_width,
      headerStyle = list(fontSize = small_header_font)
    ),
    "all22_SNP_IDs_vcf" = reactable::colDef(
      header = with_tooltip("SNP rsIDs", "Identifiers of SNPs with frequencies of at least 0.1%\nwithin the 22 bases of the sgRNA or PAM."),
      style = list(fontSize = "0.8em")
    ),
    "all22_SNP_AF_max_Kaviar" = reactable::colDef(
      header = with_tooltip("SNP frequency", "Frequency of the most common SNP\n(according to the kaviar database; https://db.systemsbiology.net/kaviar/)"),
      format = reactable::colFormat("percent" = TRUE)
    ),
    "Locations_0MM"= reactable::colDef(
      header = with_tooltip("0MM locations", "Perfect-match binding sites for the sgRNA.\nThere is usually only one, but there can be multiple for genes with closely related paralogs."),
      style = list(fontSize = "0.6em")
    ),
    "Locations_1MM"= reactable::colDef(
      header = with_tooltip("1MM locations", "Potential binding sites (featuring a PAM) with a single-base mismatch."),
      style = list(fontSize = "0.6em")
    ),
    "Color" = reactable::colDef(
      show = FALSE
    )
  )

  columns_list <- columns_list[names(columns_list) %in% names(input_df)]
  return(columns_list)
}



# Define UI ---------------------------------------------------------------

ui <- fluidPage(
  radioButtons("library_selection", "Select the CRISPR library", library_list),
  reactable::reactableOutput("table")
)





# Define server logic -----------------------------------------------------

server <- function(input, output, session) {

  GetDisplayDf <- reactive({
    if (input$library_selection == "CRISPRa") {
      results_df <- CRISPRa_sgRNA_df
      results_df[, "Plate_ID"] <- paste0("HA", results_df[, "Plate_ID"])
      names(results_df)[[which(names(results_df) == "Cut_location")]] <- "Pseudo_cut_location"
    } else if (input$library_selection == "CRISPRko") {
      results_df <- CRISPRko_sgRNA_df
      results_df[, "Plate_ID"] <- paste0("HO", results_df[, "Plate_ID"])
    }

    results_df[[2]] <- results_df[["Plate_ID"]]
    results_df[["Plate_ID"]] <- NULL
    names(results_df)[[2]] <- "Plate_ID"

    results_df[, "Well_number"] <- wells_vec[results_df[, "Well_number"]]
    names(results_df)[[which(names(results_df) == "Well_number")]] <- "Well_coords"

    remove_columns <- c("Plasmid_ID", "Sequences_1MM", "GuideScan_offtarget_category",
                        "Exon_number", "Transcript_ID", "Genomic_sequence_ID"
                        )

    results_df <- results_df[, !(names(results_df) %in% remove_columns)]
    return(results_df)
  })

  GetColumnsList <- reactive(GenerateColumnsList(GetDisplayDf()))


  output$table <- reactable::renderReactable({
    use_df <- GetDisplayDf()
    reactable::reactable(
      use_df,
      columns = GenerateColumnsList(use_df),
      defaultPageSize = 40,
      rowStyle = function(index) {
        if (use_df[index, "Color"] == 1) {
          list(backgroundColor = "#FFF2CC")
        }
      }
    )
  })
}

shinyApp(ui, server)



