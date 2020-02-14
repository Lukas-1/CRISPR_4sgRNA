### 13th February 2020 ###



# Import packages and source code -----------------------------------------

library("readxl")

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))
source(file.path(general_functions_directory, "24) Assigning genes to sublibraries.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")

CRISPRa_datasets_directory       <- file.path(CRISPR_input_directory, "CRISPR libraries", "CRISPRa")
CRISPRa_Horlbeck2016_path        <- file.path(CRISPRa_datasets_directory, "Horlbeck, Kampmann, Weissman - eLife 2016")
CRISPRa_Horlbeck2016_sgRNAs_path <- file.path(CRISPRa_Horlbeck2016_path, "2016 - Compact and highly active next-generation libraries - Table S5.xlsx")




# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))
load(file.path(general_RData_directory, "10) Compile genes that constitute the secretome - secretome_df.RData"))





# Read in data ------------------------------------------------------------

hCRISPRa_v2_df <- data.frame(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, skip = 7)[-1, ], stringsAsFactors = FALSE)
names(hCRISPRa_v2_df) <- names(read_excel(CRISPRa_Horlbeck2016_sgRNAs_path, n_max = 1))





# Select gene-sublibrary mappings from the hCRISPRa-v2 library ------------

hC2_sublibrary_df <- unique(hCRISPRa_v2_df[, c("gene", "Sublibrary")])
colnames(hC2_sublibrary_df)[colnames(hC2_sublibrary_df) == "Sublibrary"] <- "Sublibrary_code"

hC2_sublibrary_df <- hC2_sublibrary_df[hC2_sublibrary_df[["gene"]] != "negative_control", ]

if (any(duplicated(hC2_sublibrary_df[["gene"]]))) {
  stop("Duplicated symbol-sublibrary combination found!")
}

hC2_sublibrary_df[["Sublibrary"]] <- hCRISPRa_v2_sublibrary_map[hC2_sublibrary_df[["Sublibrary_code"]]]





# Map gene symbols to Entrez IDs ------------------------------------------

hC2_sublibrary_df <- data.frame(MapToEntrezs(symbols_vec = hC2_sublibrary_df[["gene"]]),
                                hC2_sublibrary_df[, c("Sublibrary_code", "Sublibrary")],
                                stringsAsFactors = FALSE
                                )
hC2_sublibrary_df <- hC2_sublibrary_df[, colnames(hC2_sublibrary_df) != "Original_entrez"]





# Resolve duplicated Entrez IDs, where possible ---------------------------

num_occurrences_vec <- table(hC2_sublibrary_df[["Entrez_ID"]])[hC2_sublibrary_df[["Entrez_ID"]]]
multiplicates_df <- hC2_sublibrary_df[(num_occurrences_vec >= 2) %in% TRUE, ]
multiplicates_df <- multiplicates_df[order(GetMinEntrez(multiplicates_df[["Entrez_ID"]])), ]

unicates_df <- hC2_sublibrary_df[(num_occurrences_vec == 1) %in% TRUE, ]

resolved_df_list <- by(multiplicates_df,
                       factor(multiplicates_df[["Entrez_ID"]], levels = unique(multiplicates_df[["Entrez_ID"]])),
                       function(x) {
                         are_original_symbol <- x[["Original_symbol"]] == ""
                         if (any(are_original_symbol)) {
                           x[are_original_symbol, ]
                         } else {
                           x
                         }
                       }, simplify = FALSE
                       )

hC2_sublibrary_df <- do.call(rbind.data.frame, c(list(unicates_df), resolved_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
hC2_sublibrary_df <- hC2_sublibrary_df[order(GetMinEntrez(hC2_sublibrary_df[["Entrez_ID"]])), ]
rownames(hC2_sublibrary_df) <- NULL




# Define the Entrez IDs from prioritized sub-libraries --------------------

TF_entrez_IDs <- TidyEntrezs(all_TF_df[["Entrez_ID"]][all_TF_df[["Is_TF"]] == "Yes"])

secretome_entrez_IDs <- TidyEntrezs(secretome_df[["Entrez_ID"]])




# Create a combined list --------------------------------------------------

sublibrary_list <- c(
  list(
    "Transcription factors" = TF_entrez_IDs,
    "Secretome" = secretome_entrez_IDs
  ),
  tapply(hC2_sublibrary_df[["Entrez_ID"]],
         factor(hC2_sublibrary_df[["Sublibrary"]], levels = hCRISPRa_v2_sublibrary_map),
         function(x) TidyEntrezs(unlist(strsplit(x, ", ", fixed = TRUE))),
         simplify = FALSE
         )
)




# Assign protein-coding genes to each of the sub-libraries ----------------

are_in_library_mat <- do.call(cbind, lapply(sublibrary_list, function(x) collected_entrez_IDs %in% x))

assignment_vec <- apply(are_in_library_mat, 1, function(x) if (!(any(x))) "None" else colnames(are_in_library_mat)[[which(x)[[1]]]])

sublibrary_df <- data.frame("Entrez_ID" = collected_entrez_IDs,
                            "Sublibrary" = factor(assignment_vec, levels = c(colnames(are_in_library_mat), "None")),
                            stringsAsFactors = FALSE
                            )




# Save data ---------------------------------------------------------------

save(list = "sublibrary_df",
     file = file.path(general_RData_directory, "11) Divide the remaining genes into sublibraries according to hCRISPRa-v2 - sublibrary_df.RData")
     )














