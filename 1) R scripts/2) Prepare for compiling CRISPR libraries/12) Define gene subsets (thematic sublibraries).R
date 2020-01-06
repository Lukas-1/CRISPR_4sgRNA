### 5th January 2019 ###



# Import packages and source code -----------------------------------------

general_functions_directory <- "~/CRISPR/1) R scripts/1) R functions"
source(file.path(general_functions_directory, "02) Translating between Entrez IDs and gene symbols.R"))





# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
CRISPR_input_directory  <- file.path(CRISPR_root_directory, "2) Input data")
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")





# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "06) Collect Entrez IDs from various sources.RData"))
load(file.path(general_RData_directory, "08) Compile a list of human transcription factors - all_TF_df.RData"))






# Read in data ------------------------------------------------------------

# Downloaded from https://kampmannlab.ucsf.edu/next-generation-shrna-libraries-datasets
# direct link: https://kampmannlab.ucsf.edu/sites/kampmannlab.ucsf.edu/files/Dataset_S4_human_sublibraries.txt
Kampmann_lib_df <- read.table(file.path(CRISPR_input_directory, "Sublibraries",
                                        "Kampmann, Horlbeck, Weissmann - PNAS 2015",
                                        "Dataset_S4_human_sublibraries.txt"
                                        ),
                              sep = "\t", quote = "", stringsAsFactors = FALSE,
                              header = FALSE, row.names = NULL
                              )







# Process the thematic sublibraries of Kampmann et al. --------------------

colnames(Kampmann_lib_df) <- c("Sublibrary_ID", "Symbol")

shRNA_sublibrary_map <- c(
  "H01" = "Apoptosis",
  "H02" = "Cancer",
  "H03" = "Drug Targets",
  "H04" = "Gene Expression",
  "H05" = "Kinase/Phosphatase 1",
  "H06" = "Kinase/Phosphatase 2",
  "H07" = "Membrane Proteins",
  "H08" = "Mitochondria",
  "H09" = "Motility",
  "H10" = "Other Cancer",
  "H11" = "Stress/Proteostasis",
  "H12" = "Trafficking",
  "H13" = "Unassigned"
)

Kampmann_lib_df[, "Sublibrary"] <- shRNA_sublibrary_map[Kampmann_lib_df[, "Sublibrary_ID"]]

Kampmann_mapped_df <- MapToEntrezs(symbols_vec = Kampmann_lib_df[, "Symbol"])

Kampmann_lib_df <- data.frame(
  Kampmann_lib_df[, c("Sublibrary_ID", "Sublibrary")],
  Kampmann_mapped_df[, c("Entrez_ID", "Gene_symbol", "Original_symbol", "Entrez_source")],
  stringsAsFactors = FALSE,
  row.names = NULL
)















