### 15th January 2020 ###



# Define functions --------------------------------------------------------

AppendIDsWithoutEntrezs <- function(entrez_IDs_list, CRISPR_df) {
  have_no_entrez <- is.na(CRISPR_df[, "Entrez_ID"]) &
                    (CRISPR_df[, "Is_control"] == "No")
  IDs_list <- entrez_IDs_list
  IDs_list[[length(IDs_list)]] <- c(IDs_list[[length(IDs_list)]],
                                    unique(extended_CRISPRko_df[have_no_entrez, "Combined_ID"])
                                    )
  return(IDs_list)
}
