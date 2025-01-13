### 2024-09-22


# Define folder paths -----------------------------------------------------

project_dir <- file.path("~", "CRISPR_4sgRNA", "shiny")
rdata_dir   <- file.path(project_dir, "3_RData")



# Load data ---------------------------------------------------------------

load(file.path(rdata_dir, "3_filter_annotation_data.RData"))



# Filter GO terms by evidence ---------------------------------------------

evidence_codes_list <- list(
  c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP"),
  c("HTP", "HDA", "HMP", "HGI", "HEP"),
  c("IBA", "IBD", "IKR", "IRD"),
  c("ISS", "ISO", "ISA", "ISM", "IGC", "RCA"),
  c("TAS", "NAS", "IC", "ND"),
  c("IEA")
)




table(GO_df[, "Evidence"] %in% high_confidence_GO_codes)
