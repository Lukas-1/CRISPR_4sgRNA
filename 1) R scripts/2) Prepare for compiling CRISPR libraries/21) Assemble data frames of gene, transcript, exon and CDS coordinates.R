### 4th August 2020 ###



# Import packages and source code -----------------------------------------

library("TxDb.Hsapiens.UCSC.hg38.knownGene")




# Define folder paths -----------------------------------------------------

CRISPR_root_directory   <- "~/CRISPR"
general_RData_directory <- file.path(CRISPR_root_directory, "3) RData files", "1) General")






# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "07) Compile TSS (transcription start site) data.RData"))
load(file.path(general_RData_directory, "17) Compile the information on gene type.RData"))
load(file.path(general_RData_directory, "20) Process the annotations from GENCODE.RData"))








# Create data frames using TxDb.Hsapiens.UCSC.hg38.knownGene --------------

TxDb_genes_df <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                     single.strand.genes.only = FALSE
                                     ))
TxDb_transcripts_df <- as.data.frame(transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene))
TxDb_exons_df <- as.data.frame(exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))
TxDb_CDS_df <- as.data.frame(cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))

stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_transcripts_df[["group_name"]]))
                    ))
stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_transcripts_df[["group_name"]]))
                    ))
stopifnot(identical(sort(unique(TxDb_genes_df[["group_name"]])),
                    sort(unique(TxDb_exons_df[["group_name"]]))
                    ))








# Create a merged data frame of gene coordinates --------------------------






groo1 <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene)
groo2 <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene)
groo3 <- cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene)

foo1 <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene)
foo2 <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")
foo3 <- cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")



length(unique(names(human_genes_GRanges)))
length(unique(as.data.frame(groo1)[["group_name"]]))
length(unique(as.data.frame(groo2)[["group_name"]]))


length(unique(as.data.frame(groo1[["group_name"]])))



id2name(TxDb.Hsapiens.UCSC.hg38.knownGene, feature.type = "tx")
id2name(TxDb.Hsapiens.UCSC.hg38.knownGene, feature.type = "exon")
id2name(TxDb.Hsapiens.UCSC.hg38.knownGene, feature.type = "cds")


















# Load data ---------------------------------------------------------------
