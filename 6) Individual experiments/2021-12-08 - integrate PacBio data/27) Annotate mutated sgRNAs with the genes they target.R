### 5th January 2022 ###



# Import packages and source code -----------------------------------------

library("devEMF")

CRISPR_root_directory       <- "~/CRISPR_4sgRNA"
general_functions_directory <- file.path(CRISPR_root_directory, "1) R scripts/1) R functions")
experiments_directory       <- file.path(CRISPR_root_directory, "6) Individual experiments")
plate1_directory            <- file.path(experiments_directory, "2020-08-29 - PacBio - first 384-well plate")
R_functions_directory       <- file.path(plate1_directory, "1) R functions")

source(file.path(general_functions_directory, "06) Helper functions for genomic ranges.R")) # For TruncateLongEntries
source(file.path(general_functions_directory, "22) Generating statistics and plots for CRISPR libraries.R"))
source(file.path(general_functions_directory, "30) Finding overlapping genes and nearby TSSs.R"))

source(file.path(R_functions_directory, "24) Finding unintended targets of mutated gRNAs.R"))




# Define folder paths -----------------------------------------------------

library_RData_directory  <- file.path(CRISPR_root_directory, "3) RData files")
general_RData_directory  <- file.path(library_RData_directory, "1) General")

s2rI_directory           <- file.path(experiments_directory, "2021-12-08 - integrate PacBio data")
s2rI_R_objects_directory <- file.path(s2rI_directory, "3) R objects")
tables_directory         <- file.path(s2rI_directory, "5) Output", "Tables",  "Targets of mutated gRNAs")
figures_directory        <- file.path(s2rI_directory, "5) Output", "Figures")
donut_plots_directory    <- file.path(figures_directory, "Targets of mutated gRNAs")
manuscript_directory     <- file.path(figures_directory, "Manuscript", "Fig. S4", "Individual plots")



# Load data ---------------------------------------------------------------

load(file.path(general_RData_directory, "20) Compile all relevant TSSs for each gene.RData"))
load(file.path(general_RData_directory, "21) Assemble data frames of gene, transcript, exon and CDS coordinates.RData"))

load(file.path(s2rI_R_objects_directory, "03) Import and process sgRNA sequences.RData"))
load(file.path(s2rI_R_objects_directory, "26) Annotate mutated sgRNAs with any perfect matches in the genome.RData"))




# Add gene information to the mutations_df data frame ---------------------

mutations_df <- AddGeneInfo(mutations_df, library_df)




# Prepare the mutations data frames ---------------------------------------

TSS_19bp_mut_list <- AnnotateMutations(mutations_df, "CRISPRa",
                                       FindNearbyTSSs, all_TSS_df
                                       )
TSS_20bp_mut_list <- AnnotateMutations(mutations_df, "CRISPRa",
                                       FindNearbyTSSs, all_TSS_df,
                                       use_20bp = TRUE
                                       )


CRISPRko_19bp_mut_list <- AnnotateMutations(mutations_df, "CRISPRko",
                                            FindOverlappingGenes,
                                            CDS_or_exon_locations_df
                                            )
CRISPRko_20bp_mut_list <- AnnotateMutations(mutations_df, "CRISPRko",
                                            FindOverlappingGenes,
                                            CDS_or_exon_locations_df,
                                            use_20bp = TRUE
                                            )



# Export data on individual reads -----------------------------------------

ExportMutatedDf(TSS_19bp_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__19bp__CRISPRa"
                )
ExportMutatedDf(CRISPRko_19bp_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__19bp__CRISPRko"
                )

ExportMutatedDf(TSS_20bp_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__20bp__CRISPRa"
                )
ExportMutatedDf(CRISPRko_20bp_mut_list[["annotated_df"]],
                "Targets_of_mutated_gRNAs__20bp__CRISPRko"
                )



# Draw doughnut/bar plots -------------------------------------------------

pdf(file = file.path(donut_plots_directory, "Donut charts - 19bp - targets of mutated gRNAs.pdf"),
    width = 6, height = 4
    )
MutationsDonutBar(CRISPRko_19bp_mut_list[["gRNA_numbers"]],
                  "CRISPRo library (target: CDS or exon)"
                  )

MutationsDonutBar(TSS_19bp_mut_list[["gRNA_numbers"]],
                  "CRISPRa library (target: within 1000 bp of the TSS)"
                  )
dev.off()



pdf(file = file.path(donut_plots_directory, "Donut charts - 20bp - targets of mutated gRNAs.pdf"),
    width = 6, height = 4
    )
MutationsDonutBar(CRISPRko_20bp_mut_list[["gRNA_numbers"]],
                  "CRISPRo library (target: CDS or exon)"
                  )

MutationsDonutBar(TSS_20bp_mut_list[["gRNA_numbers"]],
                  "CRISPRa library (target: within 1000 bp of the TSS)"
                  )
dev.off()




# Draw doughnut/bar plots for the manuscript ------------------------------

ManuscriptMutationsDonutBar <- function(gRNA_numbers, main_title, donut_label) {

  donut_labels <- c(
    "Lacks a perfect-match site\nin the human genome",
    "Target sites do not\naffect any gene",
    "Shares target genes with\nthe unmutated gRNA",
    "Affects a new\noff-target gene"
  )
  donut_colors <- c("#DDDDDD", "#88CCEE", "#332288", "#AA4499")

  DonutBars(counts_vec         = SubtractFollowing(gRNA_numbers),
            use_colors         = donut_colors,
            use_labels         = donut_labels,
            use_title          = main_title,
            show_axis          = TRUE,
            title_line         = 0.4,
            donut_label        = donut_label,
            donut_text_size    = 0.8,
            donut_radius       = 0.4,
            donut_inner_radius = 0.43,
            donut_x_mid        = 0.73,
            donut_y_mid        = 0.325,
            space              = 0.8,
            use_mai            = c(0, 1.4, 0.4, 0.15),
            use_omi            = c(0, 0, 0, 0),
            use_line_height    = 1.1,
            bar_label_line     = 0.5,
            bar_text_size      = 1,
            side_text_size     = 1,
            bar_text_font      = 1,
            donut_text_font    = 1,
            title_text_size    = 1,
            text_dark_color    = "black",
            side_space_ratio   = 0.6
            )
  return(invisible(NULL))
}



pdf(file.path(manuscript_directory, paste0("E) Doughnut plot - CRISPRa library.pdf")),
    width = 3.4, height = 2
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ManuscriptMutationsDonutBar(TSS_20bp_mut_list[["gRNA_numbers"]],
                  "CRISPRa library (target: within 1000 bp of the TSS)",
                  "Mutated\nCRISPRa\ngRNAs"
                  )
dev.off()



pdf(file.path(manuscript_directory, paste0("F) Doughnut plot - CRISPRo library.pdf")),
    width = 3.4, height = 2
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ManuscriptMutationsDonutBar(CRISPRko_20bp_mut_list[["gRNA_numbers"]],
                  "CRISPRo library (target: CDS or exon)",
                  "Mutated\nCRISPRo\ngRNAs"
                  )
dev.off()




emf(file.path(manuscript_directory, paste0("E) Doughnut plot - CRISPRa library.emf")),
    width = 3.4, height = 2, emfPlus = FALSE, coordDPI = 3000
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ManuscriptMutationsDonutBar(TSS_20bp_mut_list[["gRNA_numbers"]], "",
                  "Mutated\nT.gonfio\ngRNAs"
                  )
dev.off()



emf(file.path(manuscript_directory, paste0("F) Doughnut plot - CRISPRo library.emf")),
    width = 3.4, height = 2, emfPlus = FALSE, coordDPI = 3000
    )
par(cex = manuscript_cex, lwd = manuscript_lwd)
ManuscriptMutationsDonutBar(CRISPRko_20bp_mut_list[["gRNA_numbers"]], "",
                  "Mutated\nT.spiezzo\ngRNAs"
                  )
dev.off()





# Save data ---------------------------------------------------------------

save(list = c("TSS_19bp_mut_list", "TSS_20bp_mut_list",
              "CRISPRko_19bp_mut_list", "CRISPRko_20bp_mut_list"
              ),
     file = file.path(s2rI_R_objects_directory, "27) Annotate mutated sgRNAs with the genes they target.RData")
     )


