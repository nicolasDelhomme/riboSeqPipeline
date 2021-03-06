#' ---
#' title: "Potra V2 RiboWaltz pipeline"
#' author: "Bernard Wessels & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Installation
#' Dependencies: devtools RiboWaltz - should only be done once
#' ```{r install,eval=FALSE}
#' install.packages("devtools")
#' devtools::install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE, build_vignettes = TRUE)
#' ```

#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(here)
  library(GenomicFeatures)
  library(riboWaltz)
})

#' # Process
#' ## Annotation
#' ```{r gtf,eval=FALSE,echo=FALSE}
#' # This block will be hidden when creating the report and contains examples for Ath
#' 
#' # GTF file
#' gtf_file <- here("reference/Arabidopsis-thaliana/TAIR10/gff/TAIR10_GFF3_genes_transposons.gtf")
#' # Annotation TxDb file
#' txdbfile <- makeTxDbFromGFF(gtf_file, format=c("auto"), dataSource=NA, organism=NA, taxonomyId=NA, circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL, miRBaseBuild=NA)
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' TxDb.Athaliana.BioMart.plantsmart28
#'
#' Other species
#' Annot <- create_annotation(gtfpath = here("reference/Populus-tremula_X_Populus-tremuloides/v1.0/gtf/Potrx01-genome.gtf.gz"), txdb = NULL)
#' Annot <- create_annotation(gtfpath =  here("reference/Arabidopsis-thaliana/ARAPORT11/gtf/Araport11_GFF3_genes_transposons.201606.gtf"), txdb = NULL)
#' ```

Annot <- create_annotation(gtfpath = here("reference/Populus-tremula/v2.2/gtf/Potra02_genes.gtf.gz"), 
                                          txdb = NULL)


#' ## Data
#' * Import and read BAM files
reads_list <- bamtolist(bamfolder = here("data/Bernard_Results_Samples1-6Okt20/kallisto/Sample1_S1_L001_R1_001_sortmerna_trimmomatic"), 
                        annotation = Annot)

#' * Selection of read lengths
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom",
                               length_filter_vector = 26:34)

#' ```{r period, eval=FALSE, echo=FALSE}
#' #' * Periodicity threshold mode
#' filtered_list <- length_filter(data = reads_list, length_filter_mode = "periodicity",
#'                              periodicity_threshold = 60)
#' ```

#' * P_Site Offset
psite_offset <- psite(filtered_list, flanking = 6, start = TRUE, extremity = "auto",
                      plot = TRUE, plot_dir = here("data/Bernard_Results_Samples1-6Okt20/P_Site_images/"), 
                      plot_format = "png", cl = 99)
reads_psite_list <- psite_info(filtered_list, psite_offset)

#' * Codon Coverage
codon_coverage_example <- codon_coverage(reads_psite_list, Annot, psite = FALSE)

#' ```{r cov, eval=FALSE, echo=FALSE}
#' #' * Transcript Region Coverage
#' TranscriptID <- subset(codon_coverage_example, transcript %in% c("Potra2n1c3551.1"))
#' TranscriptID$region2 = factor(TranscriptID$region, levels=c('5utr','cds','3utr'))
#' plot <- ggplot(TranscriptID, aes(x = from_cds_start, y = pseudoalignments, colour = region, fill=region)) + geom_bar(stat = "identity") + theme_bw()
#' plot
#' ```

#' * CDS coverage
cds_coverage_example <- cds_coverage(reads_psite_list, Annot)

#' ```{r codon cov, eval=FALSE, echo=FALSE}
#' #' * Codon coverage sorted
#' codon_coverage_variance <- dcast(codon_coverage_example,  transcript ~ region, function(x){(sum(x>1)/(sum(x<=1)*var(x)))}, value.var = "pseudoalignments")
#' filtered_codon_cov <- codon_coverage_variance[!is.na(codon_coverage_variance$cds), ]
#' ```

#' # Visualisation
#' ## Heatmap
example_ends_heatmap <- rends_heat(filtered_list, Annot, sample = "pseudoalignments", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]

#' ## P-sites per region
example_psite_region <- region_psite(reads_psite_list, Annot, sample = "pseudoalignments")
example_psite_region[["plot"]]

example_frames <- frame_psite(reads_psite_list, sample = "pseudoalignments", region = "all")
example_frames[["plot"]]

example_metaprofile <- metaprofile_psite(reads_psite_list, Annot, sample = "pseudoalignments",
                                             length_range = 28, utr5l = 20, cdsl = 40,
                                             utr3l = 20, plot_title = "auto") 
example_metaprofile[["plot_pseudoalignments"]]

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
