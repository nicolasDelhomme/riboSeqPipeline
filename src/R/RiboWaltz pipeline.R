#Install
install.packages("devtools")
#Loading
library(devtools)
#Installing RiboWaltz
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE)
#Please note: to install riboWaltz generating the vignette replace the last command with:
install_github("LabTranslationalArchitectomics/riboWaltz", dependencies = TRUE, build_vignettes = TRUE)
#To load riboWaltz run
library(riboWaltz)

library(GenomicFeatures)

#GTF file
#gtf_file <- "/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/gff/TAIR10_GFF3_genes_transposons.gtf"

#Annotation TxDb file
#txdbfile <- makeTxDbFromGFF(gtf_file, format=c("auto"), dataSource=NA, organism=NA, taxonomyId=NA, circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL, miRBaseBuild=NA)

#library(TxDb.Athaliana.BioMart.plantsmart28)
#TxDb.Athaliana.BioMart.plantsmart28

#Creating annotation file
#Annot <- create_annotation(gtfpath = "/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/gff/TAIR10_GFF3_genes_transposons.gtf", txdb = NULL)
Annot <- create_annotation(gtfpath = "/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline/riboseq/data/Araport11_GFF3_genes_transposons.201606.gtf", txdb = NULL)


#Import and read BAM files
reads_list <- bamtolist(bamfolder = "/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline/riboseq/results20191008/kallisto/Galaxy300_sortmerna_trimmomatic/", annotation = Annot)
#Selection of read lengths
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom",
                               length_filter_vector = 24:34)
#Periodicity threshold mode
filtered_list <- length_filter(data = reads_list, length_filter_mode = "periodicity",
                               periodicity_threshold = 60)

#P_Site Offset
psite_offset <- psite(filtered_list, flanking = 6, start = TRUE, extremity = "auto",
                      plot = TRUE, plot_dir = "/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline/riboseq/results", plot_format = "png", cl = 99)
reads_psite_list <- psite_info(filtered_list, psite_offset)
#Codon Coverage
codon_coverage_example <- codon_coverage(reads_psite_list, Annot, psite = FALSE)
#CDS coverage
cds_coverage_example <- cds_coverage(reads_psite_list, Annot)

#read length
example_length_dist <- rlength_distr(reads_list, sample = "pseudoalignments", cl = 99)
example_length_dist[["plot"]]

#Heatmap
example_ends_heatmap <- rends_heat(reads_list, Annot, sample = "pseudoalignments", cl = 85,
                                   utr5l = 25, cdsl = 40, utr3l = 25)
example_ends_heatmap[["plot"]]

#P-sites per region

example_psite_region <- region_psite(reads_psite_list, Annot, sample = "pseudoalignments")
example_psite_region[["plot"]]

#Trinucleotide Periodicity
example_frames_stratified <- frame_psite_length(reads_psite_list, sample = "pseudoalignments",
region = "all", cl = 90)
example_frames_stratified[["plot"]]

example_frames <- frame_psite(reads_psite_list, sample = "pseudoalignments", region = "all")
example_frames[["plot"]]

example_metaprofile <- metaprofile_psite(reads_psite_list, Annot, sample = "pseudoalignments",
                                             length_range = 28, utr5l = 20, cdsl = 40,
                                             utr3l = 20, plot_title = "auto") 
example_metaprofile[["plot"]]


