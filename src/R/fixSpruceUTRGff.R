#' ---
#' title: "Fixing the spruce gff including UTR"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(magrittr)
  library(readr)
})

#' * Path
gff3_file <- here("reference/Picea-abies/v1.1/gff3/Pabies1.1-gene.artificial-100bp-UTRs.sorted.gff3.gz")

#' # Data
#' * Gff3
df <- read_tsv(gff3_file,
               comment="#",col_names=FALSE,show_col_types=FALSE)

#' * Contig length
ctig <- df %>% distinct(X1) %>% unlist(use.names=FALSE)
clength <- read_tsv(here("reference/Picea-abies/v1.0/fasta/GenomeAssemblies/Pabies01-genome-collapsed-for-STAR.fa.fai"),
                    col_names=FALSE,col_types=cols_only(X1=col_character(),X2=col_double())) %>% filter(X1 %in% ctig)

#' # Cleanup
#' * Remove the introns
df %<>% filter(X3!="intron")

#' * Extract the gene IDs
geneIDs <- df %>% 
  filter(X3=="gene") %>% 
  transmute(ID=sub("ID=","",X9)) %>% 
  unlist(use.names=FALSE)

#' * Fix the genes coordinates
genePos <- apply(
  bind_cols(
    df %>% 
      filter(X3=="five_prime_UTR") %>% 
      transmute(x=ifelse(X7=="+",X4,X5)),
    df %>% 
      filter(X3=="three_prime_UTR") %>% 
      transmute(y=ifelse(X7=="+",X5,X4))),1,range)

sel <- df$X3 == "gene"
df[sel,"X4"] <- genePos[1,]
df[sel,"X5"] <- genePos[2,]

#' * Add an mRNA
df %<>% bind_rows(df %>% 
                    filter(sel) %>% 
                    mutate(X3="mRNA",X9=paste0("ID=",geneIDs))) %>% 
  arrange(X1,X4,desc(X5))

#' * Fix the CDS
sel <- df$X3 == "CDS"
df[sel,] <- df %>% 
  filter(sel) %>% 
  mutate(X9=sub("\\.exon\\.\\d+","",X9))

#' * Fix the UTRs
sel <- df$X3=="five_prime_UTR"
df[sel,] <- df %>% 
  filter(sel) %>% 
  mutate(
    X4=ifelse(X7=="-",X4+1,X4),
    X5=ifelse(X7=="+",X5-1,X5),
    X9=paste0(X9,".5pUTR;Parent=",geneIDs))

sel <- df$X3=="three_prime_UTR"
df[sel,] <- df %>% 
  filter(sel) %>% 
  mutate(
    X4=ifelse(X7=="+",X4+1,X4),
    X5=ifelse(X7=="-",X5-1,X5),
    X9=paste0(X9,".3pUTR;Parent=",geneIDs))

#' * Update the ID
sel <- df$X3 != "gene"
df[sel,] <- df %>% 
  filter(sel) %>% 
  mutate(X9=paste0(X9,".1"))

#' * Remove wrong entries
df %<>% filter(! X5 < X4)

#' * There are no loci outside the scaffolds
stopifnot((df %>% left_join(clength,by="X1") %>% transmute(X2.y < X5) %>% sum()) == 0)

#' # Export
#' ## GFF3
gff_out <- sub("\\.gz","",gff3_file)
write("##gff-version 3",gff_out)
write_tsv(df,file=gff_out,append=TRUE)

#' ## GTF
#' Transform to gtf
gtf <- df %>% filter(X3 != "gene") %>% 
  mutate(X9=sub("\\..*;Parent=","; transcript_id ",sub("ID=","gene_id ",X9)))
sel <- gtf$X3 == "mRNA"
gtf[sel, ] <- gtf %>% filter(sel) %>% mutate(X9=paste(sub("\\.1", ";", X9),sub("gene","transcript",X9)))

gtf_out <- gsub("gff3","gtf",gff_out)
write_tsv(gtf,file=gtf_out,append=TRUE)

#' # Test 
#' ## GFF3
txdb <- GenomicFeatures::makeTxDbFromGFF(gff_out)

#' ## GTF
annot <- riboWaltz::create_annotation(gtfpath = gtf_out)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
