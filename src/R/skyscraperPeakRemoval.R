suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

prob <- 0.99

ccov <- read_tsv(here("data/codon_coverage.dms.txt")) %>% 
  column_to_rownames("X1") %>% group_by(transcript) %>% 
  mutate(max=qpois(prob,mean(pseudoalignments))) %>% 
  mutate(filter=pseudoalignments > max) %>% 
  mutate(corrected_pseudocounts=ifelse(pseudoalignments>max,-1,pseudoalignments))

ccov <- read_tsv(here("data/codon_coverage.dms.txt")) %>% 
  column_to_rownames("X1") %>% group_by(transcript) %>% 
  mutate(corrected_pseudocounts=ifelse(pseudoalignments>qpois(prob,mean(pseudoalignments)),-1,pseudoalignments))
