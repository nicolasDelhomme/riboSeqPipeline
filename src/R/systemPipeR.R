# libraries
#library(here)
#library(tidyverse)
library(systemPipeR)
library(systemPipeRdata)


# working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/riboseq-pipeline")

# create the folder structure
if (!file.exists("systemPipeRIBOseq.R")) {
  genWorkenvir(workflow = "riboseq")  
}

setwd("riboseq")

# targets
#targets <- read_tsv("targets2.txt",comment = "#")
targets <- read.delim("targets3.txt",comment = "#")

# fastq_trimming
args <- systemArgs(sysma="param/trim.param", mytargets="targets3.txt")
fctpath <- system.file("extdata", "custom_Fct.R", package="systemPipeR")

source(fctpath)
iterTrim <- ".iterTrimbatch1(fq, pattern='CTGTAGGCACCATCAAT', internalmatch=TRUE, minpatternlength=9,
                             Nnumber=1, polyhomo=50, minreadlength=16, maxreadlength=101)"
preprocessReads(args=args, Fct=iterTrim, batchsize=100000, overwrite=TRUE, compress=TRUE)
writeTargetsout(x=args, file="targets_trim.txt", overwrite=TRUE)

# fastq_report
args <- systemArgs(sysma=NULL, mytargets="targets_trim.txt")
fqlist <- seeFastq(fastq=infile1(args), batchsize=100000, klength=8)
png("./results/fastqReport.png", height=18, width=4*length(fqlist), units="in", res=72)
seeFastqPlot(fqlist)
dev.off()

# Alignment
## Bowtie
# We need to modify the params
# open param/bowtieSE.param
# in the Terminal run:
# deactive the module line
# change other parameters, e.g. the reference
args <- systemArgs(sysma="param/bowtieSE.param", mytargets="targets_trim.txt")

# Set up the resources needed for Bowtie
# walltime: DD-HH:MM:SS - absolute limit, jobs dies if it takes longer
resources <- list(walltime = "4:00:00", ntasks = 1, 
                  ncpus = cores(args), 
                  memory = "16G",
                  out.file="results/bowtie.out",
                  err.file="results/bowtie.err")

# Register the job in the queue
# Njobs = number of files to process 
reg <- clusterRun(args, conffile = ".BatchJobs.R", Njobs = nrow(targets), 
                  template = "slurm.tmpl", runid = "bowtie", 
                  resourceList = resources)

# wait for the job to be finished before proceeding
# you can also check the queue in the terminal buy running `squeue`
waitForJobs(reg = reg)

# check the status of the job
getStatus(reg = reg)


## Kallisto
# We need to modify the params
# open param/bowtieSE.param
# in the Terminal run:
# deactive the module line
# change other parameters, e.g. the reference
# adapt the -l and -s option (fragment length and its sd - we are being lenient)
args <- systemArgs(sysma="param/kallisto.param", mytargets="targets_trim.txt")
args@sysargs <- sub("\\.bam","\\.tsv",args@sysargs)


# Set up the resources needed for Bowtie
# walltime: DD-HH:MM:SS - absolute limit, jobs dies if it takes longer
resources <- list(walltime = "4:00:00", ntasks = 1, 
                  ncpus = cores(args), 
                  memory = "16G",
                  out.file="results/kallisto.out",
                  err.file="results/kallisto.err")

# Register the job in the queue
# Njobs = number of files to process 
reg <- clusterRun(args, conffile = ".BatchJobs.R", Njobs = nrow(targets), 
                  template = "slurm.tmpl", runid = "bowtie", 
                  resourceList = resources)

# wait for the job to be finished before proceeding
# you can also check the queue in the terminal buy running `squeue`
waitForJobs(reg = reg)

# check the status of the job
getStatus(reg = reg)


