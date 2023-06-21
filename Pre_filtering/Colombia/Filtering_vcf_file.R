#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-wd", "--wd", 
                    help="Path to input to files and output files")

parser$add_argument("-fd", "--fd", 
                    help="Path to function files and reference files")

parser$add_argument("-gzvcf", "--gzvcf", 
                    help="name of the gzvcf file")

parser$add_argument("-o", "--out", 
                    help="Prefix of output files")

parser$add_argument("-rkeep", "--keep_regexp", 
                    help="Regular expression to identify samples of interest")

parser$add_argument("-ebed", "--exclude_bed", 
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-gff", "--ref_gff", 
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-n", "--n", 
                    help="Number of iterations to split the vcf file")

parser$add_argument("-thres", "--thres", 
                    help="Minimun read depth to call an allele")

args = parser$parse_args()

wd = args$wd
setwd(wd)

fd = args$fd

# First filter arguments
gzvcf = args$gzvcf
foutput = paste0(args$out, '_FirstFilter')
keep_regexp = args$keep_regexp
exclude_bed = file.path(fd, args$exclude_bed)

# Second filter arguments
svcf = paste0(args$out, '_FirstFilter.recode.vcf')
ref_gff_file = file.path(fd, args$ref_gff)
soutput = paste0(args$out, '_SecondFilter')

# Load VCF to R and generate rGenome object

tvcf = paste0(args$out, '_SecondFilter.recode.vcf')
n = as.integer(args$n)
threshold = as.numeric(args$thres)

imagename = paste0(args$out, '_1st_2nd_filters.RData')


# Check packages and functions----

source(file.path(fd,'load_libraries.R'))
source(file.path(fd,'functions.R'))
sourceCpp(file.path(fd,'Rcpp_functions.cpp'))
# Section 1: Pre filtering of trusted sites ----

## Step 1: First filter of non core regions, and samples of interest----

if(!(file.exists('step1.check'))){
  start.time = Sys.time()
  run_vcftools(gzvcf = gzvcf,
               out = foutput,
               keep_regexp = keep_regexp,
               exclude_bed = exclude_bed,
               remove_filtered_all = TRUE,
               non_ref_ac_any = 1,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  system('rm samples.indv')
  
  end.time = Sys.time()

  write.table(paste0('Step 1 done,', ' it started at ', start.time, ', and it finished at ', end.time), 'step1.check', row.names = FALSE, quote = FALSE, col.names = FALSE)
  
}


## Step 2: Second filter of coding regions----

if(!(file.exists('step2.check'))){
  start.time = Sys.time()
  
  ref_gff = ape::read.gff(ref_gff_file)
  
  ref_gff = ref_gff[grepl('gene', ref_gff$type)&
                            !grepl('^Transfer',ref_gff$seqid),]
  
  ref_gff = cbind(ref_gff, as.data.frame(t(sapply(1:nrow(ref_gff), function(gene){
    attributes = strsplit(ref_gff[gene,][['attributes']], ';')[[1]]
    c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
      gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
  }))))
  
  coding_regions = ref_gff[,c('seqid', 'start', 'end')]
  
  coding_regions = coding_regions[order(coding_regions$seqid),][order(coding_regions$start),]
  
  rownames(coding_regions) = 1:nrow(coding_regions)
  
  write.table(coding_regions, 'coding_regions.bed', sep = '\t', quote = FALSE, row.names = FALSE)
  
  run_vcftools(vcf = svcf,
               out = soutput,
               bed = 'coding_regions.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  system('rm coding_regions.bed')
  
  ### Load vcf data----
  
  vcf_object = load_vcf(vcf = tvcf)
  
  ### Convert to rGenome class 
  
  rGenome_object = vcf2rGenome(vcf = vcf_object, n = n, threshold = threshold)
  
  end.time = Sys.time()
  
  rm(list = ls()[-grep('rGenome_object|vcf_object',ls())])
  
  save.image(imagename)
  write.table(paste0('Step 2 done,', ' it started at ', start.time, ', and it finished at ', end.time), 'step2.check', row.names = FALSE, quote = FALSE, col.names = FALSE)
}












