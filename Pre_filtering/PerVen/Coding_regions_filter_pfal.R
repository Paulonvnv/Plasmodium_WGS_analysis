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

# Second filter arguments
ref_gff_file = file.path(fd, args$ref_gff)
output = paste0(args$out, '_CodingRegionsOnly')

# Load VCF to R and generate rGenome object

tvcf = paste0(args$out, '_CodingRegionsOnly.recode.vcf')
n = as.integer(args$n)
threshold = as.numeric(args$thres)

imagename = paste0(args$out, '_CodingRegionsOnly.RData')


# Check packages and functions----

source(file.path(fd,'load_libraries.R'))
source(file.path(fd,'functions.R'))
#sourceCpp(file.path(fd,'Rcpp_functions.cpp'))
# Section 1: Pre filtering coding regions ----


## Step 2: Second filter of coding regions----

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

run_vcftools(vcf = gzvcf,
             out = output,
             bed = 'coding_regions.bed',
             recode = TRUE,
             recode_INFO_all = TRUE)

system('rm coding_regions.bed')


for(chrom in c(
  'Pf3D7_01_v3',
  'Pf3D7_02_v3',
  'Pf3D7_03_v3',
  'Pf3D7_04_v3',
  'Pf3D7_05_v3',
  'Pf3D7_06_v3',
  'Pf3D7_07_v3',
  'Pf3D7_08_v3',
  'Pf3D7_09_v3',
  'Pf3D7_10_v3',
  'Pf3D7_11_v3',
  'Pf3D7_12_v3',
  'Pf3D7_13_v3',
  'Pf3D7_14_v3',
  'Pf3D7_API_v3',
  'Pf3D7_MIT_v3'
)){
  
  run_vcftools(vcf = paste0(output, '.recode.vcf'),
               out = paste0(chrom, output),
               chr = chrom,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
}


### Load vcf data----

temp_vcf_object = load_vcf(vcf = paste0(chrom, output, '.recode.vcf'))

### Convert to rGenome class 

temp_rGenome_object = vcf2rGenome(vcf = temp_vcf_object, n = n, threshold = threshold)

end.time = Sys.time()

rm(list = ls()[-grep('rGenome_object|vcf_object',ls())])

save.image(imagename)
write.table(paste0('Step 2 done,', ' it started at ', start.time, ', and it finished at ', end.time), 'step2.check', row.names = FALSE, quote = FALSE, col.names = FALSE)














