#!/bin/r env

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-wd", "--wd", 
                    help="Path to input to files and output files")

parser$add_argument("-fd", "--fd", 
                    help="Path to function files and reference files")

parser$add_argument("-vcf", "--vcf", 
                    help="name of the vcf file")

parser$add_argument("-chr", "--chr", 
                    help="Chromosom to filter")

parser$add_argument("-o", "--out", 
                    help="Prefix of output files")

parser$add_argument("-n", "--n", 
                    help="Number of iterations to split the vcf file")

parser$add_argument("-thres", "--thres", 
                    help="Minimun read depth to call an allele")

parser$add_argument("-p", "--p", 
                    help="Parasite: Pf or Pv")

args = parser$parse_args()

# Working directories

wd = args$wd
setwd(wd)

fd = args$fd

# First filter arguments
vcf = args$vcf

# Select chromosome
chr = as.integer(args$chr)

parasite = args$p

## Split chromosomes----

if(parasite == 'Pf'){
  
  chromosomes =  c(
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
  )
  
}else if(parasite == 'Pv'){
  
  chromosomes =  c(
    'PvP01_01_v1',
    'PvP01_02_v1',
    'PvP01_03_v1',
    'PvP01_04_v1',
    'PvP01_05_v1',
    'PvP01_06_v1',
    'PvP01_07_v1',
    'PvP01_08_v1',
    'PvP01_09_v1',
    'PvP01_10_v1',
    'PvP01_11_v1',
    'PvP01_12_v1',
    'PvP01_13_v1',
    'PvP01_14_v1',
    'PvP01_API_v1',
    'PvP01_MIT_v1'
  )
  
}

# Output
output = paste0(args$out, '_',chromosomes[chr])

n = as.integer(args$n)
threshold = as.numeric(args$thres)

imagename = paste0(output, '.RData')

# Check packages and functions

source(file.path(fd,'load_libraries.R'))
source(file.path(fd,'functions.R'))

run_vcftools(vcf = vcf,
             bash_file = paste0(chromosomes[chr], '_run_vcf.sh'),
             out = output,
             chr = chromosomes[chr],
             recode = TRUE,
             recode_INFO_all = TRUE)

### Load vcf data----

vcf_object = load_vcf(vcf = paste0(output, '.recode.vcf'))

### Convert to rGenome class----

rGenome_object =  vcf2rGenome(vcf = vcf_object, n = n, threshold = threshold)

assign(paste0(chromosomes[chr], '_vcf_object'), vcf_object)
assign(paste0(chromosomes[chr], '_rGenome_object'), rGenome_object)

# rm(list = ls()[!grepl(
#   paste0(
#     paste0(
#       chromosomes[chr], '_rGenome_object'),
#     '|',
#     paste0(
#       chromosomes[chr], '_vcf_object')),ls())])


save.image(file = imagename)











