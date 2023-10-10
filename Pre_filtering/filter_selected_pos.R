#!/bin/r env
wd = '/gsap/garage-protistvector/PanSAmer'
setwd(wd)
source('/gsap/garage-protistvector/ColombiaData/Pviv/AnalysisTools/load_libraries.R')
source('/gsap/garage-protistvector/ColombiaData/Pviv/AnalysisTools/functions.R')

# Panama----

path_to_panama = '/gsap/garage-protistvector/lbuyon/171206_Pv.SWGA.vcf'

run_vcftools(vcf = path_to_panama,
             bash_file = 'run_vcf_Panama.sh',
             out = 'Panama_Pviv_Filtered',
             positions = 'selected_pos.bed',
             recode = TRUE,
             recode_INFO_all = TRUE)

Panama_vcf = load_vcf(vcf = 'Panama_Pviv_Filtered.recode.vcf')
Panama_rGenome =  vcf2rGenome(vcf = Panama_vcf, n = 500, threshold = 5)
# system('rm Panama_Pviv_Filtered.recode.vcf')

# Colombia----

path_to_colombia = '/gsap/garage-protistvector/ColombiaData/Pviv/DataGeneration/ColombiaGates_Pviv.JointCall.filtered.combined.snpeff.vcf.gz'

run_vcftools(gzvcf = path_to_colombia,
             bash_file = 'run_vcf_Colombia.sh',
             out = 'ColombiaGates_Pviv_Filtered',
             positions = 'selected_pos.bed',
             recode = TRUE,
             recode_INFO_all = TRUE)

Colombia_vcf = load_vcf(vcf = 'ColombiaGates_Pviv_Filtered.recode.vcf')
Colombia_rGenome =  vcf2rGenome(vcf = Colombia_vcf, n = 500, threshold = 5)

#system('rm ColombiaGates_Pviv_Filtered.recode.vcf')

# Pv4----

path_to_Pv4 = '/seq/plasmodium-archive/jtennessen/Pvivax/malariagen'

files = c('Pv4_PvP01_01_v1.vcf.gz',
          'Pv4_PvP01_02_v1.vcf.gz',
          'Pv4_PvP01_03_v1.vcf.gz',
          'Pv4_PvP01_04_v1.vcf.gz',
          'Pv4_PvP01_05_v1.vcf.gz',
          'Pv4_PvP01_06_v1.vcf.gz',
          'Pv4_PvP01_07_v1.vcf.gz',
          'Pv4_PvP01_08_v1.vcf.gz',
          'Pv4_PvP01_09_v1.vcf.gz',
          'Pv4_PvP01_10_v1.vcf.gz',
          'Pv4_PvP01_11_v1.vcf.gz',
          'Pv4_PvP01_12_v1.vcf.gz',
          'Pv4_PvP01_13_v1.vcf.gz',
          'Pv4_PvP01_14_v1.vcf.gz',
          'Pv4_PvP01_MIT_v1.vcf.gz')

Pv4_vcf = NULL

for(file in files){
  
  run_vcftools(gzvcf = file.path(path_to_Pv4, file),
               bash_file = paste0(gsub('.vcf.gz','',file),'_run_vcf.sh'),
               out = paste0(gsub('\\.vcf\\.gz','',file),'_Filtered'),
               positions = 'selected_pos.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  Pv4_vcf = rbind(Pv4_vcf,
                  load_vcf(vcf = paste0(gsub('\\.vcf\\.gz','',file),'_Filtered.recode.vcf')))
  
#  system(paste0('rm ', file.path(wd, paste0(gsub('\\.vcf\\.gz','',file),'_Filtered.recode.vcf'))))
  
}

Pv4_rGenome =  vcf2rGenome(vcf = Pv4_vcf, n = 500, threshold = 5)

# Guyana----

path_to_guyana = '/gsap/garage-protistvector/GuyanaGatesPviv/JointCalls/Oct2022/Guyana_Pviv_Oct2022.JointCall.filtered.combined.snpeff.vcf'

run_vcftools(vcf = path_to_guyana,
             bash_file = 'run_vcf_Guyana.sh',
             out = 'Guyana_Pviv_Oct2022_Filtered',
             positions = 'selected_pos.bed',
             recode = TRUE,
             recode_INFO_all = TRUE)

Guyana_vcf = load_vcf(vcf = 'Guyana_Pviv_Oct2022_Filtered.recode.vcf')
Guyana_rGenome =  vcf2rGenome(vcf = Guyana_vcf, n = 500, threshold = 5)
# system('rm Guyana_Pviv_Oct2022_Filtered.recode.vcf')


save.image('genomics_Pvivax_LAC.RData')
