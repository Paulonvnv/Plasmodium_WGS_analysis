#!/bin/bash
source /broad/software/scripts/useuse
use R-4.1
Rscript /gsap/garage-protistvector/PanSAmer/Coding_regions_filter.R \
  -wd /gsap/garage-protistvector/PanSAmer/PerVenPviv2023 \
  -fd /gsap/garage-protistvector/ColombiaData/Pviv/AnalysisTools \
  -gzvcf /seq/plasmodium/PanSAmer/JointCall/VENPER2019_2022Partial.JointCall.filtered.combined.snpeff.vcf \
  -o PerVen_Pviv_2023 \
  -gff genes.gff \
  -n 500 \
  -thres 5