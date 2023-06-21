#!/bin/bash
source /broad/software/scripts/useuse
use R-4.1
Rscript /gsap/garage-protistvector/PanSAmer/Split_chromosomes.R \
  -wd /gsap/garage-protistvector/PanSAmer/PerVenPfal2023 \
  -fd /gsap/garage-protistvector/ColombiaData/Pviv/AnalysisTools \
  -vcf /gsap/garage-protistvector/PanSAmer/PerVenPfal2023/PerVen_Pfal_2023_CodingRegionsOnly.recode.vcf \
  -chr $SGE_TASK_ID \
  -o PerVen_Pfal_2023_CodingRegionsOnly \
  -n 500 \
  -thres 5 \
  -p Pf