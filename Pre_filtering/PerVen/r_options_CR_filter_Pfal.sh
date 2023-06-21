#!/bin/bash
source /broad/software/scripts/useuse
use R-4.1
Rscript /gsap/garage-protistvector/PanSAmer/Coding_regions_filter.R \
  -wd /gsap/garage-protistvector/PanSAmer/PerVenPfal2023 \
  -fd /gsap/garage-protistvector/ColombiaData/Pviv/AnalysisTools \
  -gzvcf /seq/plasmodium/PanSAmer/JointCall/Peru.Ven.Pfal.2023.snp.indel.recalibrated.filtered.ALL.snpeff.vcf \
  -o PerVen_Pfal_2023 \
  -gff PlasmoDB-59_Pfalciparum3D7.gff \
  -n 500 \
  -thres 5