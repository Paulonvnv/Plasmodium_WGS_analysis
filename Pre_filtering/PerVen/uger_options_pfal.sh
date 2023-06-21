#!/bin/bash
qsub -l h_vmem=32G \
   -l h_rt=05:00:00 \
   -o /gsap/garage-protistvector/PanSAmer/output2/ \
   /gsap/garage-protistvector/PanSAmer/r_options_CR_filter_Pfal.sh