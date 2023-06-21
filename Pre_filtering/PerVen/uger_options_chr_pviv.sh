#!/bin/bash
qsub -l h_vmem=32G \
   -l h_rt=10:00:00 \
   -o /gsap/garage-protistvector/PanSAmer/output2/ \
   -t 1-16 \
   /gsap/garage-protistvector/PanSAmer/r_options_chr_filter_Pviv.sh