#!/bin/bash -l

SECONDS=0
#########################################################################
# Name: Eitan Vilker
# Date: $(date)
# Description: Eitan Shell Solution
#########################################################################
# Set SCC project to charge
#$ -P agedisease
# Request compute resources
#$ -pe omp 28
# Set array
#$ -t 1-2
#$ -l mem_per_core=18G
# Name job
#$ -N LargeNormElNet
#setting the error file:
#$ -m ea
#$ -o Large_generate_reports.out
#$ -e Large_Test_generate_reports_sh.err

cancerType=$(awk -F',' -v id=$SGE_TASK_ID 'NR==id+1 { gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2 }' /restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/topCancerTypesLarge.csv)
Rscript generate_single_reports.R $1 "$cancerType"
# Rscript generate_single_reports.R $1 $2
