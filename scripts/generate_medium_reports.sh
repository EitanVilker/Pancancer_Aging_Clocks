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
#$ -pe omp 16
#$ -l mem_per_core=16G
# Name job
#$ -N MediumAgeClockReport
#setting the error file:
#$ -o Medium_generate_reports.out
#$ -e Medium_Test_generate_reports_sh.err

# cancerType=$(awk -F',' -v id=$SGE_TASK_ID 'NR==id+1 { gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2 }' /restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/topCancerTypesMedium.csv)
# Rscript generate_single_reports.R $1 "$cancerType"
Rscript generate_single_reports.R $1 $2
