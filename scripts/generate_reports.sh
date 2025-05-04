#!/bin/bash -l

SECONDS=0
#########################################################################
# Name: Eitan Vilker
# Date: $(date)
# Description: Eitan Shell Solution
# Parameters: Test name
#########################################################################
# Set SCC project to charge
#$ -P agedisease
# Request compute resources
#$ -pe omp 16
# Set array
#$ -t 1-14
# Commented out -l mem_per_core=16G
# Name job
#$ -N NewCovariatesNormalElNet
# Setting the error file:
#$ -o small_generate_reports.out
#$ -e small_generate_reports_sh.err

cancerType=$(awk -F',' -v id=$SGE_TASK_ID 'NR==id+1 { gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2 }' /restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/topCancerTypes.csv)
Rscript generate_single_reports.R $1 "$cancerType"
