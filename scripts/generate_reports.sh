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
#$ -t 1-15
# Specify job limits
# Commented out -l h_rt=12:00:00
#$ -l mem_per_core=16G
# Name job
#$ -N ComboRidge005
#setting the error file:
#$ -o generate_reports.out
#$ -e Test_generate_reports_sh.err

taskinput=$(awk -F',' -v id=$SGE_TASK_ID 'NR==id+1 { gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2 }' /restricted/projectnb/agedisease/projects/pancancer_aging_clocks/results/topCancerTypes.csv)
Rscript generate_single_reports.R $1 "$taskinput"
