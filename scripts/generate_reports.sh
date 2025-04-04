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
#$ -pe omp 3
# Specify job limits
#$ -l h_rt=12:00:00
#$ -l mem_per_core=8G 
# Name job
#$ -N Test_Ridge005_generate_reports_sh 
# Send an email on job completion or failure
#$ -m ea
#$ -M evilker@bu.edu 
#setting the error file:
#$ -o generate_reports.out
#$ -e Test_generate_reports_sh.err

Rscript generate_reports.R "Ridge005"
