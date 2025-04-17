#!/bin/bash -l

SECONDS=0
#########################################################################
# Name: Eitan Vilker
# Date: $(date)
# Description: Eitan Shell Solution
# Parameters: First parameter is name of test

module load miniconda
mamba activate MontiLab # Change to your environment or comment out if not using

qsub generate_reports.sh $1
# qsub generate_large_reports.sh $1 "BRCA"
# qsub generate_medium_reports.sh $1 "LGG"
