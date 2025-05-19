#!/bin/bash

SAMPLE_LIST=$1 # List of sample names regarding the cna files
rm ./logs/*
sbatch --array=1-$(cat ${SAMPLE_LIST} | wc -l) ./cna_annotation.sbatch $SAMPLE_LIST
