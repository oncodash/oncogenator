#!/bin/bash

SAMPLE_LIST=$1 # List of sample names regarding the snv files
rm ./logs/*
sbatch --array=1-$(cat ${SAMPLE_LIST} | wc -l) ./snv_annotation.sbatch $SAMPLE_LIST
