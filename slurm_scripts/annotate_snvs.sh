#!/bin/bash

rm ./logs/*
sbatch --array=1-$(cat ${1} | wc -l) ./somatic_mutation_annotation.sbatch $1
