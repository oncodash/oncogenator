#!/bin/bash

rm ./logs/*
sbatch --export=search=${1} --array=1-$(ls ${1} | wc -l) ./somatic_mutation_annotation.sbatch $1
