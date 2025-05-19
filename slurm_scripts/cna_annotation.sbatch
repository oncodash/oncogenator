#!/bin/bash

#SBATCH --job-name=cna_annotation
#SBATCH --output=./logs/%A_%a.txt
#SBATCH --error=./logs/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

echo $1
ODANNOTATOR_PATH=""
CNA_PATH=""
ASCAT_ESTIMATES_CSV=""
OUTPUT_PATH=""
cd $ODANNOTATOR_PATH
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)

echo $sample
source $ODANNOTATOR_PATH/venv/bin/activate

srun python $ODANNOTATOR_PATH/main.py --annotator local --copy_number_alterations $CNA_PATH/${sample}.csv --ascatestimates $ASCAT_ESTIMATES_CSV --output $OUTPUT_PATH/${sample}.csv --cores 1
