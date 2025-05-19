#!/bin/bash

#SBATCH --job-name=snv_annotation
#SBATCH --output=./logs/%A_%a.txt
#SBATCH --error=./logs/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

echo $1
ODANNOTATOR_PATH=""
SNV_PATH=""
CNA_PATH=""
ASCAT_ESTIMATES_CSV=""
OUTPUT_PATH=""
cd $ODANNOTATOR_PATH
pid=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)

echo $pid
source $ODANNOTATOR_PATH/venv/bin/activate

srun python $ODANNOTATOR_PATH/local_annotator.py --somatic_variants $SNV_PATH --cn_annotations $CNA_PATH --ascatestimates $ASCAT_ESTIMATES_CSV --output $OUTPUT_PATH/${pid}.csv --pid ${pid} --cores 1
srun python $ODANNOTATOR_PATH/external_annotator.py --somatic_variants $SNV_PATH --cn_annotations $CNA_PATH --ascatestimates $ASCAT_ESTIMATES_CSV --output $OUTPUT_PATH/${pid}.csv --pid ${pid} --cores 1
