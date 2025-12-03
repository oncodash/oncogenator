#!/bin/bash

#SBATCH --job-name=cna_annotation
#SBATCH --output=./logs/%A_%a.txt
#SBATCH --error=./logs/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

echo $1
ODANNOTATOR_PATH="/mnt/storageBig8/home/aimaaral/dev/oncogenator/"
cd $ODANNOTATOR_PATH
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1)

echo $sample
source $ODANNOTATOR_PATH/venv/bin/activate

#srun python $ODANNOTATOR_PATH/main.py --annotator local --copy_number_alterations /mnt/storageSmall8/work/aimaaral/oncodash/cnGenesFolder/${sample}.cnGenes.csv --ascatestimates /mnt/storageSmall8/work/aimaaral/oncodash/combinedAscatEstimates.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/${sample}.csv --cores 1
srun python $ODANNOTATOR_PATH/local_annotator.py --copy_number_alterations /mnt/storageSmall8/work/aimaaral/oncodash/cnGenesFolder/${sample}.cnGenes.csv --ascatestimates /mnt/storageSmall8/work/aimaaral/oncodash/combinedAscatEstimates.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/local/cna/${sample}.csv --cores 1
#srun python $ODANNOTATOR_PATH/external_annotator.py --copy_number_alterations /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/${sample}.csv --oncokbcna --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output_external/${sample}.csv
