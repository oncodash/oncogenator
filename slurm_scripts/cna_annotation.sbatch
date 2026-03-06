#!/bin/bash

#SBATCH --job-name=cna_annotation
#SBATCH --output=./logs/%A_%a.txt
#SBATCH --error=./logs/%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
##SBATCH --partition=hautaniemi-epouta

echo $1
ODANNOTATOR_PATH="/mnt/storageBig8/home/aimaaral/dev/oncogenator"
cd $ODANNOTATOR_PATH
files=($(eval "ls ${search}"))
path=${files[$SLURM_ARRAY_TASK_ID-1]}
sample=$(basename $path .csv)

echo $sample
source $ODANNOTATOR_PATH/venv/bin/activate

#srun python $ODANNOTATOR_PATH/main.py --annotator local --copy_number_alterations /mnt/storageSmall8/work/aimaaral/oncodash/cnGenesFolder/${sample}.cnGenes.csv --ascatestimates /mnt/storageSmall8/work/aimaaral/oncodash/combinedAscatEstimates.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/${sample}.csv --cores 1
#python $ODANNOTATOR_PATH/src/local_annotator.py --copy_number_alterations /mnt/storageBig8/resources/processed_data/HERCULES/WGS_somatic_cnv/copyNumber_allMatched/v2.9/cnGenesFolder/${sample}.csv --ascatestimates /mnt/storageBig8/resources/processed_data/HERCULES/WGS_somatic_cnv/copyNumber_allMatched/v2.9/combinedAscatEstimates.csv --sample_info /mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/sample_info_extended.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/local/cna/${sample}.csv
python $ODANNOTATOR_PATH/src/external_annotator.py --copy_number_alterations $1/${sample}.csv --oncokbcna --output $ODANNOTATOR_PATH/output/ext/cna/${sample}.csv
