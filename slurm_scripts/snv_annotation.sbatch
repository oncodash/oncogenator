#!/bin/bash

#SBATCH --job-name=csomatic_mutation_annotation
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

#srun python $ODANNOTATOR_PATH/main.py --annotator local --somatic_variants /mnt/storageSmall8/work/aimaaral/oncodash/byPatient_CSVs/ --cn_annotations /mnt/storageSmall8/work/aimaaral/oncodash/cnGenesFolder/ --ascatestimates /mnt/storageSmall8/work/aimaaral/oncodash/combinedAscatEstimates.csv --output /mnt/storageSmall8/work/aimaaral/oncodash/somatic_mutation_annotations/${pid}.csv --pid ${pid} --cores 1
srun python $ODANNOTATOR_PATH/local_annotator.py --somatic_variants /mnt/storageBig8/resources/processed_data/HERCULES/WGS_somatic_variants/variants_v4.9/byPatient_CSVs/ --cn_annotations /mnt/storageSmall8/work/aimaaral/oncodash/cnGenesFolder --ascatestimates /mnt/storageSmall8/work/aimaaral/oncodash/combinedAscatEstimates.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/local/somatic_mutation/${sample}.csv --pid ${sample}
#srun python $ODANNOTATOR_PATH/external_annotator.py --oncokbsomatic_mutation --somatic_variants /mnt/storageBig8/home/aimaaral/dev/oncogenator/output/${sample}.csv --output /mnt/storageBig8/home/aimaaral/dev/oncogenator/output_ext/${sample}.csv

