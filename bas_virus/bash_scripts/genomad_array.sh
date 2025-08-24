#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=5
#SBATCH -t 12:00:00
#SBATCH --mem=36GB
#SBATCH --job-name="genomad_lp_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/genomad/errout/genomad_lp_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/genomad/errout/genomad_lp_array_%A_%a.err
#SBATCH --array=1-253%20


module purge all
module load mamba

source activate /projects/b1052/stefanie/software/genomad

echo "Starting genomad job"

data_dir=/projects/b1042/HartmannLab/lung_phage/assemblies/megahit_paired2
output_dir=/projects/b1042/HartmannLab/lung_phage/genomad

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
echo $param_a

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
echo $n

mkdir ${output_dir}/${param_a}

genomad end-to-end --cleanup ${data_dir}/${param_a}/${param_a}.contigs.fa ${output_dir}/${param_a} /projects/b1052/stefanie/software/genomad_db


echo "Finishing genomad job"