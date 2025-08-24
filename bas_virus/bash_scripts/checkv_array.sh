#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH -t 06:00:00
#SBATCH --mem=5GB
#SBATCH --job-name="checkv_lung_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/checkv_out/errout/checkv_lp_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/checkv_out/errout/checkv_lp_array_%A_%a.err
#SBATCH --array=1-253%20

module purge all
module load mamba

source activate /projects/b1180/software/conda_envs/checkv

echo "Starting checkv job"

data_dir=/projects/b1042/HartmannLab/lung_phage/genomad
output_dir=/projects/b1042/HartmannLab/lung_phage/checkv_out

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
echo $param_a

mkdir ${output_dir}/${param_a}

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
echo $n

checkv end_to_end ${data_dir}/${param_a}/${param_a}.contigs_summary/${param_a}.contigs_virus.fna ${output_dir}/${param_a} \
-d /projects/b1180/software/conda_envs/checkv/checkv-db-v1.5 \
--remove_tmp -t 10 

echo "Finishing checkV lung virus job"