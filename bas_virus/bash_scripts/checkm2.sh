#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=5
#SBATCH -t 3:00:00
#SBATCH --mem=20GB
#SBATCH --job-name="checkm2_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/checkm2/errout/checkm_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/checkm2/errout/checkm_%A_%a.err
#SBATCH --array=1-253%20

module purge all
module load mamba
source activate /projects/b1052/stefanie/software/checkm2


bin_dir=/projects/b1042/HartmannLab/lung_phage/semibin2
out_dir=/projects/b1042/HartmannLab/lung_phage/checkm2

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
param_b=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

mkdir ${out_dir}/${param_a}

# run checkm2
checkm2 predict --threads 5 --input ${bin_dir}/${param_a}/output_bins --output-directory ${out_dir}/${param_a} -x .fa.gz
