#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -t 3:00:00
#SBATCH --mem=20GB
#SBATCH --job-name="semibin_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/semibin2/errout/semibin_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/semibin2/errout/semibin_%A_%a.err
#SBATCH --array=1-252%20

module purge all
module load samtools
module load mamba
source activate /projects/b1052/stefanie/software/semibin


assembly_fasta=/projects/b1042/HartmannLab/lung_phage/assemblies/megahit_paired2
bam_dir=/projects/b1042/HartmannLab/lung_phage/minimap_out/assembly_aln
out_dir=/projects/b1042/HartmannLab/lung_phage/semibin2

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
param_b=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

mkdir ${out_dir}/${param_a}
#mkdir ${output_dir}/reactor_coassembly

#samtools sort bam
samtools sort -o ${bam_dir}/${param_a}_aln.sorted.bam ${bam_dir}/${param_a}_aln.bam


# run binner
SemiBin2 single_easy_bin -i ${assembly_fasta}/${param_a}/${param_a}.contigs.fa \
-b ${bam_dir}/${param_a}_aln.sorted.bam --threads 24 \
-o ${out_dir}/${param_a}