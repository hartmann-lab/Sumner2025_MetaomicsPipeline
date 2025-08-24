#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=24
#SBATCH -t 12:00:00
#SBATCH --mem=50GB
#SBATCH --job-name="comebin_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1052/stefanie/vlp_enrichment/comebin/errout/comebin_%A_%a.out
#SBATCH --error=/projects/b1052/stefanie/vlp_enrichment/comebin/errout/comebin_%A_%a.err
#SBATCH --array=1-19%1

module purge all
module load mamba
source activate /projects/b1052/stefanie/software/comebin
cd /projects/b1052/stefanie/software/comebin/bin/COMEBin/scripts


assembly_fasta=/projects/b1052/stefanie/vlp_enrichment/np_metaMBG

kneaddata_dir=/projects/b1042/HartmannLab/stefanie/mlm2/results/MGX/kneaddata
ont_dir=/projects/b1052/stefanie/vlp_enrichment/nanopore_qc/chopper

out_cov_dir=/projects/b1052/stefanie/vlp_enrichment/comebin/coverage_files/
out_bam_dir=/projects/b1052/stefanie/vlp_enrichment/comebin/bam_files
out_bin_dir=/projects/b1052/stefanie/vlp_enrichment/comebin/bin_files

param_store=/projects/b1052/stefanie/vlp_enrichment/bash_scripts/polish_sample_list.txt
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
param_b=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

mkdir ${out_cov_dir}/${param_b}
mkdir ${out_bam_dir}/${param_b}
mkdir ${out_bin_dir}/${param_b}
#mkdir ${output_dir}/reactor_coassembly


#LR assembly
# align reads to assembly
bash gen_cov_file.sh -a ${assembly_fasta}/${param_a}/contigs.fasta.gz \
-o ${out_cov_dir}/${param_b} \
-b ${out_bam_dir}/${param_b} \
--single-end ${ont_dir}/${param_a}_LR_filtered.fastq

# run binner
bash run_comebin.sh -a ${assembly_fasta}/${param_a}/contigs.fasta.gz \
-o ${out_bin_dir}/${param_b} \
-p ${out_bam_dir}/${param_b} \
-t 24