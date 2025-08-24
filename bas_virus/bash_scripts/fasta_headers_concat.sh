#!/bin/sh
#SBATCH -A p31750
#SBATCH -p short
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 00:20:00
#SBATCH --mem=1GB
#SBATCH --job-name="add samples to header"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/genomad_out/errout/update_header_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/genomad_out/errout/update_header_array_%A_%a.err
#SBATCH --array=1-252%50

## script serves 2 purposes
## 1. add sample identifier to contig headers for parsing downstream
## 2. concatenate modified fasta files

module purge all

echo "Starting update header job"

## set data directories 

genomad_dir=/projects/b1042/HartmannLab/lung_phage/genomad
output_dir=/projects/b1042/HartmannLab/lung_phage/genomad/1_virus_fasta

## give list of sample identifiers to array executor
param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

## 1. add identifier to fasta file headers

cd ${genomad_dir}/${param_a}/${param_a}.contigs_summary/
sed 's/>/>'${param_a}'_/I' ${param_a}.contigs_virus.fna > ${output_dir}/${param_a}_genomad_virus.fna

# plasmid fastas
#cd ${genomad_dir}/${param_a}/${param_a}.contigs_summary/
#sed 's/>/>'${param_a}'_/I' ${param_a}.contigs_plasmid.fna > ${output_dir}/${param_a}_genomad_plasmid.fna


## 2. concatenate vibrant, virsorter2 and genomad viral contigs into single fasta file 

#cat ${vibrant_dir}/${param_a}/VIBRANT_${param_a}.scaffolds/VIBRANT_phages_${param_a}.scaffolds/${param_a}-vibrant-combined-header.fna ${virsorter_dir}/${param_a}/${param_a}-virsorter-combined-header.fa ${genomad_dir}/${param_a}/${param_a}.scaffolds_summary/${param_a}-genomad-combined-header.fna > /projects/b1052/Stefanie/be_meta/virus_cat/${param_a}_virus_cat3.fasta

#echo "cat complete" 


echo $param_a
echo $n
echo "Finishing cdhit dereplication job"

