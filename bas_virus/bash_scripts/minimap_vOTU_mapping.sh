#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 08:00:00
#SBATCH --mem=5gb
#SBATCH --job-name="minimap"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/minimap_out/errout/minimap_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/minimap_out/errout/minimap_array_%A_%a.err
#SBATCH --array=1-254%20


module purge all

module load minimap2
module load samtools

echo "Starting minimap read to high qual vOTU mapping job"


ref_fasta=/projects/b1042/HartmannLab/lung_phage/derep_viruses/vOTUs_genomad.fna
illumina_dir=/projects/b1042/HartmannLab/lung_phage/reads
output_dir=/projects/b1042/HartmannLab/lung_phage/minimap_out/vOTU_aln

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)


minimap2 -x sr -a -t 10 ${ref_fasta} ${illumina_dir}/${param_a}_assembly_cat.fastq | \
samtools view -bS -F4 - | samtools sort -\@ 16 - -o ${output_dir}/${param_a}.bam


echo "alignment complete"

echo $param_a
echo $n
echo "Finishing minimap mapping job"

