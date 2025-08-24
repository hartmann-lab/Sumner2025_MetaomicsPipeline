#!/bin/sh
#SBATCH -A p31750
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 1:00:00
#SBATCH --mem=2gb
#SBATCH --job-name="bbmap_mapping_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/bbmap_out/errout/bbmap_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/bbmap_out/errout/bbmap_array_%A_%a.err
#SBATCH --array=1-253%20


module purge all

module load BBMap
module load samtools

echo "Starting bbmap read to virus assembly mapping job"


illumina_dir=/projects/b1042/HartmannLab/stefanie/bowtie_out
output_dir=/projects/b1042/HartmannLab/lung_phage/bbmap_out

param_store=/projects/b1052/stefanie/lung_phage/sample_sheet.tsv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

mkdir ${output_dir}/${param_a}

##must activate command from index file
cd /projects/b1042/HartmannLab/lung_phage/bbmap_out/index
/projects/b1042/HartmannLab/stefanie/bowtie_out/DNA_B01_01/DNA_B01_01.bowtie.r1.fastq.gz
#bbmap.sh -Xmx150g -t=5 in1=${illumina_dir}/${param_a}_clean_1.fastq in2=${illumina_dir}/${param_a}_clean_2.fastq ambiguous=best covstats=${output_dir}/${param_a}/${param_a}.covstats.txt covhist=${output_dir}/${param_a}/${param_a}.covhist.txt basecov=${output_dir}/${param_a}/${param_a}.basecov.txt bincov=${output_dir}/${param_a}/${param_a}.bincov.txt 
#bbmap.sh -Xmx150g -t=5 in1=${illumina_dir}/${param_a}/${param_a}.bowtie.r1.fastq.gz in2=${illumina_dir}/${param_a}/${param_a}.bowtie.r2.fastq.gz ambiguous=best outm=${output_dir}/${param_a}/${param_a}.virus.sam covstats=${output_dir}/${param_a}/${param_a}.covstats.txt covhist=${output_dir}/${param_a}/${param_a}.covhist.txt basecov=${output_dir}/${param_a}/${param_a}.basecov.txt bincov=${output_dir}/${param_a}/${param_a}.bincov.txt rpkm=${output_dir}/${param_a}/${param_a}.rpkm.txt 

samtools view -bS -F4 ${output_dir}/${param_a}/${param_a}.virus.sam > ${output_dir}/${param_a}/${param_a}.bam


#the usual wildcards don't work to remove cat with bbmap, tootbrushes and showerheads have different suffixes so when the
#first line of code fails when the array reaches the toothbrush files, it should work with the second line
#bbmap.sh in1=${illumina_dir}/${param_a}_1.fastq in2=${illumina_dir}/${param_a}_2.fastq ambiguous=toss covstats=${output_dir}/${param_a}/${param_a}.covstats.txt covhist=${output_dir}/${param_a}/${param_a}.covhist.txt basecov=${output_dir}/${param_a}/${param_a}.basecov.txt bincov=${output_dir}/${param_a}/${param_a}.bincov.txt 


echo $param_a
echo $n
echo "Finishing bbmap read to virus assembly mapping job"

