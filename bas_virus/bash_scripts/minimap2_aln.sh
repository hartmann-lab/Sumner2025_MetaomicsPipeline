#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH -t 5:00:00
#SBATCH --mem=10GB
#SBATCH --job-name="minimap2_aln_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/minimap_out/errout/minimap2_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/minimap_out/errout/minimap2_array_%A_%a.err
#SBATCH --array=1-254%20

module purge all
#module load mamba
module load minimap2
module load samtools
#source activate /projects/b1052/stefanie/software/coverm

echo "Starting minimap2 alignment"

assembly_dir=/projects/b1042/HartmannLab/lung_phage/assemblies/megahit_paired2

kneaddata_dir=/projects/b1042/HartmannLab/lung_phage/reads

output_dir=/projects/b1042/HartmannLab/lung_phage/minimap_out/assembly_aln

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

#mkdir ${output_dir}/${param_a}
#mkdir ${output_dir}/reactor_coassembly
#mkdir ${output_dir}/${param_b}


cd /projects/b1042/HartmannLab/lung_phage/reads
cat ${param_a}_paired_1.fastq ${param_a}_paired_2.fastq \
${param_a}_unmatched_1.fastq ${param_a}_unmatched_2.fastq \
${param_a}.repeats.removed.1.fastq ${param_a}.repeats.removed.2.fastq \
${param_a}.repeats.removed.unmatched.1.fastq ${param_a}.repeats.removed.unmatched.2.fastq \
${param_a}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq \
${param_a}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq \
${param_a}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq \
${param_a}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq > \
${param_a}_assembly_cat.fastq

# assembly alignment
minimap2 -x sr -a -t 10 ${assembly_dir}/${param_a}/${param_a}.contigs.fa ${param_a}_assembly_cat.fastq | \
samtools view -bS -F4 - | samtools sort -\@ 16 - -o ${output_dir}/${param_a}_aln.bam

#coverm contig -b ${output_dir}/${param_b}/${param_b}_aln.bam \
#-m mean rpkm tpm covered_fraction \
#-o /projects/b1052/stefanie/vlp_enrichment/coverm/assembly_alignment_stats/${param_b}_assembly_LR_coverm.tsv

#SR alignment
#minimap2 -x sr -a -t 10 ${SR_fasta} ${kneaddata_dir}/${param_b}/${param_b}.fastq | \
#samtools view -bS -F4 - | samtools sort -\@ 16 - -o ${output_dir}/${param_b}_SR_aln.bam

#LR alignment
#minimap2 -x map-ont -a -t 10 ${LR_fasta} ${ont_dir}/${param_a}_LR_filtered.fastq | \
#samtools view -bS -F4 - | samtools sort -\@ 16 - -o ${output_dir}/${param_b}_LR_aln.bam


echo $param_a
echo $n
echo "Finishing mapping"


#coverm contig -b /projects/b1052/stefanie/vlp_enrichment/minimap2/LR_viral_alignments/LR_reactor_alignments/NOBR_01_LR_aln.bam -m mean rpkm tpm covered_fraction