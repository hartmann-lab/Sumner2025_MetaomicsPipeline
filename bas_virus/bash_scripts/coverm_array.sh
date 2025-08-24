#!/bin/sh
#SBATCH -A p31750
#SBATCH -p short
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 00:15:00
#SBATCH --mem=5GB
#SBATCH --job-name="coverm_aln_array"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/coverm/errout/coverm_array_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/coverm/errout/coverm_array_%A_%a.err
#SBATCH --array=1-254%20

module purge all
module load mamba
source activate /projects/b1052/stefanie/software/coverm

echo "Starting coverm contig stats"

vOTU_fasta=/projects/b1042/HartmannLab/lung_phage/derep_viruses/vOTUs_genomad.fna

reads_dir=/projects/b1042/HartmannLab/lung_phage/reads

bam_dir=/projects/b1042/HartmannLab/lung_phage/minimap_out/vOTU_aln

output_dir=/projects/b1042/HartmannLab/lung_phage/coverm

param_store=/projects/b1042/HartmannLab/lung_phage/scripts/samples_MGX.csv
param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
param_b=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

#mkdir ${output_dir}/${param_a}
#mkdir ${output_dir}/reactor_coassembly

#SR alignment
#coverm contig -b ${SR_bam_dir}/${param_b}_SR_aln.bam \
#-m mean rpkm tpm covered_fraction \
#--min-read-percent-identity 0.9 --min-covered-fraction 0.4 \
#-o ${output_dir}/${param_b}_SR_coverm.tsv

#LR alignment
coverm contig -b ${bam_dir}/${param_a}.bam \
-m mean rpkm tpm covered_fraction \
--min-read-percent-identity 0.9 --min-covered-fraction 0.4 \
-o ${output_dir}/${param_a}_coverm.tsv

echo $param_a
echo $n
echo "Finishing stat calculations"