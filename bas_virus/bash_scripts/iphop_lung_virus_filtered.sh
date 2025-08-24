#!/bin/sh
#SBATCH -A p31750
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --ntasks=5
#SBATCH -t 06:00:00
#SBATCH --mem=200GB
#SBATCH --job-name="iphop_high_qual"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/iphop/errout/iphop_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/iphop/errout/iphop_%A.err

module purge all
module load mamba

source activate /projects/b1052/stefanie/software/iphop


echo "Starting iphop job"

fasta_file=/projects/b1042/HartmannLab/lung_phage/derep_viruses/vOTUs_genomad.fna
output_dir=/projects/b1042/HartmannLab/lung_phage/iphop
db_dir=/projects/b1052/stefanie/software/iphop_db/Aug_2023_atdlung0525_pub_rw


iphop predict --fa_file ${fasta_file} --db_dir ${db_dir} --out_dir ${output_dir} -t 5


echo $param_a
echo $n
echo "Finishing iphop job"
