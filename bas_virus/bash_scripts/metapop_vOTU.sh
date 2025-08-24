#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH -t 12:00:00
#SBATCH --mem=20GB
#SBATCH --job-name="metapop"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1052/stefanie/lung_phage/metapop_out/errout/metapop1.out
#SBATCH --error=/projects/b1052/stefanie/lung_phage/metapop_out/errout/metapop1.err


module purge all
module load python 
module load bcftools 

source activate /projects/b1180/software/conda_envs/metapop

echo "Starting metapop job"

#data_dir=/projects/b1052/Stefanie/be_meta/vibrant
#output_dir=/projects/b1042/HartmannLab/stefanie/be_meta/iphop_high_qual_out/iphop_cat_HQ_out
#db_dir=/projects/b1180/software/conda_envs/iphop/iphop/iphop_db/Sept_2021_pub_rw

#param_store=/projects/b1052/Stefanie/be_meta/scripts/sample_sheet.tsv
#param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')


metapop --input_samples /projects/b1052/stefanie/lung_phage/bbmap_out/viral_bam --reference /projects/b1052/stefanie/lung_phage/metapop_out/fna --norm /projects/b1052/stefanie/lung_phage/scripts/normalization_table.tsv --out /projects/b1052/stefanie/lung_phage/metapop_out --skip_preproc --no_micro


echo $n

##metapop --input_samples /projects/b1042/HartmannLab/stefanie/be_meta/bwa_cdhit_HQALL_out/final_bam --reference /projects/b1042/HartmannLab/stefanie/be_meta/high_qual_unbinned/metapop_fa --out /projects/b1042/HartmannLab/stefanie/be_meta/metapop/vOTU_high_qual --macro_only
