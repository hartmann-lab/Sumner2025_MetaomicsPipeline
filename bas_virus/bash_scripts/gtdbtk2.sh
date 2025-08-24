#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH -t 08:00:00
#SBATCH --mem=120GB
#SBATCH --job-name="gtdbtk"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/gtdbtk/errout/gtdbtk_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/gtdbtk/errout/gtdbtk_%A_%a.err


module purge all
module load mamba 
module load gtdbtk 
#source activate /projects/b1180/software/conda_envs/gtdbtk-2.1.1

echo "Starting gtdbtk job"

#data_dir=/projects/b1042/HartmannLab/stefanie/be_meta/metawrap_out/bin_refine
data_dir=/projects/b1042/HartmannLab/lung_phage/semibin2/1_quality_bins
denovo_out_dir=/projects/b1042/HartmannLab/lung_phage/gtdbtk/denovo_wf
classify_out_dir=/projects/b1042/HartmannLab/lung_phage/gtdbtk/classify_wf

#param_store=/projects/b1052/stefanie/vlp_enrichment/scripts/sample_list_bin_refine.txt
#param_a=$(cat $param_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')

#n=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

#mkdir ${output_dir}/${param_a}

#gtdbtk classify_wf --genome_dir ${data_dir} --out_dir ${classify_out_dir} --cpus 30 --force --extension fa.gz
gtdbtk de_novo_wf --genome_dir ${data_dir} --bacteria --outgroup_taxon p__Patescibacteria --out_dir ${denovo_out_dir} --cpus 30 --force --extension fa

#gtdbtk classify_wf --genome_dir ${data_dir} --out_dir ${output_dir} --cpus 20 --force --extension fa
#gtdbtk de_novo_wf --genome_dir ${data_dir} --archaea --outgroup_taxon p__Altiarchaeota --out_dir ${output_dir} --cpus 10 --force --extension fa

#gtdbtk classify_wf --genome_dir ${data_dir}/${param_a}/metawrap_70_10_bins --out_dir ${output_dir}/${param_a} --cpus 30 --force --extension fa

#gtdbtk de_novo_wf --genome_dir ${data_dir}/${param_a}/metawrap_70_10_bins --bacteria --outgroup_taxon p__Patescibacteria --out_dir ${output_dir}/${param_a} --cpus 30 --force --extension fa

#mkdir ${output_dir}/${param_a}

#cd ${data_dir}/${param_a}/metawrap_70_10_bins

#cnt=0

#for f in *.fa ;
#    do 
#    ((cnt++))
#    cp ${f} ${output_dir}/${param_a}.${f}
    
#    mkdir ${output_dir}/${param_a}/${f}
    
#      echo $cnt ;
#done

echo $param_a
echo $n
echo "Finishing gtdbtk job"
