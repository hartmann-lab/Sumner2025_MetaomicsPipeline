#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 4:00:00
#SBATCH --mem=120GB
#SBATCH --job-name="iphop_add_to_db"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1052/stefanie/vlp_enrichment/iphop/errout/iphop_VLP_atd_%A_%a.out
#SBATCH --error=/projects/b1052/stefanie/vlp_enrichment/iphop/errout/iphop_VLP_atd_%A_%a.err



module purge all
module load mamba 

source activate /projects/b1052/stefanie/software/iphop

echo "Starting iphop job"

#iphop add_to_db --fna_dir /path/to/dir/containing/all/MAGs/to/be/added/metawrap_70_10  --gtdb_dir /path/to/gtdbtk/de_novo_wf --out_dir /path/to/new/db/in/iphop/db/folder --db_dir /path/to/original/db

#must have fna extension
#awk '{print "cp "$2 ".fa " $2 ".fna"}' /projects/b1042/HartmannLab/lung_phage/Rproject/bas_virus/clean_data/semibin2_qual_bins.tsv |bash

iphop add_to_db --fna_dir /projects/b1042/HartmannLab/lung_phage/iphop/atd_bins --gtdb_dir /projects/b1042/HartmannLab/lung_phage/gtdbtk/denovo_wf --out_dir /projects/b1052/stefanie/software/iphop_db/Aug_2023_atdlung0525_pub_rw --db_dir /projects/b1052/stefanie/software/iphop_db/Aug_2023_pub_rw

echo $param_a
echo $n
echo "Finishing iphop job"
