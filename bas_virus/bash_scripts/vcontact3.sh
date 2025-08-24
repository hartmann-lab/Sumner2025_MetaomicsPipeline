#!/bin/sh
#SBATCH -A p31750
#SBATCH -p normal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 02:00:00
#SBATCH --mem=40GB
#SBATCH --job-name="vcontact3"
#SBATCH --mail-user=stefaniehuttelmaier2024@u.northwestern.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/projects/b1042/HartmannLab/lung_phage/vcontact3/errout/vcontact3_%A_%a.out
#SBATCH --error=/projects/b1042/HartmannLab/lung_phage/vcontact3/errout/vcontact3_%A_%a.err

module purge all
module load mamba
source activate /projects/b1052/stefanie/software/vcontact3

vcontact3 run --nucleotide /projects/b1042/HartmannLab/lung_phage/derep_viruses/lung_vOTUs_filtered.fna \
--output /projects/b1042/HartmannLab/lung_phage/vcontact3 --db-domain "prokaryotes" --db-version 220 \
--db-path /projects/b1052/stefanie/db -e cytoscape

#vcontact3 run --nucleotide /projects/b1052/stefanie/vlp_enrichment/vOTUs/vlpe_vOTUs.fna \
#--output /projects/b1052/stefanie/vlp_enrichment/vcontact3/vOTUs --db-domain "prokaryotes" --db-version 220 \
#--db-path /projects/b1052/stefanie/db -e newick

#vcontact3 run --nucleotide /projects/b1052/stefanie/vlp_enrichment/checkv/genomad/VLP_2_1_3/viruses.fna \
#--output /projects/b1052/stefanie/vlp_enrichment/vcontact3/VLP_2 --db-domain "prokaryotes" --db-version 220 \
#--db-path /projects/b1052/stefanie/db -e newick


#vcontact3 run --nucleotide /projects/b1052/stefanie/vlp_enrichment/checkv/genomad/VLP_3_1_2/viruses.fna \
#--output /projects/b1052/stefanie/vlp_enrichment/vcontact3/VLP_3 --db-domain "prokaryotes" --db-version 220 \
#--db-path /projects/b1052/stefanie/db -e newick

#vcontact3 run --nucleotide /projects/b1052/stefanie/vlp_enrichment/vOTUs/vlpe_vOTUs.fna \
#--output /projects/b1052/stefanie/vlp_enrichment/vcontact3/vOTUs --db-domain "prokaryotes" --db-version 220 \
#--db-path /projects/b1052/stefanie/db -e cytoscape