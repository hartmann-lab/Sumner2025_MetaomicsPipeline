# Import viral data
# Purpose: Import and clean raw genomad and checkv files

# Load libraries ----
library(tidyverse)
library(here)
library(skimr)

# Import genomad files ---
# Import genomad for taxonomy
get_genomad <- function(sample){
  genomad_summary <- read_tsv(sample, col_names = TRUE) |> 
    rename(feature = seq_name ) %>% 
    mutate(sample = sample) %>% 
    separate(col=sample, into=c(NA,NA,NA,NA,NA,NA,NA,NA, NA,"sample"), sep="/") %>%
    mutate(sample = str_sub(sample, end = -27)) %>% 
    select(c(feature, everything()))
  return(genomad_summary)
}

genomad_file_list <- list.files("/projects/b1042/HartmannLab/lung_phage/genomad/", 
                                pattern="*_virus_summary.tsv", 
                                recursive = TRUE,
                                full.names=TRUE)

genomad_table_list <- lapply(genomad_file_list, get_genomad)

genomad_fulltable <- bind_rows(lapply(genomad_table_list, function(x) if(nrow(x) == 0) NULL else x))

genomad_fulltable <- genomad_fulltable |> 
  mutate('contig_id' = paste0(sample, "_", feature)) |> 
  select(c(feature, contig_id, everything()))

count(genomad_fulltable, length>=2500)

skimr::skim(genomad_fulltable)

# filter genomad for high quality virus ----

genomad_filter <- genomad_fulltable %>% 
  filter(length >= 2500, n_hallmarks >= 1, virus_score >= 0.8 )

min(genomad_filter$virus_score)

skimr::skim(genomad_filter)

write.table(genomad_filter, file='clean_data/genomad_quality_filtered.tsv', quote=FALSE,row.names = F, sep='\t')

# use quality filtered table to parse viral contigs for dereplication to vOTUs
# use awk on command line

# genomad taxa table ---

vOTU_list <- read.csv(file = "/projects/b1042/HartmannLab/lung_phage/derep_viruses/vOTUs.txt", col.names = c('feature'), header = FALSE ) %>% 
  as_tibble()


vOTU_taxa_table_all <- genomad_fulltable %>% 
  select(contig_id, taxonomy, sample) %>% 
  rename(feature = contig_id) %>% 
  separate(col = taxonomy, into = c("classified" ,"realm", "kingdom", "phylum", "class", "order", "family"), sep = ';')

vOTU_taxa_table <- vOTU_list %>% 
  left_join(vOTU_taxa_table_all) %>% 
  select(-c(classified))

skim(vOTU_taxa_table)

write.table(vOTU_taxa_table, file='clean_data/vOTU_taxa_table.tsv', quote=FALSE,row.names = F, sep='\t')

# Import checkv files ----

get_checkv <- function(sample){
  checkv_summary <- read_tsv(sample, col_names = TRUE) |> 
    rename(feature = contig_id) %>% 
    mutate(sample = sample) %>% 
    separate(col=sample, into=c(NA,NA,NA,NA,NA,NA,"sample", NA), sep="/") %>%
    select(c(feature, everything()))
  return(checkv_summary)
}

checkv_file_list <- list.files("/projects/b1042/HartmannLab/lung_phage/checkv_out", 
                                pattern="quality_summary.tsv", 
                                recursive = TRUE,
                                full.names=TRUE)

checkv_table_list <- lapply(checkv_file_list, get_checkv)

checkv_fulltable <- bind_rows(checkv_table_list)

checkv_fulltable <- checkv_fulltable |> 
  mutate('contig_id' = paste0(sample, "_", feature)) |> 
  select(c(feature, contig_id, everything()))

skimr::skim(checkv_fulltable)

count(checkv_fulltable, contig_length>=2000)
count(checkv_fulltable, checkv_quality=='Complete')
min(checkv_fulltable$contig_length)
unique(checkv_fulltable$sample)
max(checkv_fulltable$contig_length)

virus_quality_filter <- checkv_fulltable %>% 
  filter(checkv_quality=='Medium-quality'| checkv_quality=='High-quality'| checkv_quality=='Complete') %>% 
  filter(contig_length>=2000)

unique(virus_quality_filter$sample)

# read in coverm abundance tables ----
get_coverm <- function(sample){
  coverm_summary <- read_tsv(sample, skip = 1,
                             col_names = c("feature", "mean", "rpkm", "tpm", "covered_fraction")) |> 
    mutate(sample = sample) %>% 
    separate(col=sample, into=c(NA,NA,NA,NA,NA, NA, "sample"), sep="/") %>%
    mutate(sample = str_sub(sample, end = -12)) %>% 
    select(c(feature, everything()))
  return(coverm_summary)
}

coverm_file_list <- list.files("/projects/b1042/HartmannLab/lung_phage/coverm", 
                               pattern="*_coverm.tsv", 
                               recursive = TRUE,
                               full.names=TRUE)

coverm_table_list <- lapply(coverm_file_list, get_coverm)

coverm_fulltable <- bind_rows(coverm_table_list)

skimr::skim(coverm_fulltable)

write.table(coverm_fulltable, file='clean_data/vOTU_coverm_depth.tsv', quote=FALSE,row.names = F, sep='\t')

## rpkm abundance table ----

votu_feature_table <- coverm_fulltable %>% 
  select(feature, sample, rpkm) %>% 
  pivot_wider(names_from = 'sample', values_from = 'rpkm' )

write.table(votu_feature_table, file='clean_data/vOTU_feature_table.tsv', quote=FALSE,row.names = F, sep='\t')

# read in iphop tables ----
# separate host genus into new columns to sort by 
#'h_' flag is meant to distinguish host from phage taxa downstream
iphop_genus <- read_csv(here("/projects/b1042/HartmannLab/lung_phage/iphop/derep_vOTUs/Host_prediction_to_genus_m90.csv")) %>% 
  rename(feature = 'Virus', 
         host_genus = 'Host genus',
         aai_rafah = 'AAI to closest RaFAH reference',
         confidence_score = 'Confidence score',
         methods = 'List of methods') |> 
  separate(col = host_genus, into = c("h_domain" ,"h_phylum", "h_class", "h_order", "h_family", "h_genus"), sep = ';') %>% 
  mutate(h_domain = str_sub(h_domain, start = 4),
         h_phylum = str_sub(h_phylum, start = 4),
         h_class = str_sub(h_class, start = 4),
         h_order = str_sub(h_order, start = 4),
         h_family = str_sub(h_family, start = 4),
         h_genus = str_sub(h_genus, start = 4))

saveRDS(iphop_genus, here('clean_data/iphop_genus_table.RDS'))
write.table(iphop_genus, file='clean_data/iphop_vOTU_to_genus.tsv', quote=FALSE,row.names = F, sep='\t')
# used this table to make the phage host network


iphop_genome <- read_csv(here("/projects/b1042/HartmannLab/lung_phage/iphop/derep_vOTUs/Host_prediction_to_genome_m90.csv")) %>% 
  rename(feature = 'Virus', 
         host_taxonomy = 'Host taxonomy',
         host_genome = 'Host genome',
         confidence_score = 'Confidence score') |> 
  separate(col = host_taxonomy, into = c("h_domain" ,"h_phylum", "h_class", "h_order", "h_family", "h_genus", "h_species"), sep = ';') %>% 
  mutate(h_domain = str_sub(h_domain, start = 4),
         h_phylum = str_sub(h_phylum, start = 4),
         h_class = str_sub(h_class, start = 4),
         h_order = str_sub(h_order, start = 4),
         h_family = str_sub(h_family, start = 4),
         h_genus = str_sub(h_genus, start = 4),
         h_species = str_sub(h_species, start = 4))


saveRDS(iphop_genome, here('clean_data/iphop_genome_table.RDS'))
write.table(iphop_genome, file='clean_data/iphop_vOTU_to_bin.tsv', quote=FALSE,row.names = F, sep='\t')




