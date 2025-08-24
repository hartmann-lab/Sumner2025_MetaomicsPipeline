# MGX binning import data 

# load libraries ----
library(tidyverse)
library(here)
library(skimr)

# import checkm results ---
get_checkm <- function(sample){
  checkm_summary <- read_tsv(sample, col_names = TRUE) |> 
    mutate(sample = sample) %>% 
    separate(col=sample, into=c(NA,NA,NA,NA,NA,NA,NA,"sample", NA), sep="/") 
  return(checkm_summary)
}

checkm_file_list <- list.files("/projects/b1042/HartmannLab/lung_phage/checkm2/", 
                                pattern="quality_report.tsv", 
                                recursive = TRUE,
                                full.names=TRUE)

checkm_table_list <- lapply(checkm_file_list, get_checkm)

checkm_fulltable <- bind_rows(lapply(checkm_table_list, function(x) if(nrow(x) == 0) NULL else x))

skimr::skim(checkm_fulltable)

count(checkm_fulltable, Completeness >= 90, Contamination <= 10)
count(checkm_fulltable, Completeness == 100)

write.table(checkm_fulltable, file='clean_data/checkm_fulltable.tsv', quote=FALSE,row.names = F, sep='\t')

# Filter quality bins

hq_bins <- checkm_fulltable %>% 
  filter(Completeness >= 90, Contamination <= 5) |>
  mutate(bin_id = paste0(sample, "_HQ_", Name)) |>
  mutate(bin_id = str_sub(bin_id, end = -4)) %>% 
  select(Name, bin_id, sample)

mq_bins <- checkm_fulltable %>% 
  filter(Completeness >= 50, Completeness < 90, Contamination <= 10) |>
  mutate(bin_id = paste0(sample, "_MQ_", Name)) |>
  mutate(bin_id = str_sub(bin_id, end = -4)) %>% 
  select(Name, bin_id, sample)

count(checkm_fulltable, str_detect(sample, "ZYMO"))

skim(qual_bins)

write.table(qual_bins, file='clean_data/semibin2_qual_bins.tsv', quote=FALSE,row.names = F, sep='\t')

# use the qual bin table to parse good quality bins into a separate folder using awk:
#awk '{ print "cp /projects/b1042/HartmannLab/lung_phage/semibin2/" $3 "/output_bins/" $1 ".gz \
#/projects/b1042/HartmannLab/lung_phage/semibin2/1_quality_bins/"$2".fa.gz"}' \
#/projects/b1042/HartmannLab/lung_phage/Rproject/bas_virus/clean_data/semibin2_qual_bins.tsv | bash
