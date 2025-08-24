# phage host pairing
# parse iphop data for network analysis

# load libraries ----
library(tidyverse)
library(here)
library(skimr)

# load data ----
iphop_genome <- readRDS(here("clean_data/iphop_genome_table.RDS"))
iphop_genus <- readRDS(here("clean_data/iphop_genus_table.RDS"))


# filter out host predictions made by reference database
iphop_genome_bin_filtered <- iphop_genome %>% 
  filter(str_detect(host_genome, 'SemiBin')) 

# no species level identification my from bins
skim(iphop_genome_bin_filtered)

# iphop genus table split and compare phage-based vs host-based identifications

host_based_id <- iphop_genus %>% 
  filter(!str_detect(methods, 'RaFAH'))

phage_based_id <- iphop_genus %>% 
  filter(str_detect(methods, 'RaFAH'))

# try pivot wider instead 

iphop_genus %>% 
  select(feature, methods) %>% 
  mutate(rafah = str_detect(methods, 'RaFAH'),
         host_based = !str_detect(methods, 'RaFAH')) %>% 



