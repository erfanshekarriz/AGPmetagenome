# renv:install.packages("tidyverse")
# renv:install.packages("phyloseq")

library(tidyverse)
library(phyloseq)


df <- read.csv("data/tabular/metaphlan_out_SGB.txt", sep = "\t", comment.char = "#")

otudf <- df %>% 
  mutate(SGB_id = str_extract(clade_name, "t__SGB.+$"), 
         SGB_id = str_remove(SGB_id, "t__")) %>% 
  filter(!str_detect(clade_name, "Euk")) %>% 
  column_to_rownames(var = "SGB_id")

taxdf <- otudf %>% 
  select(clade_name) %>% 
  separate(
    col = clade_name,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
    sep = "\\|",
    fill = "right" # Fill missing values with NA
  ) %>% 
  mutate_all(~str_remove_all(., "^.*__")) %>% 
  mutate(Species = str_replace(Species, "_", " ")) 


otumat <- otudf %>% 
  select(-clade_name) %>% 
  as.matrix()

taxmat <- as.matrix(taxdf)


OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX)
saveRDS(physeq, "data/rds/physeq.rds")

