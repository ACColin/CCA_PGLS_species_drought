library(ape)
library(tidyverse)
library(ggtree)

metadata <- read.csv("../../../phyloGWAS CurrencyCreek/36BUSCO_phylo/AC_BUSCO_samples_metadata.csv", header = T)
set1_tree <- read.tree("../../../../set1_species.tree")
set2_tree <- read.tree("../../../../set2_species.tree")
set3_tree <- read.tree("../../../../set3_species.tree")

set1 <- ggtree(set1_tree, layout='circular') %<+% metadata + 
  geom_tippoint(aes(fill=Series, color="green")) +
  geom_tiplab(aes(color=Series), hjust = -.1)
print(set1)
