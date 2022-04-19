### Renaming my trees ###

library(ape)
library(tidyverse)
library(TreeTools)

# The ASTRAL tree leaves the branch length of terminal branches empty.
# Some tools for visualization and tree editing do not like this (e.g., ape).
# In FigTree, if you open the tree several times, it eventually opens up (at least on our machines).
# In ape, if you ask it to ignore branch lengths all together, it works.
# In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, as in this script: "scr/add-bl.py"


ID_table <-read.csv("../../../../tree_ID_table_reduced.csv")
tree1 <- read.tree("../../../../set1_species.tree")
tree1$tip.label <- ID_table$Binomial[TreeTools::match(tree1$tip.label, ID_table$FieldID)]
class(tree1)

plot.phylo(tree1)

write.tree(tree1, "../../../../set1_speciesnames.tree")


tree2 <- read.tree("../../../../set2_species.tree")
tree2$tip.label <- ID_table$Binomial[TreeTools::match(tree2$tip.label, ID_table$FieldID)]
plot.phylo(tree2[[2]])
write.tree(tree2, "../../../../set2_speciesnames.tree")


tree3 <- read.tree("../../../../set3_species.tree")
tree3$tip.label <- ID_table$Binomial[TreeTools::match(tree3$tip.label, ID_table$FieldID)]
plot.phylo(tree3[[2]])
write.tree(tree3, "../../../../set3_speciesnames.tree")



ID_table <-read.csv("../../../../tree_ID_table_reduced.csv")
tree1 <- read.tree("../../../../set2_RAxML-style_partition.txt_reduced.treefile")
tree1$tip.label <- ID_table$Binomial[TreeTools::match(tree1$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree1,"../../../../set2_RAxML-style_partition.txt_reduced_spnames.treefile")

plot.phylo(tree1)

tree2 <- read.tree("../../../../set3_RAxML-style_partition.txt_reduced.treefile")
tree2$tip.label <- ID_table$Binomial[TreeTools::match(tree2$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree2,"../../../../set3_RAxML-style_partition.txt_reduced_spnames.treefile")


tree3 <- read.tree("../../../../set1_RAxML-style_partition.txt_reduced.treefile")
tree3$tip.label <- ID_table$Binomial[TreeTools::match(tree3$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree3,"../../../../set1_RAxML-style_partition.txt_reduced_spnames.treefile")


tree2 <- read.tree("../../../../loci_all_BUSCOs.treefile")
tree2$tip.label <- ID_table$Binomial[TreeTools::match(tree2$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree2,"../../../../loci_all_BUSCOs_MARE-IQTREE_spnames.treefile")

tree3 <- read.tree("../../../../all_BUSCOs_concat.treefile")
tree3$tip.label <- ID_table$Binomial[TreeTools::match(tree3$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree3,"../../../../concat_all_BUSCOs_MARE-IQTREE_spnames.treefile")

tree3 <- read.tree("../../../../all_BUSCOs_ASTRAL.treefile.trees")
tree3$tip.label <- ID_table$Binomial[TreeTools::match(tree3$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree3,"../../../../all_BUSCOs_ASTRAL_spnames.treefile")

tree3 <- read.tree("../../../../supertree.ufbt.treefile")
tree3$tip.label <- ID_table$Binomial[TreeTools::match(tree3$tip.label, ID_table$FieldID)]
class(tree1)
write.tree(tree3,"../../../../supertree.ufbt_spnames.treefile")
