## Changing the tree with appropriate tip labels:

librarian::shelf(ape, ggtree, geiger, tidyverse, phylobase, phytools, treeplyr)

thornhill_labels <- read.csv("data/tree_labels.csv") # TIP LABELS DATA
euc_tree = read.nexus("data/Eucalypts_ML2_dated_r8s.tre") # ORIGINAL TREE DATA
euc_data = read.csv(file = "data/DroughtImpact_SummaryBy_Completetraits.csv") # TRAITS AND ENV DATA

euc_tree$tip.label = thornhill_labels$name5[TreeTools::match(euc_tree$tip.label, thornhill_labels$original)]
par(mfrow=c(1,2))
plot(euc_tree)

