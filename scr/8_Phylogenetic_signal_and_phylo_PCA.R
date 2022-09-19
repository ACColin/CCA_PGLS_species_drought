librarian::shelf(adephylo, tidyverse, ape, ggtree, treeplyr, phylobase, geiger)

## Multivariate analysis with phylogenetic principal component analysis (pPCA) with ppca {adephylo}
euc_tree = read.tree("plots/final_tree_uniqueID.tre") 
plot(euc_tree)
euc_dat = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv",
                   header = T) %>% 
  dplyr::select(UniqueID, Drought2020Scale.mean, AICorrected, AP, height_max_m) %>% 
  drop_na() %>% 
  filter(Drought2020Scale.mean !="") %>% 
  distinct(UniqueID, .keep_all = T)
head(euc_dat)
rownames(euc_dat) = NULL
rownames(euc_dat) = euc_dat$UniqueID
pPCA_euctree = euc_tree # custom tree for analysis
check = name.check(phy = pPCA_euctree, data = euc_dat, # Check whether the names match in the data and the tree
                    data.names = pPCA_euctree$tip.label)
check
eucstuff = treeplyr::make.treedata(tree = euc_tree, data = euc_dat, # Combine and match the tree and data
                          name_column = "UniqueID")
eucstuff$phy # look at tree
glimpse(eucstuff$dat) # Look at the data
eucstuff$dat$tiplabel = eucstuff$phy$tip.label # Make a new column called tiplabel with the tip labels in it
mydata = as.data.frame(eucstuff$dat) # Force mydata to be a data frame
# Now we have a list of values with associated seedlot.
dat = phylo4d(euc_tree, mydata, missing.data = "OK")
res = ppca(dat, scannf = T, nfposi=1, nfnega=1, method="Abouheif")

## Phylogenetic pattern of quantitative traits with orthogram {adephylo}
# find code from p. 214 of Emmanuel Paradis book on APE with R
adephylo::orthogram(trait, tree)
