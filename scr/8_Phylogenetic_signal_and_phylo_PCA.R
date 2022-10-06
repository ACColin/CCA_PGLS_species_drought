librarian::shelf(adephylo, tidyverse, ape, ggtree, treeplyr, phylobase, geiger, TreeTools, phytools)

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

# Warning messages:
#  1: In checkTree(object) : Labels are not unique.
#  2: In checkTree(object) : Labels are not unique.
# There are no duplicates in the tip labels. Each unique ID appears only once in the tree.
# I get this error message with phylo4d which I think is caused because the tip labels (UniqueIDs) are numbers like the nodes. I just added "CCA" in front of each tip label for the tree and for the data see here:

CCAeuc_tree = read.tree("plots/final_tree_uniqueID.tre") # now let's swap the UniqueID with CCAUniqueID
CCCA_euc_tree_labels = read.csv("data/final_tree_UniqueID_tip_labels_as_df.csv", header = T)
glimpse(CCCA_euc_tree_labels)
new_CCAtree = CCAeuc_tree
new_CCAtree$tip.label = CCCA_euc_tree_labels$tree_CCAUniqueIDs[TreeTools::match(CCAeuc_tree$tip.label, CCCA_euc_tree_labels$tree_UniqueIDs)]
par(mfrow=c(1,2))
plot(new_CCAtree)
# LET'S TRY AGAIN NOW
mydata = read.csv("data/UniqueID_tree_tips_data_for_ASR.csv", header = T)
dat = phylo4d(new_CCAtree, mydata, missing.data = "OK")

# Side note I am going to resolve polytomies on the tree, because it is not dichotomous like that it will probably solve some issues and I know it can be problematic.
# The {ape} package has a function called multi2di to resolve polytomies by adding branches of zero length (while preserving the mappings) in a tree or set of trees. Let's try it.

di_CCAeuc_tree = ape::multi2di(new_CCAtree)
is.ultrametric(di_CCAeuc_tree)
is.ultrametric(di_CCAeuc_tree) #great..
plot(di_CCAeuc_tree)

# maybe now I can fix the tree with the tutorial I found...
# >>>>>>>>>>>>>>>>>>>>>
N = Ntip(new_CCAtree)
root_node = N + 1 
root_to_tip = dist.nodes(new_CCAtree)[1:N, root_node] # compute root to tip distances for all tips 
var(root_to_tip)                                      # compute variance using those distances
# [1] 9.670593e-13
sqrt(.Machine$double.eps) # compare to default tolerance in ape (sqrt of R's machine epsilon)
                          # to be considered ultrametric,
                          # the variance has to be smaller than the tolerance in ape
# [1] 1.490116e-08
# 1e-08 > 9e-13 therefore,
# the euc tree is ultrametric according to the variance method
tre_extend = new_CCAtree                                  # have a copy of the tree
age_difference = max(root_to_tip) - root_to_tip       # compute the diff. from each root-to-tip distance to their max
tip_edges = tre_extend$edge[, 2] <= Ntip(tre_extend)  # grab the edges from matrix that corresponds to tips
# edges in $edge.label corresponds to the row numbers in $edge

tre_extend$edge.length[tip_edges] = tre_extend$edge.length[tip_edges] + age_difference
is.ultrametric(tre_extend) # [1] TRUE (yeyy!)
# <<<<<<<<<<<<<<<<<<<<<
is.rooted(tre_extend)
is.ultrametric(tre_extend) #great..
write.tree(tre_extend, "plots/final_tree_CCAUniqueID_tips_dichotomous_ultra-extended.tre")

## Phylogenetic pattern of quantitative traits with orthogram {adephylo}
# find code from p. 214 of Emmanuel Paradis book on APE with R

dat = phylo4d(tre_extend, mydata, missing.data = "OK")
res = ppca(dat)
adephylo::orthogram(mydata, tree)

### PGLS
bm.euc = corBrownian(phy = tre_extend)
DF.euc = data.frame(mydata)
library(nlme)
m1= gls(AICorrected ~ Drought2020Scale.mean, DF.euc, correlation = bm.euc)
summary(m1)
ou.euc = corMartins(1, tre_extend)
m2 = gls(AICorrected ~ Drought2020Scale.mean, DF.euc, correlation = ou.euc)
summary(m2)
