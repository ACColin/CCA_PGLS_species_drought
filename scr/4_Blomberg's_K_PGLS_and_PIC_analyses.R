#' @title: 4: Phylogenetic signal with Pagel's lambda and Blomberg's K
#' @author: A.C. Colin
#' 
#' # Analysis decription:
#' *Background:* Some traits from Victoria's predictive models have been found to be good predictors of the impact of drought on trees at CCA. Those traits include functional traits measured in the lab, and life-history traits extracted from databases. From life-history traits, good predictors include species maximum height, resprouting strategy, aridity index and mean annual precipitation of the climate of origin. Wood density is the only functional trait found to predict drought impact at CCA.
#' *Objectives:* Here the objective is to control for the statistical dependence among traits resulting from evolution of traits along the phylogenetic tree. If there is a significant phylogenetic pattern, the objective is to perform phylogenetic independant contrasts for maximum height, aridity index,mean annual precipitation, wood density and drought impact to obtain corrected values to run in Victoria's models.
#' Here we use a pre-acquired dated ultrametric phylogeny of Eucalypt species that was updated according to the Classification of Eucalypts to assess phylogenetic relatedness of all species at the Currency Creek Arboretum
#' *Questions:*
#'   - Is variation in drought impact observed across all trees at CCA independant from the phylogenetic relationships of species?
#'   - Is there a significant phylogenetic signal for the traits described above?
#'   *Plan:* First, for all continuous traits, I will obtain Pagel's λ and Blomberg's K values of phylogenetic independance to determine whether phylogenetic signal is significant or not and perform a PIC analysis for each trait.
#' *Sidenote:* From what I found, the packages ape, phangorn, phytools, picante, caper, Geiger and phylolm implement methods for tests of phylogenetic dependence and trait:trait associations with phylogenetic correction.
#' 
#' 
#' # Getting ready:
#' Librarian will do everything for you (install, load, get dependencies, etc.). check out `??librarian` for info
#' To get librarian to work run those commands in the console:
#' `install.packages("librarian")`
#' `librarian::lib_startup(librarian, magrittr, lib = "C:/ACPrograms/R/R-3.6.3/library", global = TRUE)` #change to your library directory
#' 
librarian::shelf(ape, geiger, nlme, tidyverse, treeplyr, phytools, caper, tidytree)
#'
#' ## Use Blomberg’s K to assess the degree of phylogenetic signal
#' 
#' *Reference:* S. P. Blomberg, T. Garland Jr, A. R. Ives, Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution 57, 717–745 (2003).
#' K might be usefully thought of as a measure of the partitioning of variance. If K>1 then variance tends to be among clades; while if K<1 then variance is within clades (with BM as reference).
#' 
#' 
#' # Analysis:
#' The first step is to estimate and test the phylogenetic signal against a null model of no phylogenetic signal. Here we use Blomberg's K as test statistics for phylogenetic signal.
#' 
#' ## Load the data:
#' PGLS analyses do not deal with missing data. Here I prune the tree to collapse only to tips with data for the PGLS for each trait.
#' Note: there is a value for drought impact (DI) for each tip of the tree with Dean's survey at the arboretum so no need to prune for DI.

#' Read in trait and tree data
euctree = read.tree("plots/final_tree_uniqueID.tre") # has 1839 tips
eucdata = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv") # index from 0-4
eucdata = read.csv("data/DroughtImpact_SummaryBy_Completetraits_mod.csv") # index from 1-5
eucdata_D = read.csv("data/Phylo_AllCCAtraits.csv") # has 
### here include the dataframe from Victoria with a drought impact value for each individual tree in the phylogeny

#' Check if loaded correctly
#' We want everything to be true
str(euctree)
is.rooted(euctree)
is.ultrametric(euctree) # actually it is but not the way ape expects. don't mind ape, it's not always right

#this code is for doing a PGLS. keep for later, might be useful
# wingL <- geodata$wingL
# tarsusL <- geodata$tarsusL
# DF.geospiza <- data.frame(wingL,tarsusL,row.names=row.names(geodata))
# DF.geospiza <- DF.geospiza[geospiza13.tree$tip.label, ]
# DF.geospiza
# 
# bm.geospiza <- corBrownian(phy=geospiza13.tree)
# bm.gls <- gls(wingL~tarsusL,correlation=bm.geospiza,data=DF.geospiza)
# summary(bm.gls)
#' 
#' 
#' 
##############################################################
## Tests of phylogenetic signal for CCA drought sensitivity ##
##############################################################
#' 
#' Let's combine and match tree and data:
#' 
#' Create Impact dataset containing just drought impact values
Impact = read.csv("data/DroughtImpactSummary_data_for_Blombergs_final.csv", header = T) %>%
  dplyr::select(UniqueID, Drought2020Scale.mean)
#' Look at the first few rows
head(Impact)
#' 
#' Create the custom tree for analysis
K_euctree <- euctree
#' 
# Check whether the names match in the data and the tree
check <- name.check(phy = euctree, data = Impact, 
                    data.names = Impact$UniqueID)
check
#'
#' Combine and match the tree and data
eucstuff <- make.treedata(tree = euctree, data = Impact, 
                           name_column = "UniqueID")
#' 
#' Look at the tree
eucstuff$phy

# Look at the data
glimpse(eucstuff$dat)

# Make a new column called tiplabel with the tip labels in it
eucstuff$dat$tiplabel <- eucstuff$phy$tip.label
# Force mydata to be a data frame
mydata <- as.data.frame(eucstuff$dat)

# Create DI containing just drought impact values
DI <- pull(mydata, Drought2020Scale.mean)

# Notice that this is currently just a long list of numbers. We can then name these values with the species names from mydata using the function names. Note that this requires the trait data is in the same order as the tree tip labels, but luckily make.treedata does this automatically.

# Give DI names = species names at the tips of the phylogeny
names(DI) <- mydata$tiplabel
# Look at the first few rows
head(DI)

# Now we have a list of values with associated species names.

# Estimate Pagel's lambda
lambdaDI <- phylosig(K_euctree, DI, method = "lambda", test = T)

#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaDI

# > Phylogenetic signal lambda : 0.865945 
# > logL(lambda) : -1526.04 
# > LR(lambda=0) : 796.463 
# > P-value (based on LR test) : 3.16936e-175

# The λ estimate for DI is around 0.865.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).

# Here P < 0.001. We interpret this as λ being significantly different from 0, i.e. there is significant phylogenetic signal in drought sensitivity.

# Now let's estimate Blomberg’s K. To do so we also use phylosig but with method = K.
# Estimate Blomberg’s *K*

K_Impact <- phylosig(K_euctree, DI, method = "K", test = TRUE, nsim = 1000)
det(DI)

###################################################
## Test of phylogenetic signal for Aridity Index ##
###################################################
#'
# Get the data in and ready first:
euctree = read.tree("plots/final_tree_uniqueID.tre") # has 1839 tips
AI = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T) %>% 
  dplyr::select(UniqueID, AICorrected)
head(AI)
K_euctree = euctree
check = name.check(phy = euctree, data = AI, 
                    data.names = AI$UniqueID)
check
eucstuff = make.treedata(tree = euctree, data = AI, 
                          name_column = "UniqueID")
eucstuff$phy
glimpse(eucstuff$dat)
eucstuff$dat$tiplabel = eucstuff$phy$tip.label
mydata = as.data.frame(eucstuff$dat)
AI_data = pull(mydata, AICorrected)
names(AI_data) = mydata$tiplabel
head(AI_data)
lambdaAI = phylosig(K_euctree, AI_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaAI
# > Phylogenetic signal lambda : 0.954991 
# > logL(lambda) : 350.082 
# > LR(lambda=0) : 1094.28 
# > P-value (based on LR test) : 5.78454e-240 
# The λ estimate for DI is ~0.955.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in Aridity Index.

K_Impact <- phylosig(K_euctree, AI_data, method = "K", test = TRUE, nsim = 1000)
det(AI_data)

####################################################
## Test of phylogenetic signal for Maximum height ##
####################################################
#'
# Get the data in and ready first:
euctree = read.tree("plots/final_tree_uniqueID.tre") # has 1839 tips
max_h = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T) %>% 
  dplyr::select(UniqueID, height_max_m)
head(max_h)
K_euctree = euctree
check = name.check(phy = euctree, data = max_h, 
                   data.names = max_h$UniqueID)
check
eucstuff = make.treedata(tree = euctree, data = max_h, 
                         name_column = "UniqueID")
eucstuff$phy
glimpse(eucstuff$dat)
eucstuff$dat$tiplabel = eucstuff$phy$tip.label
mydata = as.data.frame(eucstuff$dat)
MH_data = pull(mydata, height_max_m)
names(MH_data) = mydata$tiplabel
head(MH_data)
lambdaMH = phylosig(K_euctree, MH_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMH
# > Phylogenetic signal lambda : 0.990818 
# > logL(lambda) : -4885.21 
# > LR(lambda=0) : 850.542
# > P-value (based on LR test) : 5.54273e-187
# The λ estimate for DI is ~0.991.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact <- phylosig(K_euctree, MH_data, method = "K", test = TRUE, nsim = 1000)
det(AI_data)


###############################################################
## Test of phylogenetic signal for Mean Annual Precipitation ##
###############################################################
#'
# Get the data in and ready first:
euctree = read.tree("plots/final_tree_uniqueID.tre") # has 1839 tips
MAT = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T) %>% 
  dplyr::select(UniqueID, MAT)
head(MAT)
K_euctree = euctree
check = name.check(phy = euctree, data = MAT, 
                   data.names = MAT$UniqueID)
check
eucstuff = make.treedata(tree = euctree, data = MAT, 
                         name_column = "UniqueID")
eucstuff$phy
glimpse(eucstuff$dat)
eucstuff$dat$tiplabel = eucstuff$phy$tip.label
mydata = as.data.frame(eucstuff$dat)
MAT_data = pull(mydata, MAP)
names(MAT_data) = mydata$tiplabel
head(MAT_data)
lambdaMAT = phylosig(K_euctree, MAT_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMAT
# > Phylogenetic signal lambda : 0.990818 
# > logL(lambda) : -4885.21 
# > LR(lambda=0) : 850.542
# > P-value (based on LR test) : 5.54273e-187
# The λ estimate for DI is ~0.991.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact <- phylosig(K_euctree, MH_data, method = "K", test = TRUE, nsim = 1000)
det(AI_data)
