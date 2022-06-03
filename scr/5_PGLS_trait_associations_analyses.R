#' @title 5: Phylogenetic Generalised Least Squares (PGLS) statistical tests
#' @author A.C. Colin
#' 
#' 
#' # Analysis decription:
#' PGLS is a tool for trait:trait associations that implements GLS with residual covariances defined by a model of evolution.
#' *Background:* Some traits from Victoria's predictive models have been found to be good predictors of the impact of drought on trees at CCA. Those traits include functional traits measured in the lab and life-history traits extracted from databases. From life-history traits, good predictors include species maximum height, resprouting strategy, aridity index and mean annual precipitation of the climate of origin. Wood density is the only functional trait found to predict drought impact at CCA.
#' *Objective:* Victoria has a a bunch of good predictors of drought tolerance at CCA. Here the objective is to look at the association of those traits with the drought impact at CCA in 2020, to see if after fitting a GLS model, the residual covariances follow a model of evolution.
#' 
#' *Questions:*
#'   - Accounting for phylogenetic signal, do we observe the same trait associations between the drought incidence at CCA and the other variables across all species?
#'  *Plan:* I perform a Phylogenetic Generalized Least Squares analysis on drought impact against all other traits across all trees at CCA, accounting for phylogenetic signal using the modified tree made previously from Thornhill et al 2016.

#' *Preparation of 
#' First let's load what we need:
librarian::shelf(ape, geiger, tidyverse, treeplyr, caper)
#' 
#' Read in the data:
eucdata <- read.csv("data/DroughtImpactSummary_data_for_Blombergs_final.csv")
eucdata$UniqueID <- gsub("^", "PGLS_", eucdata$UniqueID) # I have to do this otherwise comparative.data() is angry with me for having duplicate names between node labels and tip labels.
glimpse(eucdata) # has 1839 observations
euctree <- read.tree("plots/final_tree_uniqueID.tre")
euctree$tip.label <- gsub("^", "PGLS_", euctree$tip.label) #same reason, they have to match between data and tree now
str(euctree) # has 1839 tips
#' 
#' Now let's check if the species names match up in the tree and the data.
#' This should reveal any typos and taoxnomic missmatches. In our case, it is the UniqueID so we shouldn't have any typo issues.
check <- name.check(phy = euctree, data = eucdata,
                    data.names = eucdata$UniqueID)
check
#' Only Egilii1 is missing so that's alright.
#' 
#' Let's combine and match the tree and data now:
eucstuff <- make.treedata(tree = euctree, data = eucdata,
                          name_column = "UniqueID")
#' Check out the tree and data:
eucstuff$phy
glimpse(eucstuff$dat)
#' Noice
#' 
#' Make a new column called tiplabel with the tip labels in it
eucstuff$dat$tiplabel <- eucstuff$phy$tip.label
#' Force mydata to be a data frame
mydata <- as.data.frame(eucstuff$dat)
#' And the tree
mytree <- eucstuff$phy
#' 
#' 
#' Now we are ready for the analysis!
#' 
#' Let's start by investigating the relationship between drought impact and wood density. From Victoria we know that that they are highly correlated already.
#' 
#' Plot drought impact against WD, coloured by habit, I wanted by section but it's not in it.
#' 
ggplot(data = mydata, aes(x = Drought2020Scale.mean,
                          y = Wooddensity.g.cm3..mean,
                          colour = Habit)) +
  geom_point() +
  theme_bw()
#' 
#' Unfortunately without the taxonomic section we can't see that close relatives are more similar than distant relatives but let's pretend that we can. Also, the results of Pagel's λ statistical test revealed that drought impact is strongly significantly influenced by phylogenetic signal.
#' We need to account for phylogenetic non-independance, because of the statistical issues caused by it and also because there is a better way to model the biological reality of this.
#' We know that euc species evolve from other euc species, and that close relatives for evolutionary reasons will therefore be similar, so we should add this into our models!
#' There are several ways of accounting for phylogenetic non-independence in your analyses. Here we will use phylogenetic generalized least squares (PGLS). Another popular earlier method is independent contrasts (PIC). This method is really similar to PGLS, in fact it is just a special kind of PGLS where   λ is equal to 1.
#' PGLS offers some important advantages over independent contrasts. The model of trait evolution can be more flexible i.e., it can depart from a strict Brownian motion process ( λ or K = 1). Different scaling parameters ( λ,  κ, and  δ) can be incorporated in the analysis, which can significantly improve the fit of the data to the model and thus also improve the estimation of the trait correlation. Another advantage of PGLS is that the intercept of the regression is not forced to be zero. See the Primer for more details on the theory underlying PICs and PGLS.
#' 
#' *Fitting the PGLS*
#' To perform PGLS models in R, caper requires you to first combine the phylogeny and data into one object using the function comparative.data. This is similar to what we did with make.treedata, but it does some stuff that is particular to how caper works so we still need to do this here.
#' 
#' Note that vcv = TRUE stores a variance covariance matrix of your tree (you will need this for the pgls function). na.omit = FALSE stops the function from removing species without data for all variables. warn.dropped = TRUE will tell you if any species are not in both the tree and the data and are therefore dropped from the comparative data object. Here we won’t drop any species because we already did this using make.treedata.
#' 
#' 
euc <- comparative.data(phy = mytree, data = mydata, 
                        names.col = tiplabel, vcv = TRUE, 
                        na.omit = FALSE, warn.dropped = TRUE)
