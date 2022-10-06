#' @title Seedlot-level Analyses - PICs and PGLS
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

#' PREPARING PACKAGES:
librarian::shelf(ape, geiger, tidyverse, treeplyr, caper, ggpubr)
#' 
#' 
#' 
#' READ IN DATA:
eucdata = read.csv("data/DroughtImpactSummary_data_for_Blombergs_final.csv")
glimpse(eucdata)
eucdata$UniqueID = gsub("^", "PGLS_", eucdata$UniqueID) # I have to do this otherwise comparative.data() is angry at me for having duplicate names between node labels and tip labels.
glimpse(eucdata) # has 1839 observations
euctree = read.tree("plots/final_tree_with_CCAUniqueIDs_instead_of_UniqueIDs_dichotomous.tre")
euctree$tip.label = gsub("^CCA", "PGLS_", euctree$tip.label) #same reason, they have to match between data and tree now
euctree$tip.label = gsub("^A", "PGLS_A", euctree$tip.label)
euctree$tip.label = gsub("^E", "PGLS_E", euctree$tip.label)

# euctree$tip.label = gsub("^PGLS_E", "E", euctree$tip.label)
# euctree$tip.label = gsub("^PGLS_A", "A", euctree$tip.label)
str(euctree) # has 1839 tips
# Remember to check the tree is dichotomous, i.e. has no polytomies, rooted, and ultrametric.
is.binary.tree(euctree) # [1] TRUE
is.rooted(euctree)      # [1] TRUE
is.ultrametric(euctree) # [1] FALSE
#' 
#' We need to adjust the root to tip branch lengths for the tree to be ultrametric:
#' Let's use the variance method to check the tree ultrametricity
N = Ntip(euctree)
root_node = N + 1 
root_to_tip = dist.nodes(euctree)[1:N, root_node] # compute root to tip distances for all tips 
min_tip = min(root_to_tip)
max_tip = max(root_to_tip)
(max_tip - min_tip) / max_tip # [1] 6.482982e-08
                              # considering the tolerance in ape is 1e-08,
                              # immediately we see why the tree stopped being ultrametric
                              # this relative difference is larger than the default tolerance
# Let's scale the root-to-tip distance to see what it looks like:
scaled_root_to_tip = root_to_tip * 1000
var(scaled_root_to_tip)                 # [1] 9.570864e-07
min_tip = min(scaled_root_to_tip)
max_tip = max(scaled_root_to_tip)
(max_tip - min_tip) / max_tip           # [1] 6.482982e-08
# even though none of the branch lengths have changed relative to each other
# the tree is no longer ultrametric per the variance statistic
#' 
# One solution for trees that are not quite there,
# is to extend the tips of the tree until the root-to-tip distances are completely equal.

# This is implemented in R as `BioGeoBEARS::extend_tips_to_ultrametricize`
# and `phytools::force.ultrametric(method = "extend")`

tre_extend = euctree                                  # have a copy of the tree
age_difference = max(root_to_tip) - root_to_tip       # compute the diff. from each root-to-tip distance to their max
tip_edges = tre_extend$edge[, 2] <= Ntip(tre_extend)  # grab the edges from matrix that corresponds to tips
# edges in $edge.label corresponds to the row numbers in $edge

tre_extend$edge.length[tip_edges] = tre_extend$edge.length[tip_edges] + age_difference
is.ultrametric(tre_extend) # [1] TRUE (yeyy!)
is.rooted(tre_extend)      # [1] TRUE
is.binary(tre_extend)      # [1] TRUE
#'
#'
#'
diff_edge_lengths = function(phy, phy2) { # function to compare two phylogenies with identical topologies
  # but differing branch lengths
  diffs = phy2$edge.length - phy$edge.length
  cols = sign(diffs)
  cols[cols == 1] = "#7fbc41"
  cols[cols == -1] = "#de77ae"
  cols[cols == 0] = NA
  plot(phy, show.tip.label = FALSE, no.margin = TRUE)
  edgelabels(pch = 15, col = cols)
  sprintf("%i longer branches, %i shorter branches", sum(diffs > 0), sum(diffs < 0))
}

diff_edge_lengths(euctree, tre_extend) # [1] "1838 longer branches, 0 shorter branches"
# using this fix basically increased the size of almost all the final
# branches lengths which were problematic especially in the very young clades
#' 
#' 
#' 
#' Now let's check if the species names match up in the tree and the data.
#' This should reveal any typos and taXonomic missmatches. In our case, it is the UniqueID so we shouldn't have any typo issues.
check = name.check(phy = tre_extend, data = eucdata,
                    data.names = eucdata$UniqueID)
check
#' Only Egilii1 is missing so that's alright.
#' 
#' Let's combine and match the tree and data now:
eucstuff = make.treedata(tree = euctree, data = eucdata,
                          name_column = "UniqueID")
#' Check out the tree and data:
eucstuff$phy
glimpse(eucstuff$dat)
#' Noice
#' 
#' Make a new column called tiplabel with the tip labels in it
eucstuff$dat$tiplabel = eucstuff$phy$tip.label
#' Force mydata to be a data frame
mydata = as.data.frame(eucstuff$dat)
summary(mydata)
#' Summary for `Drought2020Scale.mean`:
  # Drought2020Scale.mean
  # Min.   :0.0000
  # 1st Qu.:0.0000
  # Median :0.0000
  # Mean   :0.5577
  # 3rd Qu.:0.7500
  # Max.   :4.0000
  # NA's   :481
#' Summary for `Wooddensity.g.cm3..mean`:
  # Wooddensity.g.cm3..mean
  # Min.   :0.5100
  # 1st Qu.:0.6375
  # Median :0.6858
  # Mean   :0.6815
  # 3rd Qu.:0.7300
  # Max.   :0.9050
  # NA's   :1620

#' And the tree
mytree = eucstuff$phy
#' 
#' 
#' Now we are ready for the analysis!
#' 
#' Let's start by investigating the relationship between drought impact and wood density. From Victoria we know that that they are highly correlated already.
#' 
#' Plot drought impact against WD, coloured by habit, I wanted by section but it's not in it.
#' 
fit = glm(mydata$height_max_m~mydata$AICorrected)
co = coef(fit)
ggplot(data = mydata, aes(x = log(height_max_m),
                          y = log(AICorrected),
                          colour = Habit)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  facet_wrap(~Habit) +
  theme_bw()

ggplot(data = mydata, aes(x = log(Drought2020Scale.mean),
                          y = log(Wooddensity.g.cm3..mean),
                          colour = Habit)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +
  facet_wrap(~Habit) +
  theme_bw()


#' 
#' Unfortunately without the taxonomic section we can't see that close relatives are more similar than distant relatives but let's pretend that we can. Also, the results of Pagel's λ statistical test revealed that drought impact is strongly significantly influenced by phylogenetic signal.
#' We need to account for phylogenetic non-independance, because of the statistical issues caused by it and also because there is a better way to model the biological reality of this.
#' We know that euc species evolve from other euc species, and that close relatives for evolutionary reasons will therefore be similar, so we should add this into our models!
#' There are several ways of accounting for phylogenetic non-independence in your analyses. Here we will use phylogenetic generalized least squares (PGLS). Another popular earlier method is independent contrasts (PIC). This method is really similar to PGLS, in fact it is just a special kind of PGLS where λ is equal to 1.
#' PGLS offers some important advantages over independent contrasts. The model of trait evolution can be more flexible i.e., it can depart from a strict Brownian motion process (λ or K = 1).
#' Different scaling parameters (λ,  κ, and  δ) can be incorporated in the analysis, which can significantly improve the fit of the data to the model and thus also improve the estimation of the trait correlation. Another advantage of PGLS is that the intercept of the regression is not forced to be zero. See the Primer for more details on the theory underlying PICs and PGLS.
#' 
#' *Fitting the PGLS*
#' To perform PGLS models in R, caper requires you to first combine the phylogeny and data into one object using the function comparative.data. This is similar to what we did with make.treedata, but it does some stuff that is particular to how caper works so we still need to do this here.
#' 
#' Note that vcv = TRUE stores a variance covariance matrix of your tree (you will need this for the pgls function). na.omit = FALSE stops the function from removing species without data for all variables. warn.dropped = TRUE will tell you if any species are not in both the tree and the data and are therefore dropped from the comparative data object. Here we won’t drop any species because we already did this using make.treedata.
#' 
#' 
mytree$node.label=NULL # IMPORTANT: remove node labels otherwise comparative.data throws error
                       # Labels duplicated between tips and nodes in phylogeny
euc = comparative.data(phy = mytree, data = mydata, 
                        names.col = tiplabel, vcv = TRUE, 
                        na.omit = FALSE, warn.dropped = TRUE)
str(euc)
# If you do need to drop species, this function will give a warning telling you that some species have been dropped. You can view the dropped species using:
euc$dropped$tips # character(0)
euc$dropped$unmatched.rows # character(0)
# Always make sure you check the list of dropped species is what you expected, it often reveals typos in your species names, or mismatches in taxonomies used etc. 
#' 
#' 
# The function for PGLS analyses in caper is pgls. To fit a model which uses the Maximum Likelihood (ML) estimate of λ we use the following code:
# Fit a PGLS model
model.glm = glm(log(Drought2020Scale.mean)~log(AICorrected), data = euc, family = )
model.pgls = pgls(log(Drought2020Scale.mean) ~ log(AICorrected),
                  data = euc, lambda = "ML",
                  bounds = list(lambda=c(0.001,1), kappa=c(1e-6,3), delta=c(1e-6,3)))
#' I keep having the following errors:
#' Error in optim(optimPar, fn = pgls.likelihood, method = "L-BFGS-B", control = control,  : L-BFGS-B needs finite values of 'fn'
#' When the variables are log transformed
#' and Error in solve.default(V, tol = .Machine$double.eps) : system is computationally singular: reciprocal condition number = 7.16501e-18
#' when they are not log transformed
#' 
librarian::shelf(fitdistrplus)
fitdis
x = mydata$Drought2020Scale.mean %>% 
  na.omit() %>% 
  as.vector()
normal_dist = fitdist(x, "norm")
plot(normal_dist)
descdist(x)
descdist(x, boot = 1000)
density.x = density.default(log(x), bw = 0.5)
plot(density.x, xlab = "N Bandwidth = 0.5", ylab = " Density")
