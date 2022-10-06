

########################## Setup - Load Packages ##########################

librarian::shelf(ape, ggtree, geiger, tidyverse, phylobase, phytools, treeplyr, caper, ggpubr, MASS, countreg, adephylo)


########################## Setup - Reading Data ##########################

thornhill_labels = read.csv("data/tree_labels.csv") %>% # TIP LABELS DATA
  glimpse()
UniqueID_to_RowTree = read.csv("data/CCA_all_arboretum_rawdata.csv", header = T) %>% 
  glimpse()
taxa_to_ID_labels = read.csv("data/NewTipsMatrixUniqueID.csv") %>%  # SWAP BETWEEN TAXO NAMES AND UNIQUEIDS 
  glimpse()
euc_tree = read.tree("plots/final_tree_species_tips_ultra-extended.tre") %>% # ORIGINAL TREE DATA
  glimpse() # has labels as "kjellbergiodendron_celebicum" (lowercase and underscore for \s)
euc_data = read.csv(file = "data/DroughtImpact_SummaryBy_Completetraits.csv") %>% # TRAITS AND ENV DATA
  glimpse()
drought_data = read.csv("outputs/CCADroughtdataMarch_2020_DN_AC_percentage_aff_trees.csv", header = T) %>%
  dplyr::select(UniqueID, percentage_affected) %>%
  dplyr::distinct(UniqueID, .keep_all = T) %>% 
  glimpse()
LHT_data = read.csv("data/life_history_CCA.csv", header = T) %>% 
  glimpse()
first = left_join(LHT_data, UniqueID_to_RowTree, "RowTree")
first$DNNumber = as.numeric(first$DNNumber)
drought_data$UniqueID = as.numeric(drought_data$UniqueID)
second = left_join(first, drought_data, c("DNNumber" = "UniqueID"))
euc_data$UniqueID = as.double(euc_data$UniqueID)
third = left_join(euc_data, second, c("UniqueID" = "DNNumber"))
euc_data_all = read.csv("outputs/CCADroughtdataMarch_2020_DN_AC_percentage_aff_trees_shifted.csv") %>% 
  dplyr::select(- per_aff_spavg) %>% 
  rename(per_aff_spavg = per_aff_spavg_shift) %>% 
  glimpse()
fourth = left_join(third, euc_data_all, "UniqueID") %>% 
  distinct(UniqueID, .keep_all = T)
euc_data_all = fourth %>% 
  dplyr::select(-contains(c('.x.x', '.y', '.x.y'))) %>% 
  dplyr::distinct(UniqueID, .keep_all = T) %>% 
  dplyr::distinct(Thornhill_bionomial_tip_match2.x, .keep_all = T)


#euc_data_all = third %>% 
#  distinct(UniqueID, .keep_all = T) %>% 
#  group_by(Thornhill_bionomial_tip_match2) %>%
#  dplyr::select(UniqueID, Yearplanted.x, Taxonsimplified_DN, Provenance.x,
#                latitude, longitude, AICorrected, MAT, ISOT,
#                AP, PDryM, MaxTWarmM, PrecDeficit, per_aff_spavg_shift, 
#                height_max_m, LeafLenghtavg_mm, LeafWidth_avg_mm, LeafAreaEst_cm2,
#                Thornhill_bionomial_tip_match2, percentage_affected) %>% 
#  dplyr::filter(!is.na(percentage_affected)) %>% 
#  mutate(
#    per_aff_spavg = mean(percentage_affected),
#    AICorrected_spavg = mean(AICorrected),
#    MAT_spavg = mean(MAT),
#    MAP_spavg = mean(AP),
#    PDryM_spavg = mean(PDryM),
#    MaxTWarmM_spavg = mean(MaxTWarmM),
#    ISOT = mean(ISOT),
#    PrecDeficit_spavg = mean(PrecDeficit),
#    height_max_m_spavg = mean(height_max_m),
#    LeafLenghtavg_mm_spavg = mean(LeafLenghtavg_mm),
#    LeafWidth_avg_mm_spavg = mean(LeafWidth_avg_mm),
#    LeafAreaEst_cm2_spavg = mean(LeafAreaEst_cm2)) %>% 
#  distinct(Thornhill_bionomial_tip_match2, .keep_all = T) %>% 
#  glimpse()

duplicated(euc_data_all$Thornhill_bionomial_tip_match2) # all FALSE, no dup. values that can be problematic

euc_data_all = read.csv("outputs/CCADroughtdataMarch_2020_DN_AC_percentage_aff_trees_shifted.csv") %>% 
  dplyr::select(- per_aff_spavg) %>% 
  rename(per_aff_spavg = per_aff_spavg_shift) %>% 
  glimpse()


# this is the complete df with the additional column for the average percentage
# of drought impact per SPECIES


str(euc_tree) # the euc_tree$tip.label and euc_data_all$Thornhill_bionomial_tip_match2 seem to match well
# if some have terminal taxa have missing data they will be dropped with make.treedata()
is.binary(euc_tree)         # [1] FALSE
di_euc_tree = multi2di(euc_tree)
is.binary(di_euc_tree)      # [1] TRUE
is.rooted(di_euc_tree)      # [1] TRUE
is.ultrametric(di_euc_tree) # [1] TRUE
plot(di_euc_tree)




########################## Data Prep - Model and Assumption Tests ##########################

plot = ggplot(euc_data_all, aes(per_aff_spavg)) +
  geom_density()

# The data is highly over-dispersed: the variance is much larger than the mean
# fit a negative binomial model (over dispersion of data by number of zeros in data)
nbGLM = glm.nb(log(per_aff_spavg) ~ log(AICorrected_spavg), data = euc_data_all)
summary(nbGLM)

          # Call:
          #   glm.nb(formula = log(per_aff_spavg) ~ log(AICorrected_spavg), 
          #          data = euc_data_all, init.theta = 7190.172416, link = log)
          # 
          # Deviance Residuals: 
          #     Min       1Q   Median       3Q      Max  
          # -1.1814  -0.5341  -0.3999   0.8445   1.7425  
          # 
          # Coefficients:
          #                          Estimate Std. Error z value Pr(>|z|)    
          #              (Intercept)  -0.5389     0.1608  -3.352 0.000802 ***
          #   log(AICorrected_spavg)   1.0012     0.1619   6.185 6.21e-10 ***
          #   ---
          #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          # 
          # (Dispersion parameter for Negative Binomial(7190.172) family taken to be 1)
          # 
          #     Null deviance: 225.05  on 335  degrees of freedom
          # Residual deviance: 184.96  on 334  degrees of freedom
          # (6 observations deleted due to missingness)
          # AIC: 308.42
          # 
          # Number of Fisher Scoring iterations: 1
          # 
          # 
          #     Theta:  7190 
          # Std. Err.:  43668 
          # Warning while fitting theta: iteration limit reached 
          # 
          # 2 x log-likelihood:  -302.424
# AI has a coef of 1.00 which is statistically highly significant.
# It means that for each one-unit increase in AI, the expected log count of trees impacted by drought increases by 1.00.


## Check model assumption:

# NB model assume the conditional mean is not equal to the conditional variance.
# This inequality is captured by estimating a dispersion parameter that is held constant in a Poisson model.
# Thus, the Poisson model is actually nested in the negative binomial model.
# We can then use a likelihood ratio test to compare these two and test this model assumption.
# To do this, we will run our model as Poisson:

# compared to a Poisson model (assumes  mean ~ variance):
pGLM = glm(log(per_aff_spavg) ~ log(AICorrected_spavg), data = euc_data_all, family = poisson)
summary(pGLM)
          # Call:
          # glm(formula = per_aff_spavg ~ AICorrected_spavg, family = poisson, 
          #    data = euc_data_all)
          # 
          # Deviance Residuals: 
          #     Min       1Q   Median       3Q      Max  
          # -1.4595  -0.6427  -0.6087   0.3328   1.3331  
          #
          # Coefficients:
          #                     Estimate Std. Error z value Pr(>|z|)    
          #         (Intercept)  -1.8030     0.1489 -12.107  < 2e-16 ***
          #   AICorrected_spavg   1.2707     0.1846   6.885 5.77e-12 ***
          #   AIC: Inf
pchisq(2 * (logLik(nbGLM) - logLik(pGLM)), df = 1, lower.tail = FALSE)
# it seems that the negative binomial model is just as appropriate as the poisson model for estimating the dispersion parameter.
# compare test fits:
countreg::rootogram(pGLM)
countreg::rootogram(nbGLM)
# A bar hanging below 0 indicates underfitting.
# A bar hanging above 0 indicates overfitting.
# The counts have been transformed with a square root transformation to prevent smaller counts from getting obscured and overwhelmed by larger counts.
# We see a massive overfitting for counts 0 and 2; and slightly overfitting for all remaining counts.
# However the negative binomial model has a much lower AIC (308.42) which is not hard considering
# the poisson model has AIC: Inf...
plot = ggplot(euc_data_all, aes(height_max_m_spavg)) +
  geom_density()
plot = ggplot(euc_data_all, aes(LeafLenghtavg_mm_spavg)) +
  geom_density()

########################## Data Prep - Formatting for PGLS ##########################

# this is a species level analysis so we have to take the species level tree, not the UniqueID one
# here we are using "Thornhill_bionomial_tip_match2" to match to the tree

glimpse(euc_data_all) 
check = name.check(phy = di_euc_tree, data = euc_data_all,
                   data.names = euc_data_all$Thornhill_bionomial_tip_match2)
check # we have 391 tips with no data, this is because we don't have drought assessment for all species
      # and it is the crucial factor for which NAs have been removed
      # however only one observation in our data has no tip so the matching 
      # by Australian Taxonomic Criteria worked well
eucstuff = make.treedata(tree = di_euc_tree, data = euc_data_all,
                         name_column = "Thornhill_bionomial_tip_match2.x")
# name_column specifies col in data that matches the names in the tree

# check out the tree and data:
eucstuff$phy
plot(eucstuff$phy) # the tree looks good
glimpse(eucstuff$dat) # the data looks good as well

# Make a new column called tiplabel with the tip labels in it
eucstuff$dat$tiplabel = eucstuff$phy$tip.label
# Force eucstuff$dat to be a data frame
mydata = as.data.frame(eucstuff$dat)
summary(mydata) # only 5 rows have missing data for the environment variables because we don't have the
                # we don't have the seedlot location data
# And now we create the tree as a separate object
# the difference between mytree and di_euc_tree is that by passing the tree through the make.treedata()
# function, we dropped the tips that don't have data points in the matching dataframe
mytree = eucstuff$phy
 


########################## Data Prep - Are we missing non-matching data? ##########################

# Here the crucial part is to make sure that we have the minimum possible number of observations
# (averaged at the binomial data) that is not matching with the tree by comparing the data 
# transformed by make.treedata() with the dataset pre-transformation (euc_data_all).
glimpse(mydata)
summary(mydata)
# Now we are ready for the analysis!




########################## PCA - all traits ##########################

pPCA_euctree = di_euc_tree # custom tree for analysis
check = name.check(phy = pPCA_euctree, data = mydata, # Check whether the names match in the data and the tree
                   data.names = pPCA_euctree$tip.label)
check

## Phylogenetic pattern of quantitative traits with orthogram {adephylo}
# find code from p. 214 of Emmanuel Paradis book on APE with R

dat = phylo4d(pPCA_euctree, mydata, missing.data = "OK")
res = ppca(dat)
res = ppca(dat, scannf = T, nfposi=1, nfnega=1, method="Abouheif")
adephylo::orthogram(mydata, tree)




########################## Phylogenetic signal - Drought Impact ##########################

# We already prepared the data in the right format earlier, now we are just going to make a copy of it for
# the phylogenetic signal analysis.
# Make a new column called tiplabel with the tip labels in it
eucstuff2 = eucstuff
eucstuff2$dat$tiplabel = eucstuff2$phy$tip.label
# Force mydata to be a data frame
mydata2 = as.data.frame(eucstuff2$dat)

# Create DI containing just drought impact values
DI = pull(mydata2, per_aff_spavg)

# Notice that this is currently just a long list of numbers. We can then name these values with the species names from mydata using the function names. Note that this requires the trait data is in the same order as the tree tip labels, but luckily make.treedata does this automatically.

# Give DI names = species names at the tips of the phylogeny
names(DI) = mydata2$tiplabel
# Look at the first few rows
head(DI)

# Now we have a list of values with associated species names.

# Estimate Pagel's lambda
lambdaDI = phylosig(K_euctree, DI, method = "lambda", test = T)

#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaDI

      # Phylogenetic signal lambda :    0.719254 
      #               logL(lambda) :    -75.4051 
      #               LR(lambda=0) :     91.4199 
      # P-value (based on LR test) : 1.16201e-21

# The λ estimate for DI is around 0.719.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).

# Here P < 0.001. We interpret this as λ being significantly different from 0, i.e. there is significant phylogenetic signal in drought sensitivity.

# Now let's estimate Blomberg’s K. To do so we also use phylosig but with method = K.
# Estimate Blomberg’s *K*

K_Impact = phylosig(K_euctree, DI, method = "K", test = TRUE, nsim = 1000)
K_Impact
      #                  Phylogenetic signal K : 0.0247906 
      # P-value (based on 1000 randomizations) :     0.002




########################## Phylogenetic signal - Aridity Index ##########################

AI_data = pull(mydata2, AICorrected_spavg)
names(AI_data) = mydata2$tiplabel
head(AI_data)
lambdaAI = phylosig(K_euctree, AI_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaAI

        # Phylogenetic signal lambda :     0.69848 
        #               logL(lambda) :    -29.6423 
        #               LR(lambda=0) :     85.7544 
        # P-value (based on LR test) : 2.03726e-20 

# The λ estimate for DI is ~0.700.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in Aridity Index.

K_Impact = phylosig(K_euctree, AI_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
        #                  Phylogenetic signal K : 0.027016 
        # P-value (based on 1000 randomizations) :    0.001 




########################## Phylogenetic signal - Maximum Height ##########################

MH_data = pull(mydata2, height_max_m_spavg)
names(MH_data) = mydata2$tiplabel
head(MH_data)
lambdaMH = phylosig(K_euctree, MH_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMH
          # Phylogenetic signal lambda :    0.502976 
          #               logL(lambda) :    -993.049 
          #               LR(lambda=0) :     17.8726 
          # P-value (based on LR test) : 2.36195e-05 

# The λ estimate for DI is ~0.503.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MH_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
      #                  Phylogenetic signal K : 0.0166002 
      # P-value (based on 1000 randomizations) :     0.144 



########################## Phylogenetic signal - Mean Annual Precipitation ##########################

MAP_data = pull(mydata2, MAP_spavg)
names(MAP_data) = mydata2$tiplabel
head(MAP_data)
lambdaMAP = phylosig(K_euctree, MAP_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMAP
        # Phylogenetic signal lambda :    0.733893 
        #               logL(lambda) :    -1761.98 
        #               LR(lambda=0) :     92.7468 
        # P-value (based on LR test) : 5.94312e-22 

# The λ estimate for DI is ~0.734.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MAP_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
        #                  Phylogenetic signal K : 0.0304961 
        # P-value (based on 1000 randomizations) :     0.001 




########################## Phylogenetic signal - Mean Annual Temperature ##########################

MAT_data = pull(mydata2, MAT_spavg)
names(MAT_data) = mydata2$tiplabel
head(MAT_data)
lambdaMAT = phylosig(K_euctree, MAT_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMAT
          # Phylogenetic signal lambda :    0.732944 
          #               logL(lambda) :    -616.727 
          #               LR(lambda=0) :     85.2582 
          # P-value (based on LR test) : 2.61841e-20 
# The λ estimate for DI is ~0.734.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MAT_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0307412 
# P-value (based on 1000 randomizations) :     0.001 




########################## Phylogenetic signal - PDryM ##########################

PDM_data = pull(mydata2, PDryM_spavg)
names(PDM_data) = mydata2$tiplabel
head(PDM_data)
lambdaPDM = phylosig(K_euctree, PDM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaPDM
# Phylogenetic signal lambda :    0.786641 
#               logL(lambda) :    -1006.68 
#               LR(lambda=0) :     130.845 
# P-value (based on LR test) : 2.67657e-30 
# The λ estimate for DI is ~0.787.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, PDM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0319226 
# P-value (based on 1000 randomizations) :     0.001 




########################## Phylogenetic signal - MaxTWarmM ##########################

MTWM_data = pull(mydata2, MaxTWarmM_spavg)
names(MTWM_data) = mydata2$tiplabel
head(MTWM_data)
lambdaMTWM = phylosig(K_euctree, MTWM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMTWM
# Phylogenetic signal lambda :    0.775332 
#               logL(lambda) :    -658.284 
#               LR(lambda=0) :     96.3987 
# P-value (based on LR test) : 9.39288e-23 
# The λ estimate for DI is ~0.775.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MTWM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.026974 
# P-value (based on 1000 randomizations) :     0.001 




########################## Phylogenetic signal - LeafWidth_mean ##########################

LWM_data = pull(mydata2, LeafWidth_avg_mm_spavg)
names(LWM_data) = mydata2$tiplabel
head(LWM_data)
lambdaLWM = phylosig(K_euctree, LWM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWM
# Phylogenetic signal lambda :    0.413435 
#               logL(lambda) :    -881.044 
#               LR(lambda=0) :     24.0894 
# P-value (based on LR test) :  9.1963e-07 
# The λ estimate for DI is ~0.413.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.0160404 
# P-value (based on 1000 randomizations) :      0.098 




########################## Phylogenetic signal - LeafWidth_max ##########################

LWMax_data = pull(mydata2, leaf_width_mm_max)
names(LWMax_data) = mydata2$tiplabel
head(LWMax_data)
lambdaLWMax = phylosig(K_euctree, LWMax_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWMax
# Phylogenetic signal lambda :    0.651092 
#               logL(lambda) :    -1284.91 
#               LR(lambda=0) :      35.357 
# P-value (based on LR test) : 2.74479e-09 
# The λ estimate for DI is ~0.651.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWMax_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.0273574 
# P-value (based on 1000 randomizations) :      0.001 




########################## Phylogenetic signal - LeafWidth_min ##########################

LWMin_data = pull(mydata2, leaf_width_mm_min)
names(LWMin_data) = mydata2$tiplabel
head(LWMin_data)
lambdaLWMin = phylosig(K_euctree, LWMin_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWMin
# Phylogenetic signal lambda :    0.308525 
#               logL(lambda) :     -994.97 
#               LR(lambda=0) :     6.50056 
# P-value (based on LR test) :   0.0107841 
# The λ estimate for DI is ~0.309.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWMin_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :    0.02503 
# P-value (based on 1000 randomizations) :      0.004 




########################## PIC - Drought Impact ##########################

# We already prepared the data in the right format earlier, now we are just going to make a copy of it for
# the phylogenetic signal analysis.
# Make a new column called tiplabel with the tip labels in it
eucstuff2 = eucstuff
eucstuff2$dat$tiplabel = eucstuff2$phy$tip.label
# Force mydata to be a data frame
mydata2 = as.data.frame(eucstuff2$dat)

# Create DI containing just drought impact values
DI = pull(mydata2, per_aff_spavg)

# Notice that this is currently just a long list of numbers. We can then name these values with the species names from mydata using the function names. Note that this requires the trait data is in the same order as the tree tip labels, but luckily make.treedata does this automatically.

# Give DI names = species names at the tips of the phylogeny
names(DI) = mydata2$tiplabel
# Look at the first few rows
head(DI)

# Now we have a list of values with associated species names.

# Estimate Pagel's lambda
lambdaDI = phylosig(K_euctree, DI, method = "lambda", test = T)

#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaDI

# Phylogenetic signal lambda :    0.719254 
#               logL(lambda) :    -75.4051 
#               LR(lambda=0) :     91.4199 
# P-value (based on LR test) : 1.16201e-21

# The λ estimate for DI is around 0.719.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).

# Here P < 0.001. We interpret this as λ being significantly different from 0, i.e. there is significant phylogenetic signal in drought sensitivity.

# Now let's estimate Blomberg’s K. To do so we also use phylosig but with method = K.
# Estimate Blomberg’s *K*

K_Impact = phylosig(K_euctree, DI, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0247906 
# P-value (based on 1000 randomizations) :     0.002




########################## PIC - Aridity Index ##########################

AI_data = pull(mydata2, AICorrected_spavg)
names(AI_data) = mydata2$tiplabel
head(AI_data)
lambdaAI = phylosig(K_euctree, AI_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaAI

# Phylogenetic signal lambda :     0.69848 
#               logL(lambda) :    -29.6423 
#               LR(lambda=0) :     85.7544 
# P-value (based on LR test) : 2.03726e-20 

# The λ estimate for DI is ~0.700.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in Aridity Index.

K_Impact = phylosig(K_euctree, AI_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.027016 
# P-value (based on 1000 randomizations) :    0.001 




########################## PIC - Maximum Height ##########################

MH_data = pull(mydata2, height_max_m_spavg)
names(MH_data) = mydata2$tiplabel
head(MH_data)
lambdaMH = phylosig(K_euctree, MH_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMH
# Phylogenetic signal lambda :    0.502976 
#               logL(lambda) :    -993.049 
#               LR(lambda=0) :     17.8726 
# P-value (based on LR test) : 2.36195e-05 

# The λ estimate for DI is ~0.503.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MH_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0166002 
# P-value (based on 1000 randomizations) :     0.144 



########################## PIC - Mean Annual Precipitation ##########################

MAP_data = pull(mydata2, MAP_spavg)
names(MAP_data) = mydata2$tiplabel
head(MAP_data)
lambdaMAP = phylosig(K_euctree, MAP_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMAP
# Phylogenetic signal lambda :    0.733893 
#               logL(lambda) :    -1761.98 
#               LR(lambda=0) :     92.7468 
# P-value (based on LR test) : 5.94312e-22 

# The λ estimate for DI is ~0.734.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MAP_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0304961 
# P-value (based on 1000 randomizations) :     0.001 




########################## PIC - Mean Annual Temperature ##########################

MAT_data = pull(mydata2, MAT_spavg)
names(MAT_data) = mydata2$tiplabel
head(MAT_data)
lambdaMAT = phylosig(K_euctree, MAT_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMAT
# Phylogenetic signal lambda :    0.732944 
#               logL(lambda) :    -616.727 
#               LR(lambda=0) :     85.2582 
# P-value (based on LR test) : 2.61841e-20 
# The λ estimate for DI is ~0.734.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MAT_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0307412 
# P-value (based on 1000 randomizations) :     0.001 




########################## PIC - PDryM ##########################

PDM_data = pull(mydata2, PDryM_spavg)
names(PDM_data) = mydata2$tiplabel
head(PDM_data)
lambdaPDM = phylosig(K_euctree, PDM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaPDM
# Phylogenetic signal lambda :    0.786641 
#               logL(lambda) :    -1006.68 
#               LR(lambda=0) :     130.845 
# P-value (based on LR test) : 2.67657e-30 
# The λ estimate for DI is ~0.787.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, PDM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K : 0.0319226 
# P-value (based on 1000 randomizations) :     0.001 




########################## PIC - MaxTWarmM ##########################

MTWM_data = pull(mydata2, MaxTWarmM_spavg)
names(MTWM_data) = mydata2$tiplabel
head(MTWM_data)
lambdaMTWM = phylosig(K_euctree, MTWM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaMTWM
# Phylogenetic signal lambda :    0.775332 
#               logL(lambda) :    -658.284 
#               LR(lambda=0) :     96.3987 
# P-value (based on LR test) : 9.39288e-23 
# The λ estimate for DI is ~0.775.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, MTWM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.026974 
# P-value (based on 1000 randomizations) :     0.001 




########################## PIC - LeafWidth_mean ##########################

LWM_data = pull(mydata2, LeafWidth_avg_mm_spavg)
names(LWM_data) = mydata2$tiplabel
head(LWM_data)
lambdaLWM = phylosig(K_euctree, LWM_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWM
# Phylogenetic signal lambda :    0.413435 
#               logL(lambda) :    -881.044 
#               LR(lambda=0) :     24.0894 
# P-value (based on LR test) :  9.1963e-07 
# The λ estimate for DI is ~0.413.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWM_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.0160404 
# P-value (based on 1000 randomizations) :      0.098 




########################## PIC - LeafWidth_max ##########################

LWMax_data = pull(mydata2, leaf_width_mm_max)
names(LWMax_data) = mydata2$tiplabel
head(LWMax_data)
lambdaLWMax = phylosig(K_euctree, LWMax_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWMax
# Phylogenetic signal lambda :    0.651092 
#               logL(lambda) :    -1284.91 
#               LR(lambda=0) :      35.357 
# P-value (based on LR test) : 2.74479e-09 
# The λ estimate for DI is ~0.651.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWMax_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :  0.0273574 
# P-value (based on 1000 randomizations) :      0.001 




########################## PIC - LeafWidth_min ##########################

LWMin_data = pull(mydata2, leaf_width_mm_min)
names(LWMin_data) = mydata2$tiplabel
head(LWMin_data)
lambdaLWMin = phylosig(K_euctree, LWMin_data, method = "lambda", test = T)
#test = TRUE specifies that we want to run a likelihood ratio test to determine if λ is significantly different from 0. To look at the output we just type in the name of the model; lambdaDI in this case.
lambdaLWMin
# Phylogenetic signal lambda :    0.308525 
#               logL(lambda) :     -994.97 
#               LR(lambda=0) :     6.50056 
# P-value (based on LR test) :   0.0107841 
# The λ estimate for DI is ~0.309.
# logL is the log-likelihood, LR(lambda=0) is the log-likelihood for λ of 0, and P-value is the p value from a likelihood ratio test testing whether  λ is significantly different from 0 (no phylogenetic signal).
# P < 0.001 meaning λ significantly different from 0, i.e. there is significant phylogenetic signal in maximum height

K_Impact = phylosig(K_euctree, LWMin_data, method = "K", test = TRUE, nsim = 1000)
K_Impact
#                  Phylogenetic signal K :    0.02503 
# P-value (based on 1000 randomizations) :      0.004 




########################## PGLS model - DroughtImpact:Aridity Index ##########################

# Let's start by investigating the relationship between drought impact (percentage_aff_spavg) and AI.
# From Victoria we know that that they are highly correlated already. And we know from the phylo signal analysis
# that their values are both highly correlated with the phylogenetic tree.

# Plot drought impact against AI, coloured by habit, I wanted by section but it's not in it (yet).
ggplot(data = mydata, aes(x = log(per_aff_spavg),
                           y = log(AICorrected_spavg),
                          colour = section)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  stat_regline_equation(aes(label = ..eq.label..)) +
  stat_regline_equation(aes(label = ..rr.label..)) +  # R-squared value of 0.37 that is so exciting !!
  theme_bw()                                          # R-squared = 0.17 with per_aff_spavg
ggplot(mydata, aes(x = log(per_aff_spavg), 
                   y = log(AICorrected_spavg),
                   colour = height_max_m_spavg)) +
  geom_point() +
  theme_bw()

# To perform PGLS models in R, caper requires you to first combine the phylogeny and data into one object using the function comparative.data. This is similar to what we did with make.treedata, but it does some stuff that is particular to how caper works so we still need to do this here.
# 
# Note that vcv = TRUE stores a variance covariance matrix of your tree (you will need this for the pgls function).
# na.omit = FALSE stops the function from removing species without data for all variables.
# warn.dropped = TRUE will tell you if any species are not in both the tree and the data and are therefore dropped from the comparative data object.
# Here we won’t drop any species because we already did this using make.treedata.

# !! Here we want to use "mytree" and "mydata" that we just created !!


mytree$node.label=NULL # IMPORTANT: remove node labels otherwise comparative.data throws error

euc = comparative.data(phy = mytree, data = mydata,          # this part always surprises me when it
                       names.col = tiplabel, vcv = TRUE,     # actually works...
                       na.omit = FALSE, warn.dropped = TRUE)
str(euc)

# If you do need to drop species, this function will give a warning telling you that some species have been dropped.
#You can view the dropped species using:
euc$dropped$tips # character(0)
euc$dropped$unmatched.rows # character(0)

# Always make sure you check the list of dropped species is what you expected, it often reveals typos in your species names, or mismatches in taxonomies used etc. 


# The function for PGLS analyses in caper is pgls.
# To fit a model which uses the Maximum Likelihood (ML) estimate of λ we use the following code:
# Fit a PGLS model

model.pgls = pgls(log(per_aff_spavg) ~ log(AICorrected_spavg),
                  data = euc, lambda = "ML")
# !!! I keep having the following errors: !!!
# "Error in optim(optimPar, fn = pgls.likelihood, method = "L-BFGS-B", control = control,  : L-BFGS-B needs finite values of 'fn'"
# !! SOLUTION: !! the problem comes from the way optim() or other internal packages transform the original data
# the over representation of 0 in the dataset is problematic and results in 'fn' being Inf
# shifting all values by n+1 solves that issue

par(mfrow=c(2,2))
plot(model.pgls) # that is absolutely incredible, so the data is extremly phylogenetically constrained

# Let's look at the model outputs
anova(model.pgls)
    # Analysis of Variance Table
    # Sequential SS for pgls: lambda = 0.60, delta = 1.00, kappa = 1.00
    # 
    # Response: log(per_aff_spavg)
    #                         Df  Sum Sq  Mean Sq F value    Pr(>F)
    # log(AICorrected_spavg)   1 0.08133 0.081331  51.853 4.008e-12 ***
    #   Residuals            333 0.52230 0.001568
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# It’s always good to look at the output using anova first.
# This uses sequential sum of squares to tell you whether a model including your predictor variable(s) is a better fit than a model without your predictor variable(s).
# For a complex model with lots of predictors this is the easiest way to find out the answer to the question you were asking (this will become more obvious in the next example using two predictor variables).

# Here we asked “is there a significant effect of log(AICorrected_spavg) on log(per_aff_spavg)?”.
# Reporting this result in the manuscript:
# There was a significant effect of Aridity Index on the percentage of affected trees at CCA
# (PGLS: F = 51.853, df = 1, p < 0.001, λ = 0.60).


# We might also be interested in the model coefficients, i.e. the intercept and slope.
# To do this, just like we do for lm, we use summary:

# Let's look at the model coefficients:
summary(model.pgls)

      # Call:
      #   pgls(formula = log(per_aff_spavg) ~ log(AICorrected_spavg), data = euc, 
      #        lambda = "ML")
      # 
      # Residuals:
      #       Min        1Q    Median        3Q       Max 
      # -0.105300 -0.024399  0.004042  0.026206  0.130267 
      # 
      # Branch length transformations:
      #   
      # kappa [Fix] : 1.000
      # lambda [ ML] : 0.599
      # lower bound : 0.000, p = 4.1199e-07
      # upper bound : 1.000, p = < 2.22e-16
      # 95.0% CI : (0.341, 0.786)
      # delta [Fix] : 1.000
      # 
      # Coefficients:
      #                          Estimate Std. Error t value  Pr(>|t|)    
      #   (Intercept)            0.453424   0.097196  4.6651 4.475e-06 ***
      #   log(AICorrected_spavg) 0.151305   0.021012  7.2009 4.008e-12 ***
      #   ---
      #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      # 
      # Residual standard error: 0.0396 on 333 degrees of freedom
      # (5 observations deleted due to missingness)
      # Multiple R-squared: 0.1347,	Adjusted R-squared: 0.1321 
      # F-statistic: 51.85 on 1 and 333 DF,  p-value: 4.008e-12


# Reporting results in manuscript :
# There was a significant positive relationship between AI and drought impact
# (PGLS: slope ± SE = 0.151 ± 0.021, t = 7.2, df = 333, p < 0.001, λ = 0.599).

# Note that as well as the standard regression outputs, the summary output includes the estimated ML value of λ (0.976) and p values from likelihood ratio tests showing whether the ML λ is significantly different from 0 or 1.
# !! You may have also noticed κ and δ in the PGLS output.κ and δ are also tree transformations which can improve the fit of the data to the tree. It is possible to use pgls to optimise κ or δ (using kappa = “ML” or delta = “ML” instead of lambda = “ML” in the code above). We will not cover this here. Optimizing more than one of these parameters at the same time is not advisable because it would be impossible to interpret the results!


# Let's plot the results:
ggplot(mydata, aes(x = log(AICorrected_spavg), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls)[2], 
              intercept = coefficients(model.pgls)[1]) +
  theme_bw()
# Note that coefficients(model.pgls) gives us the intercept coefficients(model.pgls)[2], and slope coefficients(model.pgls)[2] of the line, allowing us to use geom_abline to fit the line.

# check for the model coefficients (intercept and slope)
coefficients(model.pgls)
      #      (Intercept) log(AICorrected_spavg)
      #        0.4534236              0.1513047






########################## PGLS model - DroughtImpact:MAT ##########################
model.pgls2 = pgls(log(per_aff_spavg) ~ log(MAT_spavg),
                  data = euc, lambda = "ML")
# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls2) # GORG'

# Model outputs:
anova(model.pgls2)
# There was a significant effect of MAT on the percentage of affected trees at CCA
# (PGLS: F = 5.272, df = 1, p < 0.05, λ = 0.72).

# Model coefficients:
summary(model.pgls2)
# There was a significant positive relationship between MAT and drought impact
# (PGLS: slope ± SE = -0.151 ± 0.066, t = -2.3, df = 333, p < 0.05, λ = 0.716).

# Plot the results:
ggplot(mydata, aes(x = log(MAT_spavg), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls2)[2], 
              intercept = coefficients(model.pgls2)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls2)
    #     (Intercept) log(MAT_spavg) 
    #       0.7350491     -0.1517047




########################## PGLS model - DroughtImpact:MAP ##########################

model.pgls3 = pgls(log(per_aff_spavg) ~ log(MAP_spavg),
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls3) 

# Model outputs:
anova(model.pgls3)
# There was a significant effect of MAP on the percentage of affected trees at CCA
# (PGLS: F = 43.761, df = 1, p < 0.001, λ = 0.63).

# Model coefficients:
summary(model.pgls3)
# There was a significant positive relationship between MAT and drought impact
# (PGLS: slope ± SE = 0.183 ± 0.027, t = 6.6, df = 333, p < 0.001, λ = 0.632).

# Plot the results:
ggplot(mydata, aes(x = log(MAP_spavg), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls3)[2], 
              intercept = coefficients(model.pgls3)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls3)
#     (Intercept) log(MAP_spavg) 
#      -0.8966740      0.1831992




########################## PGLS model - DroughtImpact:PDryM ##########################

model.pgls4 = pgls(log(per_aff_spavg) ~ PDryM_spavg,
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls4) 

# Model outputs:
anova(model.pgls4)
# There was a significant effect of PDryM on the percentage of affected trees at CCA
# (PGLS: F = 13.695, df = 1, p < 0.001, λ = 0.70).

# Model coefficients:
summary(model.pgls4)
# There was a significant positive relationship between PDryM and drought impact
# (PGLS: slope ± SE = 0.003 ± 0.001, t = 3.700, df = 333, p < 0.001, λ = 0.699).

# Plot the results:
ggplot(mydata, aes(x = PDryM_spavg, 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls4)[2], 
              intercept = coefficients(model.pgls4)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls4)
#     (Intercept)    PDryM_spavg 
#     0.213831386    0.003268861




########################## PGLS model - DroughtImpact:MaxTWarmM ##########################

model.pgls5 = pgls(log(per_aff_spavg) ~ MaxTWarmM_spavg,
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls5) 

# Model outputs:
anova(model.pgls5)
# There was a significant effect of PDryM on the percentage of affected trees at CCA
# (PGLS: F = 28.008, df = 1, p < 0.001, λ = 0.66).

# Model coefficients:
summary(model.pgls5)
# There was a significant positive relationship between PDryM and drought impact
# (PGLS: slope ± SE = -0.018 ± 0.003, t = -5.292, df = 333, p < 0.001, λ = 0.660).

# Plot the results:
ggplot(mydata, aes(x = MaxTWarmM_spavg, 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls5)[2], 
              intercept = coefficients(model.pgls5)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls5)
#     (Intercept) MaxTWarmM_spavg 
#      0.83418319     -0.01771587




########################## PGLS model - DroughtImpact:MaxHeight ##########################

model.pgls6 = pgls(log(per_aff_spavg) ~ log(height_max_m_spavg),
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls6) 

# Model outputs:
anova(model.pgls6)
# There was a significant effect of species maximum height on the percentage of affected trees at CCA
# (PGLS: F = 14.479, df = 1, p < 0.001, λ = 0.69).

# Model coefficients:
summary(model.pgls6)
# There was a significant positive relationship between PDryM and drought impact
# (PGLS: slope ± SE = 0.067 ± 0.018, t = 3.805, df = 333, p < 0.001, λ = 0.691).

# Plot the results:
ggplot(mydata, aes(x = log(height_max_m_spavg), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls6)[2], 
              intercept = coefficients(model.pgls6)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls6)
#     (Intercept) log(height_max_m_spavg) 
#       0.1209815               0.0674866




########################## PGLS model - DroughtImpact:LeafWidth_mean ##########################

model.pgls7 = pgls(log(per_aff_spavg) ~ log(LeafWidth_avg_mm_spavg),
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls7) 

# Model outputs:
anova(model.pgls7)
# There was a significant effect of species average leaf width on the percentage of affected trees at CCA
# (PGLS: F = 4.313, df = 1, p < 0.05, λ = 0.71).

# Model coefficients:
summary(model.pgls7)
# There was a significant positive relationship between PDryM and drought impact
# (PGLS: slope ± SE = 0.053 ± 0.025, t = 2.077, df = 333, p < 0.03, λ = 0.715).

# Plot the results:
ggplot(mydata, aes(x = log(LeafWidth_avg_mm_spavg), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls7)[2], 
              intercept = coefficients(model.pgls7)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls7)
#     (Intercept) log(LeafWidth_avg_mm_spavg) 
#       0.1209815               0.0674866



########################## PGLS model - DroughtImpact:LeafWidth_max ##########################

model.pgls8 = pgls(log(per_aff_spavg) ~ log(leaf_width_mm_max),
                   data = euc, lambda = "ML")

# Let's first check the model
par(mfrow=c(2,2))
plot(model.pgls8) 

# Model outputs:
anova(model.pgls8)
# There was a significant effect of species average leaf width on the percentage of affected trees at CCA
# (PGLS: F = 4.313, df = 1, p < 0.05, λ = 0.71).

# Model coefficients:
summary(model.pgls8)
# There was a significant positive relationship between PDryM and drought impact
# (PGLS: slope ± SE = 0.053 ± 0.025, t = 2.077, df = 333, p < 0.03, λ = 0.715).

# Plot the results:
ggplot(mydata, aes(x = log(leaf_width_mm_max), 
                   y = log(per_aff_spavg))) +
  geom_point() +
  geom_abline(slope = coefficients(model.pgls8)[2], 
              intercept = coefficients(model.pgls8)[1]) +
  theme_bw()
# Model coefficients (intercept and slope)
coefficients(model.pgls8)
#     (Intercept) log(LeafWidth_avg_mm_spavg) 
#       0.1209815               0.0674866




########################## ASR - DroughtImpact ##########################
########################## ASR - Aridity Index ##########################
########################## ASR - Maximum Height ##########################
########################## ASR - Resprouting strategy ##########################
########################## ASR - MAT ##########################
########################## ASR - MAP ##########################
