#' @title: 6: Principal Component Analysis with species life-history traits and environmental variables associated with drought
#' @author: A.C. Colin
#' 
#' # Analysis description:
#' 
#'
#'PACKAGES SETUP:
#'
librarian::shelf(ape, geiger, nlme, tidyverse, treeplyr, phytools, caper, tidytree, readxl, ggfortify, FactoMineR)
#'
#'
#' 
#'#### DATA SETUP:
#' >>>>>>>>>>> Data prep (don't run again)
eucs = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T) #%>% 
  dplyr::select(UniqueID, Drought2020Scale.mean, AICorrected, AP, PDryQ, PrecDeficit, height_max_m,
                MAT, PWetM, LeafAreaEst_cm2)
str(eucs)
head(eucs)
#'
#' Let's first add the classification so we can color samples by lineages:
#' This is Kevin's clean taxonomic metadata for CCA:
meta_taxo = read_xlsx("../../../Eucalyptus_CCA_PhyloGWAS/Data/kevin_murray_data/DNTaxonomyCleaned.xlsx", skip = 1, col_names = T)
str(meta_taxo)
#' Now we need to link back with our data using the CCA voucher metadata:=
meta_voucher = read_xlsx("../../../Eucalyptus_CCA_PhyloGWAS/Data/kevin_murray_data/CCATreesCleaned.xlsx", skip = 1, col_names = T)
str(meta_voucher)
#'
#' Now let's merge the two metas:
full_meta = left_join(meta_voucher, meta_taxo, by = c("CurrentName" = "Binomial"))
str(full_meta)
#' there is a bit of missing data for the taxonomic ranking but nothing we can't deal with
#'
#'Let's add this metadata to the dataset of traits and environmental variables:
all_data = left_join(eucs, full_meta, by = c("UniqueID" = "DNNumber"))
str(all_data)
#' good good
#' Now let's go for a crude first PCA to see what it looks like:
all_data = all_data %>% 
  dplyr::select(Taxonsimplified_DN, UniqueID, Drought2020Scale.mean, AICorrected, AP, PDryQ, PrecDeficit, height_max_m, RowTree, UniqueID, Genus, Subgenus, Section, Series, Subseries, Species, Subspecies, Wooddensity.g.cm3..mean, Diamthickness.mm..mean, SLA_m2.kg.mean, d13C.mean, totalN.mean, Habit, RegenStrategysimplified, LeafAreaEst_cm2, Aridclass, latitude, longitude)
str(all_data)
#'
#'Now we need to change the discrete categories of Habit, AridityClass and RegenStrategy in numerical categories:
# write.csv(all_data, "outputs/all_eucs_PCA_dataframe.csv")
#' <<<<<<<<<<<<<<<<<< Data prep (don't run again)
#' 
#' 
#' Load data here:
all_dat = read.csv("outputs/all_eucs_PCA_dataframe.csv", header = T) %>% 
  dplyr::select(UniqueID, Drought2020Scale.mean,AP,AICorrected,PDryQ,
                PrecDeficit,height_max_m,Habit_num,
                RegenStrategysimplified_num,LeafAreaEst_cm2,Aridclass_num, Section,Genus) %>% 
  na.omit()
str(all_dat)
#'
#'
#'
#'#### PRINCIPAL COMPONENT ANALYSIS:
#'
eucs.pca = prcomp(all_dat[,c("Drought2020Scale.mean","AP","AICorrected","PDryQ","PrecDeficit","height_max_m","Habit_num","RegenStrategysimplified_num","LeafAreaEst_cm2","Aridclass_num")], center = TRUE, scale. = TRUE)
#'
#'Summary of the prcomp object
summary(eucs.pca)
# Importance of components:
#                        PC1    PC2     PC3     PC4
# Standard deviation     2.3768 1.1429 0.90579 0.83264
# Proportion of Variance 0.5649 0.1306 0.08205 0.06933
# Cumulative Proportion  0.5649 0.6955 0.77758 0.84691
#                        PC5     PC6     PC7    PC8
# Standard deviation     0.82695 0.61793 0.56190 0.3391
# Proportion of Variance 0.06838 0.03818 0.03157 0.0115
# Cumulative Proportion  0.91529 0.95348 0.98505 0.9966
#                        PC9      PC10
# Standard deviation     0.18575 6.162e-05
# Proportion of Variance 0.00345 0.000e+00
# Cumulative Proportion  1.00000 1.000e+00
#'
#' PC1 explains 56% of the total variance i.e. more than half of the information in the data set can be encapsulated by just that one PC. PC2 explains 13% of the total variance.
#' 
str(eucs.pca)
#' Using auto lot to start with:
eucs.pca.plot = autoplot(eucs.pca, data = all_dat, colour = "Section") #, geom_point(aes(shape=Genus))) +
  scale_shape_manual(values=c(3, 4, 16, 17))
eucs.pca.plot
#' There is definitely some clustering according to the taxonomy
#' But there are also some interesting clusters.
#' The ugly mustard represents
biplot.eucs.pca = biplot(eucs.pca) #+
  stat_ellipse(aes(group=all_dat$Section))
biplot.eucs.pca
ggbiplot(eucs.pca)


#### PCA WITH FACTOMINER:
#' Compute PCA:
euc.pca = PCA(all_dat, quanti.sup = c(1:6,9), quali.sup = c(7,8,10), graph = F)
#' 
#' Extract PC1 and 2 for plotting
PC1 = euc.pca$ind$coord[,1]
PC2 = euc.pca$ind$coord[,2]
labs = rownames(euc.pca$ind$coord)
PCs = data.frame(cbind(PC1,PC2))
rownames(PC1) = labs
