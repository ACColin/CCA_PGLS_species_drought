#' @title: 6: Principal Component Analysis with species life-history traits and environmental variables associated with drought
#' @author: A.C. Colin
#' 
#' # Analysis description:
#' 
#'
#'Packages setup: 
librarian::shelf(ape, geiger, nlme, tidyverse, treeplyr, phytools, caper, tidytree, readxl)
#'
#'
#' 
#'Data setup:
eucs = read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T) #%>% 
  dplyr::select(UniqueID, Drought2020Scale.mean, AICorrected, AP, PDryQ, PrecDeficit, height_max_m,
                MAT, PWetM, LeafAreaEst_cm2)
str(eucs)
  head(eucs)
  
  #' Let's first add the classification so we can color samples by lineages:
#' This is Kevin's clean taxonomic metadata for CCA:
meta_taxo = read_xlsx("../../../Eucalyptus_CCA_PhyloGWAS/Data/kevin_murray_data/DNTaxonomyCleaned.xlsx", skip = 1, col_names = T)
str(meta_taxo)
#' Now we need to link back with our data using the CCA voucher metadata:
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
  dplyr::select(UniqueID, Drought2020Scale.mean, AICorrected, AP, PDryQ, PrecDeficit, height_max_m, RowTree, UniqueID, Genus, Subgenus, Section, Series, Subseries, Species, Subspecies, Wooddensity.g.cm3..mean, Diamthickness.mm..mean, SLA_m2.kg.mean, d13C.mean, totalN.mean, Habit, RegenStrategysimplified, LeafAreaEst_cm2, Aridclass, latitude, longitude)
str(all_data)
#'
#'Now we need to change the discrete categories of Habit, AridityClass and RegenStrategy in numerical categories:

#'
eucs.pca = prcomp(all_data[,c()], center = TRUE, scale. = TRUE)
#'
#'Summary of the prcomp object
summary(iris.pca)