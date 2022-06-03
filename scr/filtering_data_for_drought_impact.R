shelf(ape, geiger, nlme, tidyverse, treeplyr, phytools, caper, ggtree, sqldf)
# 
# all <- read.csv("data/Phylo_AllCCAtraits.csv", header = T)
# no_NA <- read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv", header = T)
# 
# all_df <- all[,c("UniqueID","Drought2020Scale.mean")] %>% 
#   arrange(all_df,"UniqueID")
# noNA_df <- no_NA[,c("UniqueID","Drought2020Scale.mean")] %>% 
#   arrange(noNA_df,"UniqueID")
# 
# diff1 <- sqldf('SELECT * FROM all_df EXCEPT SELECT * FROM noNA_df') #rows from all_df which are not in noNA_df
# is.na(!diff1$Drought2020Scale.mean) 
# # it should all be TRUE meaning NA
# # this is meant to be only NA values but I see a few with a numerical value
# # ok so we have a bunch of cleaning to do.
# 
# diff2 <- sqldf('SELECT * FROM noNA_df EXCEPT SELECT * FROM all_df') #rows from noNA_df which are not in all_df
# # there shouldn't be any value assigned to diff2, let's see what they are.
# # I *think* they are unique IDs with an assigned numerical value *after* removing one rep with NA.
# 
# all_df <- all[,c("UniqueID","Drought2020Scale.mean")] 
# all_df$UniqueID <-  as.logical(all_df$UniqueID)
# all_df$Drought2020Scale.mean <-  as.logical(all_df$Drought2020Scale.mean)
# 
# glimpse(all)
# all$UniqueID <- as.logical(all$UniqueID)
# all$Drought2006Scale.mean <- as.logical(all$UniqueID)
# 
# all_df_reduced <- all_df %>% 
#   mutate(is.character(as.logical())) %>% 
#   glimpse()
# #  group_by(UniqueID) %>% 
# all_df_reduced <- all_df %>% 
#   filter(all_df, !duplicated(UniqueID) & is.na(Drought2020Scale.mean))



##### None of the above worked #####
#     New approach here is to use the summarized data, and any missing UniqueID will receive a NA value for Drought Impact because if absent it means no numerical value of DI so meaning it was previously dead or not usable.
#####

# First, load the data
# This data is a modified version of `DroughtImpact_Summary` in which I manually filtered the duplicates of UniqueID to remove some that had both NAs and a numerical value and all is left to do is summarize by Unique ID to pool together the NAs for the remaining duplicates of UniqueIDs across the year blocks.

data <- read.csv("data/DroughtImpact_SummaryBy_Completetraits_for_Blombergs_analysis.csv", header = T)

# Now let's get rid of UniqueID duplicates (they will have NAs in different year categories) so I have only one single UniqueID with one single value of Drought Impact (hopefully)

nodup_data <- data %>%                               # Summary by group using dplyr
  group_by(UniqueID) %>% 
  summarize(min = mean(Drought2006Scale.mean))

# before summarizing I had 1465 rows and now 1424 so I believe it worked, let's see
duplicated(nodup_data$UniqueID) # no duplicates, cool.

# now let's bring in the full dataset to have all the UniqueIDs, and summarize by UniqueID because there are so many duplicates in there

other_data <- read.csv("data/Phylo_AllCCAtraits.csv", header = T)

nodup_other_data <- other_data
nodup_other_data <- nodup_other_data[,c("UniqueID")]
  distinct(UniqueID, .keep_all = TRUE)  # I see 1919 rows which is more than what is at CCA, let's have a look
glimpse(nodup_other_data)
glimpse(nodup_data)

duplicated(nodup_other_data$UniqueID)

# Let's see what is in this that is absent from the other
ID_data <- nodup_data[, 1, drop = F]
ID_other_data <- nodup_other_data[, 1, drop = F]
ID_other_data <- nodup_other_data %>% 
  select(UniqueID)

diff1 <- sqldf('SELECT * FROM ID_other_data EXCEPT SELECT * FROM ID_data') #rows from nodup_other_data which are not in nodup_data

all_UniqueIDs_join = full_join(nodup_data, nodup_other_data, by = "UniqueID") %>% 
  filter(!duplicated("UniqueID"))

# Ok now let's get the UniqueIDs from the tree:

tree = read.tree("plots/final_tree_uniqueID.tre")
tree_UniqueIDs = as_tibble(tree$tip.label)

# now let's compare which ones are not in the tree and which ones are not in the data and reduce the dataset to the ones in the tree

in_tree_not_in_data <- sqldf('SELECT * FROM tree_UniqueIDs EXCEPT SELECT * FROM all_ID')
in_data_not_in_tree <- sqldf('SELECT * FROM all_ID EXCEPT SELECT * FROM tree_UniqueIDs')

final_DI_data = all_UniqueIDs_join %>% 
  filter(UniqueID %in% tree_UniqueIDs$value)

final_ID = final_DI_data[,1,drop=F]
in_final_data_not_in_tree <- sqldf('SELECT * FROM final_ID EXCEPT SELECT * FROM tree_UniqueIDs') # perfect
in_tree_not_in_final_data <- sqldf('SELECT * FROM tree_UniqueIDs EXCEPT SELECT * FROM final_ID') # Egilii1 only

final_DI_data = final_DI_data %>%
  distinct(UniqueID, .keep_all=T)

write.csv(final_DI_data, "data/DroughtImpactSummary_data_for_Blombergs_final.csv") # ok let's try this

