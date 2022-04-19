library(ape)
library(tidyverse)

full_tree <- read.tree("tree_inprogress_check_taxanames_in_double.tre")

# 
# tibble_tree <- as_tibble(full_tree) %>% # putting the full_tree$tips and full_tree$branchlength to the right format to merge
#   as_data_frame() %>%
#   glimpse()
# 
# branch_length_matrix <- left_join(new_tips_data, tibble_tree, by = c("phylo.tip.label.match" = "label")) %>% # merging
#   as_data_frame() %>%
#   select("new.tip_unique.ID", "node") # keeping only the important stuff
# 
# branch_length_matrix %>% # as matrix otherwise the function doesn't work
#   as.matrix() %>%
#   glimpse()

branch_length_matrix <- read_csv("bind.tip_tip_matrix.csv") %>% 
  as.matrix() %>%
  glimpse()

test_tree <- lapply(full_tree, bind.tip,
                    full_tree, length(branch_length_matrix[,1]),
                    length(branch_length_matrix[,2]))

test_tree <- bind.tip(full_tree, branch_length_matrix[,1],edge.length= , where=branch_length_matrix[,2])
plot(test_tree)

test_tree <- bind.tip(full_tree, tip.label = "angophora_costata_subsp._costata_2073", where=10)
write.tree(test_tree, "angoph_costata_test.tree")
plot(test_tree)

                                #####################
#####
# TRYING TO ADD ALL THE NAMES WITH THE FOR LOOP
#####


all_tree <- pruned_tree
new_tips_data <- read.csv("data/NewTips(onlyduplicates).csv") # only duplicates = only contains the tips from the current tree that will be receive new tips. not the tips that will remain unchanged (which I manually removed from the dataframe).

tibble_tree <- as_tibble(all_tree) %>% # putting the full_tree$tips and
  as_data_frame() %>%                   # full_tree$branchlength to the right format to merge
  glimpse()

branch_length_matrix <- left_join(new_tips_data, tibble_tree, by = c("phylo.tip.label.match" = "label")) %>% # merging
  as_data_frame() %>%
  select("uniqueID", "node") %>% # keeping only the important stuff
  as.matrix() %>% 
  glimpse()

#write_csv(branch_length_matrix, "bind.tip_tip_matrix.csv") # saved it just in case


head(all_tree)
head(branch_length_matrix)
?bind.tip

bldf <- as.data.frame(branch_length_matrix)
bldf$node <- as.numeric(bldf[,2])
head(bldf)


test_tree2 <- all_tree

for (i in 1:length(bldf$uniqueID)) { # the for loop takes a minute to run
  labs <- bldf$uniqueID[i]
  nodes <- which(test_tree2$tip.label == all_tree$tip.label[bldf$node[i]])
  test_tree2 <- bind.tip(test_tree2, tip.label = labs, edge.length= 0, where=nodes)
}

plot(test_tree2, cex = 0.25)

write.tree(test_tree2, "plots/with_new_tips_as_unique_ID_full_tree.tre")

                                  ###################################


head(full_tree)
head(branch_length_matrix)
?bind.tip

bldf <- as.data.frame(branch_length_matrix)
bldf$node <- as.numeric(bldf[,2])
head(bldf)


test_tree2 <- full_tree

for (i in 1:length(bldf$new.tip_unique.ID)) { # the for loop takes a minute to run
  labs <- bldf$new.tip_unique.ID[i]
  nodes <- which(test_tree2$tip.label == full_tree$tip.label[bldf$node[i]])
  test_tree2 <- bind.tip(test_tree2, tip.label = labs, edge.length= 0, where=nodes)
}

plot(test_tree2, cex = 0.25)

write.tree(test_tree2, "plots/with_new_tips_full_tree.tre")



