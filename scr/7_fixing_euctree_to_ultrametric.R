
librarian::shelf(ape, phytools, picante, ggtree, BioGeoBEARS)


######### Tree setup for species level traits ###########

thornhill_labels = read.csv("data/tree_labels.csv")
tree = read.nexus("data/Eucalypts_ML2_dated_r8s.tre")
ggtree(tree,layout = "fan")
new_tree = tree
new_tree$tip.label = thornhill_labels$name5[TreeTools::match(tree$tip.label, thornhill_labels$original)]
par(mfrow=c(1,2))
plot(tree)
plot(new_tree)
euctree = new_tree
is.rooted(euctree)
is.ultrametric(euctree) # we have a problem. according to ape, the tree is not ultrametric
# let's check in more details ...

## Following this to fix the tree and make it ultrametric:
## https://www.r-bloggers.com/2021/07/three-ways-to-check-and-fix-ultrametric-phylogenies/



###################################################
################ FINDING THE ISSUE ################
###################################################


################ Method 1 - check for variance: ################ >>>>>>>>>>


N = Ntip(euctree)
root_node = N + 1 
root_to_tip = dist.nodes(euctree)[1:N, root_node] # compute root to tip distances for all tips 
var(root_to_tip)                                  # compute variance using those distances
                                                  # [1] 9.570864e-13
sqrt(.Machine$double.eps) # compare to default tolerance in ape (sqrt of R's machine epsilon)
                          # to be considered ultrametric,
                          # the variance has to be smaller than the tolerance in ape
                          # [1] 1.490116e-08
# 1e-08 > 9e-13 therefore,
# the euc tree is ultrametric according to the variance method

################# <<<<<<<<<<<<<



################ Method 2 - check for relative difference: ################ >>>>>>>>>>


# Now let's check the relative difference of the minimum and maximum root-to-tip distances.
# This method has been implemented in is.ultrametric() ape v.4.0
# and therefore affects the result of the boolean function

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

################# <<<<<<<<<<<<<



################ Method 3 - check for node ages: ################ >>>>>>>>>>


# The previous two methods all relied on comparing some aspect of root-to-tip distances.
# Here, we’ll actually compare the distances to the tips of all nodes.
# This is the method used in DendroPy and the original impetus for this investigation.
# I’ll reimplement enough of this method in R to illustrate this technique.

# First, reorder the tree for postorder traversal, and set up some convenience variables.
tre_node_adjust = reorder(euctree, "postorder")
e1 = tre_node_adjust$edge[, 1] # parent node
e2 = tre_node_adjust$edge[, 2] # child node
EL = tre_node_adjust$edge.length

# Also set up an ages variable that will hold internal calculations for how old a node should be.
ages = numeric(N + tre_node_adjust$Nnode)

# Start iterating:
for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) { # If we haven’t already stored an age for the parent node,
                           # go ahead and compute that now from the (left)1 child node and the current edge length.
    ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
  } else { # Otherwise, retrieve the stored age for the parent node,
           # and re-compute what the age should be based on the (right)1 child node.
    recorded_age <- ages[e1[ii]]
    new_age <- ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) { # Now test whether those ages differ.
                                   # I could actually use either the variance or the relative difference method,
                                   # but here I’ll just check for absolute difference.
      cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
    }
  }
}

# Many internal nodes have ages that differ depending on
# whether you use the left or right child node to compute the age.
# This metric allows you to pinpoint exactly where in your phylogeny
# precision issues could be causing ultrametricity issues.
# I can imagine this being quite useful when you are grafting subtrees
# onto a backbone phylogeny and trying to figure out if you did the math on your branch lengths correctly.

################# <<<<<<<<<<





######################################################
################## FIXING THE ISSUE ##################
######################################################


################ Fix 1 - extending the tips ################ >>>>>>>>>>>>>>>

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


diff_edge_lengths <- function(phy, phy2) { # function to compare two phylogenies with identical topologies
                                           # but differing branch lengths
  diffs <- phy2$edge.length - phy$edge.length
  cols <- sign(diffs)
  cols[cols == 1] <- "#7fbc41"
  cols[cols == -1] <- "#de77ae"
  cols[cols == 0] <- NA
  plot(phy, show.tip.label = FALSE, no.margin = TRUE)
  edgelabels(pch = 15, col = cols)
  sprintf("%i longer branches, %i shorter branches", sum(diffs > 0), sum(diffs < 0))
}

diff_edge_lengths(euctree, tre_extend) # using this fix basically increased the size of almost all the final
                                       # branches lengths which were problematic especially in the very young clades
is.rooted(tre_extend)
is.ultrametric(tre_extend)

# Now we should save it and try to perform the ancestral state reconstruction and PGLS analyses with it.
write.tree(tre_extend, "plots/final_tree_species_tips_ultra-extended.tre")

################ <<<<<<<<<<<<<<<<

################## end tree setup #######################


############### LET'S GIVE IT A SHOT NOW ##############

thornhill_labels = read.csv("data/tree_labels.csv")
tree = read.tree("plots/final_tree_species_tips_ultra-extended.tre")
ggtree(tree,layout = "fan")
tree = multi2di(tree)
new_tree = tree
new_tree$tip.label = thornhill_labels$name5[TreeTools::match(tree$tip.label, thornhill_labels$original)]
par(mfrow=c(1,2))
plot(tree)
plot(new_tree)
euctree = new_tree
is.rooted(euctree)
is.ultrametric(euctree)

eucdat = read.csv("outputs/all_eucs_PCA_dataframe.csv", row.names = NULL) #%>% 
                    dplyr::select(Taxonsimplified_DN, UniqueID, Drought2020Scale.mean,AP,AICorrected,PDryQ,
                                  PrecDeficit,height_max_m,Habit_num,
                                  RegenStrategysimplified_num,LeafAreaEst_cm2,Aridclass_num, Section,Genus) %>% 
                    na.omit()
eucdat_sp = eucdat %>%
  dplyr::select(Taxonsimplified_DN, RegenStrategysimplified_num,
                Habit_num, Aridclass_num,
                height_max_m, Section, Genus) %>%
  dplyr::distinct(Taxonsimplified_DN, .keep_all = T)

eucdat = eucdat[euctree$tip.label, ]

######## Reconstructing ancestral state for discrete variables #########


## Estimating ancestral character states for discrete characters under a continuous-time Markov chain


habit = eucdat$Habit_num
names(habit) = row.names(eucdat)

ERreconstruction = ace(habit, euctree, type = "discrete", model = "ER")
ERreconstruction
SYMreconstruction = ace(habit, euctree, type = "discrete", model = "SYM")
ARDreconstruction = ace(habit, euctree, type = "discrete", model = "ARD")

