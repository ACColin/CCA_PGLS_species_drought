
librarian::shelf(ape, phytools, picante, ggtree, BioGeoBEARS, tidyverse, geiger)


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
e1 = tre_node_adjust$edge[, 1]                  # parent node
e2 = tre_node_adjust$edge[, 2]                  # child node
EL = tre_node_adjust$edge.length

# Also set up an ages variable that will hold internal calculations for how old a node should be.
ages = numeric(N + tre_node_adjust$Nnode)

# Start iterating:
for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) {          # If we haven’t already stored an age for the parent node,
                                    # go ahead and compute that now from the (left)1 child node and the current edge length.
    ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
  } else {                          # Otherwise, retrieve the stored age for the parent node,
                                    # and re-compute what the age should be based on the (right)1 child node.
    recorded_age <- ages[e1[ii]]
    new_age <- ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) {  # Now test whether those ages differ.
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



################ Fix 2 - non-negative least squares ################ >>>>>>>>>>>>>>>


# This is the default approach used in the R function phytools::force.ultrametric and the topic of a phytools blog post.

# For a given tree, this function can find the set of edge lengths with
# implied distances with minimum sum-of-squared differences to the true distances
# - in this case the patristic distances on our phylogeny.

tre_nnls = phangorn::nnls.tree(cophenetic(euctree), euctree, rooted = TRUE)
is.ultrametric(tre_nnls)

# Since it’s trying to minimize differences among the pairwise tip distance matrix,
# you’d expect many branches to be adjusted. Plotting the differences show that this is indeed the case:
diff_edge_lengths(euctree, tre_nnls)

################ <<<<<<<<<<<<<<<<


################ Fix 3 - node adjustment  ################ >>>>>>>>>>>>>>>

# This is the approach optionally used in DendroPy,4 and is how TACT fixes ultrametricity issues if asked.
# Whenever a node’s age differs between its left and right children,
# correct one of the branch lengths so the node’s age will be calculated the same regardless of
# whether you’re using the left descendants or the right descendants.1

# Change the prior loop that checks node ages to additionally adjust the branch length:

for (ii in seq_along(EL)) {
  if (ages[e1[ii]] == 0) {
    ages[e1[ii]] = ages[e2[ii]] + EL[ii]
  } else {
    recorded_age = ages[e1[ii]]
    new_age = ages[e2[ii]] + EL[ii]
    if (recorded_age != new_age) {
      cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
      + EL[ii] = recorded_age - ages[e2[ii]]
    }
  }
}

################## end tree setup #######################

############### LET'S GIVE IT A SHOT NOW ##############

thornhill_labels = read.csv("data/tree_labels.csv")
tree = read.tree("plots/final_tree_species_tips_ultra-extended.tre")
ggtree(tree,layout = "fan")
tree = multi2di(tree)
plotTree(tree, fsize=.6)
new_tree = tree
#new_tree$tip.label = thornhill_labels$name4[TreeTools::match(tree$tip.label, thornhill_labels$original)]
par(mfrow=c(1,2))
plot(tree)
plot(new_tree)
euctree = new_tree
is.rooted(euctree)
is.ultrametric(euctree)

eucdat = read.csv("outputs/all_eucs_PCA_dataframe_sp_avg.csv",header = T) %>% 
                    dplyr::select(binomial, height_max_m,) %>% 
                    na.omit()
eucdat

MH_sp_avg = aggregate(eucdat,by = list(binomial = eucdat$binomial), FUN = mean)
row.names(MH_sp_avg) = MH_sp_avg$binomial

d = eucdat %>%
  dplyr::select(binomial)
rownames(eucdat) = NULL
data = cbind(rownamesd, eucdat)
# I kept only one continuous and one discrete variable

eucdat_sp = eucdat %>%
  dplyr::select(binomial, Habit_num, height_max_m) %>%
  dplyr::distinct(binomial, .keep_all = T) %>%
  glimpse()
is.na(eucdat_sp) # we are good to go it seems..
eucdat_sp
######## Reconstructing ancestral state for discrete variables #########


## Estimating ancestral character states for discrete characters under a continuous-time Markov chain

# habit = eucdat_sp$Habit_num
# height = eucdat_sp$height_max_m
# DF.traits = data.frame(habit, height, row.names = row.names(eucdat_sp))
# eucdat = eucdat_sp[euctree$tip.label, ]

names(habit) = row.names(eucdat)

check = name.check(phy = test, data = MH_sp_avg,
                    data.names = MH_sp_avg$binomial)
check

# removing in tree what is not in data

list_names_1 = c("acmena_smithii", "agonis_flexuosa", "allosyncarpia_ternata", "amomyrtus_luma", "angophora_bakeri", "angophora_exul", "angophora_robur", "angophora_woodsiana", "arillastrum_gummiferum", "cloezia_floribunda", "corymbia_abbreviata", "corymbia_aparrerinja", "corymbia_arafurica", "corymbia_arenaria", "corymbia_aspera", "corymbia_blakei", "corymbia_bleeseri", "corymbia_cadophora", "corymbia_catenaria", "corymbia_clarksoniana", "corymbia_collina", "corymbia_confertiflora", "corymbia_dallachiana", "corymbia_dampieri", "corymbia_dendromerinx", "corymbia_dolichocarpa", "corymbia_dunlopiana", "corymbia_erythrophloia", "corymbia_ferruginea", "corymbia_flavescens", "corymbia_foelscheana", "corymbia_gilbertensis", "corymbia_grandifolia", "corymbia_greeniana", "corymbia_gummifera", "corymbia_henryi", "corymbia_hylandii", "corymbia_lamprophylla", "corymbia_latifolia", "corymbia_lenziana", "corymbia_nesophila", "corymbia_novaguineesis", "corymbia_opacula", "corymbia_papuana", "corymbia_pauciseta", "corymbia_petalophylla", "corymbia_plena", "corymbia_polysciada", "corymbia_porphyritica", "corymbia_porrecta", "corymbia_ptychocarpa", "corymbia_stockeri", "corymbia_tessellaris", "corymbia_torelliana", "corymbia_torta", "corymbia_zygophylla", "decaspermum_humile", "eucalyptopsis_alauda", "eucalyptopsis_papuana", "eucalyptus_abdita", "eucalyptus_acroleuca", "eucalyptus_agglomerata", "eucalyptus_alaticaulis", "eucalyptus_alba", "eucalyptus_argyphea", "eucalyptus_aspersa", "eucalyptus_balanites")
list_names_2 = c("eucalyptus_balanopelex", "eucalyptus_bensonii", "eucalyptus_brachyandra", "eucalyptus_brachycorys", "eucalyptus_brachyphylla", "eucalyptus_burdettiana", "eucalyptus_captiosa", "eucalyptus_carnabyi", "eucalyptus_castrensis", "eucalyptus_ceracea", "eucalyptus_chartaboma", "eucalyptus_communalis", "eucalyptus_conglomerata", "eucalyptus_crenulata", "eucalyptus_croajingolensis", "eucalyptus_cuspidata", "eucalyptus_cylindriflora", "eucalyptus_decurva", "eucalyptus_deglupta", "eucalyptus_depauperata", "eucalyptus_diptera", "eucalyptus_disclusa", "eucalyptus_dolorosa", "eucalyptus_dorotoxylon", "eucalyptus_dunnii", "eucalyptus_erectifolia", "eucalyptus_famelica", "eucalyptus_fastigata", "eucalyptus_fulgens", "eucalyptus_gigantangion", "eucalyptus_gregoriensis", "eucalyptus_houseana", "eucalyptus_ignorabilis", "eucalyptus_imlayensis", "eucalyptus_impensa", "eucalyptus_johnsoniana", "eucalyptus_kingsmillii", "eucalyptus_koolpinensis", "eucalyptus_kybeanensis", "eucalyptus_latens", "eucalyptus_lateritica", "eucalyptus_limitaris", "eucalyptus_lirata", "eucalyptus_medialis", "eucalyptus_microschema", "eucalyptus_miniata", "eucalyptus_minigwalica", "eucalyptus_mooreana", "eucalyptus_nerthicola", "eucalyptus_nitens", "eucalyptus_nobilis", "eucalyptus_oblonga", "eucalyptus_odontocarpa", "eucalyptus_pachyloma", "eucalyptus_paliformis", "eucalyptus_paludicola", "eucalyptus_patens", "eucalyptus_pellita", "eucalyptus_pendens", "eucalyptus_phoenicea")
list_names_3 = c("eucalyptus_pilligaensis", "eucalyptus_pilularis", "eucalyptus_placita", "eucalyptus_planchoniana", "eucalyptus_platyphylla", "eucalyptus_plauta", "eucalyptus_pluricaulis", "eucalyptus_prominens", "eucalyptus_pruiniramis", "eucalyptus_pulchella", "eucalyptus_quinniorum", "eucalyptus_ralla", "eucalyptus_rameliana", "eucalyptus_recurva", "eucalyptus_regnans", "eucalyptus_retinens", "eucalyptus_risdonii", "eucalyptus_saligna", "eucalyptus_semota", "eucalyptus_similis", "eucalyptus_splendens", "eucalyptus_suberea", "eucalyptus_surgens", "eucalyptus_talyuberlup", "eucalyptus_tectifica", "eucalyptus_tephrodes", "eucalyptus_tetrodonta", "eucalyptus_tintinnans", "eucalyptus_tortilis", "eucalyptus_transcontinentalis", "eucalyptus_trivalvis", "eucalyptus_umbra", "eucalyptus_urophylla", "eucalyptus_verrucata", "eucalyptus_wetarensis", "eucalyptus_wilcoxii", "eucalyptus_williamsiana", "eucalyptus_yangoura", "eucalyptus_yarraensis", "hypocalymma_linifolium", "kjellbergiodendron_celebicum", "kunzea_ericoides", "lophostemon_confertus", "lophostemon_suaveolens", "myrtus_communis", "stockwellia_quadrifida", "syncarpia_glomulifera", "syncarpia_hillii", "syzygium_angophoroides", "tepualia_stipularis")
pruned_euctree = drop.tip(euctree, list_names_1, trim.internal = T)
pruned_euctree = drop.tip(pruned_euctree, list_names_2, trim.internal = T)
pruned_euctree = drop.tip(pruned_euctree, list_names_3, trim.internal = T)

check = name.check(phy = pruned_euctree, data = MH_sp_avg,
                   data.names = MH_sp_avg$binomial)
check

#now adding tips for data that is not in tree
#update: can't add tips because they don't have another tip to be anchored to somehow.

is.rooted(pruned_euctree)
is.ultrametric(pruned_euctree)

tre_extend = pruned_euctree                           # have a copy of the tree
N = Ntip(tre_extend)
root_node = N + 1 
root_to_tip = dist.nodes(tre_extend)[1:N, root_node]

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

diff_edge_lengths(pruned_euctree, tre_extend) # using this fix basically increased the size of almost all the final
# branches lengths which were problematic especially in the very young clades
is.rooted(tre_extend)
is.ultrametric(tre_extend)

# Now we should save it and try to perform the ancestral state reconstruction and PGLS analyses with it.
write.tree(tre_extend, "plots/final_tree_species_tips_ultra-extended_for_MH_ACR.tre")

test = read.tree("plots/final_tree_species_tips_ultra-extended_for_MH_ACR.tre")
plot(test)
test_root = multi2di(test)
par(mfrow=c(1,2))
plot(test_root)
is.ultrametric(test_root)



ERreconstruction = ace(MH_sp_avg, test_root, type = "continuous", model = "ER")
ERreconstruction

list1 = c("angophora_leiocarpa", "corymbia clarksoniana", "corymbia_capricornia", "corymbia_chlorolampra", "corymbia_serendipita", "corymbia_setosa", "corymbia_umbonata", "eucalyptus_'arctata'_ms", "eucalyptus_'dorsiventralis'_ms", "eucalyptus_aenea", "eucalyptus_ammophila", "eucalyptus_angulosa_(Sa_variant)", "eucalyptus_annettae", "eucalyptus_arachneae", "eucalyptus_archeri", "eucalyptus_arenicola", "eucalyptus_aridomontana", "eucalyptus_armillata", "eucalyptus_aurifodina", "eucalyptus_baiophylla", "eucalyptus_bancroftii", "eucalyptus_banksii", "eucalyptus_baudiniana", "eucalyptus_bicostata", "eucalyptus_biterranea", "eucalyptus_brandiana", "eucalyptus_brownii", "eucalyptus_cajuputea", "eucalyptus_cambageana", "eucalyptus_camfieldii", "eucalyptus_capitanea", "eucalyptus_chapmaniana", "eucalyptus_chloroclada", "eucalyptus_clelandiorum", "eucalyptus_coccifera", "eucalyptus_codonocarpa", "eucalyptus_connexa", "eucalyptus_corynodes", "eucalyptus_cullenii", "eucalyptus_dendromorpha", "eucalyptus_doratoxylon", "eucalyptus_dorrienii", "eucalyptus_dorrigoensis", "eucalyptus_drepanophylla", "eucalyptus_ecdysiastes", "eucalyptus_ecostata", "eucalyptus_efflorescens", "eucalyptus_elegans", "eucalyptus_erosa", "eucalyptus_expressa", "eucalyptus_exserta", "eucalyptus_falciformis", "eucalyptus_flavida", "eucalyptus_fracta", "eucalyptus_glomericassis", "eucalyptus_grasbyi", "eucalyptus_gregoryensis", "eucalyptus_helidonica", "eucalyptus_hypostomatica", "eucalyptus_imitans", "eucalyptus_improcera", "eucalyptus_jensenii", "eucalyptus_johnstonii", "eucalyptus_kabiana", "eucalyptus_lactea", "eucalyptus_laevis", "eucalyptus_leucophylla", "eucalyptus_longissima", "eucalyptus_lunata", "eucalyptus_macta", "eucalyptus_maidenii", "eucalyptus_mediocris", "eucalyptus_microcodon", "eucalyptus_molyneuxii", "eucalyptus_nebulosa", "eucalyptus_notabilis", "eucalyptus_notactites", "eucalyptus_nubilis", "eucalyptus_obstans", "eucalyptus_odorata", "eucalyptus_omissa", "eucalyptus_ordiana", "eucalyptus_pallida", "eucalyptus_peninsularis", "eucalyptus_planipes")
list2 = c("eucalyptus_platypus", "eucalyptus_plumula", "eucalyptus_prominula", "eucalyptus_provecta", "eucalyptus_psammitica", "eucalyptus_pseudoglobulus", "eucalyptus_punctata", "eucalyptus_relicta", "eucalyptus_repullulans", "eucalyptus_revelata", "eucalyptus_robusta", "eucalyptus_rossii", "eucalyptus_rowleyi", "eucalyptus_rupestris", "eucalyptus_sabulosa", "eucalyptus_selachiana", "eucalyptus_semiglobosa", "eucalyptus_sheathiana", "eucalyptus_shirleyi", "eucalyptus_sinuensis", "eucalyptus_sinuosa", "eucalyptus_sp._cable_Haul_Road", "eucalyptus_sp._Dartmoor", "eucalyptus_sp._Dunbar_Road", "eucalyptus_sp._esperance", "eucalyptus_sp._Lake_Magenta", "eucalyptus_sp._Little_Sandy_Desert", "eucalyptus_sp._North_Balladonia", "eucalyptus_sp._Queen_Victoria_Spring", "eucalyptus_sp._South_Newdegate_lakes", "eucalyptus_sp._Southern_wheatbelt", "eucalyptus_spectatrix", "eucalyptus_sporadica", "eucalyptus_staigeriana", "eucalyptus_striaticalyx", "eucalyptus_subcaerulea", "eucalyptus_suffulgens", "eucalyptus_sweedmaniana", "eucalyptus_tardecidens", "eucalyptus_taurina", "eucalyptus_tephroclada", "eucalyptus_terrica", "eucalyptus_trivalva", "eucalyptus_umbrawarrensis", "eucalyptus_urnigera", "eucalyptus_varia", "eucalyptus_virella", "eucalyptus_virginea", "eucalyptus_volcanica", "eucalyptus_wimmerensis", "eucalyptus_woollsiana", "angophora_leiocarpa", "corymbia clarksoniana", "corymbia_capricornia", "corymbia_chlorolampra", "corymbia_serendipita", "corymbia_setosa", "corymbia_umbonata", "eucalyptus_'arctata'_ms", "eucalyptus_'dorsiventralis'_ms", "eucalyptus_aenea", "eucalyptus_ammophila", "eucalyptus_angulosa_(Sa_variant)", "eucalyptus_annettae", "eucalyptus_arachneae", "eucalyptus_archeri", "eucalyptus_arenicola", "eucalyptus_aridomontana", "eucalyptus_armillata", "eucalyptus_aurifodina", "eucalyptus_baiophylla", "eucalyptus_bancroftii", "eucalyptus_banksii", "eucalyptus_baudiniana", "eucalyptus_bicostata", "eucalyptus_biterranea", "eucalyptus_brandiana", "eucalyptus_brownii", "eucalyptus_cajuputea", "eucalyptus_cambageana", "eucalyptus_camfieldii", "eucalyptus_capitanea", "eucalyptus_chapmaniana", "eucalyptus_chloroclada", "eucalyptus_clelandiorum", "eucalyptus_coccifera", "eucalyptus_codonocarpa", "eucalyptus_connexa", "eucalyptus_corynodes", "eucalyptus_cullenii", "eucalyptus_dendromorpha", "eucalyptus_doratoxylon", "eucalyptus_dorrienii", "eucalyptus_dorrigoensis", "eucalyptus_drepanophylla", "eucalyptus_ecdysiastes", "eucalyptus_ecostata", "eucalyptus_efflorescens", "eucalyptus_elegans", "eucalyptus_erosa", "eucalyptus_expressa", "eucalyptus_exserta", "eucalyptus_falciformis", "eucalyptus_flavida", "eucalyptus_fracta", "eucalyptus_glomericassis", "eucalyptus_grasbyi")
list3 = c("eucalyptus_gregoryensis", "eucalyptus_helidonica", "eucalyptus_hypostomatica", "eucalyptus_imitans", "eucalyptus_improcera", "eucalyptus_jensenii", "eucalyptus_johnstonii", "eucalyptus_kabiana", "eucalyptus_lactea", "eucalyptus_laevis", "eucalyptus_leucophylla", "eucalyptus_longissima", "eucalyptus_lunata", "eucalyptus_macta", "eucalyptus_maidenii", "eucalyptus_mediocris", "eucalyptus_microcodon", "eucalyptus_molyneuxii", "eucalyptus_nebulosa", "eucalyptus_notabilis", "eucalyptus_notactites", "eucalyptus_nubilis", "eucalyptus_obstans", "eucalyptus_odorata", "eucalyptus_omissa", "eucalyptus_ordiana", "eucalyptus_pallida", "eucalyptus_peninsularis", "eucalyptus_planipes", "eucalyptus_platypus", "eucalyptus_plumula", "eucalyptus_prominula", "eucalyptus_provecta", "eucalyptus_psammitica", "eucalyptus_pseudoglobulus", "eucalyptus_punctata", "eucalyptus_relicta", "eucalyptus_repullulans", "eucalyptus_revelata", "eucalyptus_robusta", "eucalyptus_rossii", "eucalyptus_rowleyi", "eucalyptus_rupestris", "eucalyptus_sabulosa", "eucalyptus_selachiana", "eucalyptus_semiglobosa", "eucalyptus_sheathiana", "eucalyptus_shirleyi", "eucalyptus_sinuensis", "eucalyptus_sinuosa", "eucalyptus_sp._cable_Haul_Road", "eucalyptus_sp._Dartmoor", "eucalyptus_sp._Dunbar_Road", "eucalyptus_sp._esperance", "eucalyptus_sp._Lake_Magenta", "eucalyptus_sp._Little_Sandy_Desert", "eucalyptus_sp._North_Balladonia", "eucalyptus_sp._Queen_Victoria_Spring", "eucalyptus_sp._South_Newdegate_lakes", "eucalyptus_sp._Southern_wheatbelt", "eucalyptus_spectatrix", "eucalyptus_sporadica", "eucalyptus_staigeriana", "eucalyptus_striaticalyx", "eucalyptus_subcaerulea", "eucalyptus_suffulgens", "eucalyptus_sweedmaniana", "eucalyptus_tardecidens", "eucalyptus_taurina", "eucalyptus_tephroclada", "eucalyptus_terrica", "eucalyptus_trivalva", "eucalyptus_umbrawarrensis", "eucalyptus_urnigera", "eucalyptus_varia", "eucalyptus_virella", "eucalyptus_virginea", "eucalyptus_volcanica", "eucalyptus_wimmerensis", "eucalyptus_woollsiana")

MH_sp_avg_red = MH_sp_avg[ ! row.names(MH_sp_avg) %in% list1, ]
MH_sp_avg_red = MH_sp_avg_red[ ! row.names(MH_sp_avg_red) %in% list2, ]
MH_sp_avg_red = MH_sp_avg_red[ ! row.names(MH_sp_avg_red) %in% list3, ]

ERreconstruction = ace(MH_sp_avg_red, test_root, type = "continuous", model = "ER")
ERreconstruction

check = name.check(phy = test_root, data = MH_sp_avg_red,
                   data.names = MH_sp_avg_red$binomial)
check


SYMreconstruction = ace(habit, euctree, type = "discrete", model = "SYM")
ARDreconstruction = ace(habit, euctree, type = "discrete", model = "ARD")



### continuous trait with http://www.phytools.org/eqg/Exercise_5.2/
set.seed(1234)
library(phytools)
x = fastBM(test_root)
C = vcvPhylo(test_root)
pp = rep(1, test_root$Nnode)