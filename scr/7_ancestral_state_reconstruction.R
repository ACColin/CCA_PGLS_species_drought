

librarian::shelf(ape, phytools, picante)


######### Tree setup for species level traits ###########


thornhill_labels <- read.csv("data/tree_labels.csv")
tree = read.nexus("data/Eucalypts_ML2_dated_r8s.tre")
ggtree(tree,layout = "fan")
new_tree <- tree
new_tree$tip.label <- thornhill_labels$name5[TreeTools::match(tree$tip.label, thornhill_labels$original)]
par(mfrow=c(1,2))
plot(tree)
plot(new_tree)
is.rooted(euctree)
is.ultrametric(euctree) 

print(new_tree$edge.length)
n = length(new_tree$edge.length)
new_tree$edge.length = new_tree$edge.length[1:n-1] + .0000001
print(new_tree$edge.length)
plotTree(new_tree)

tree_ultra = chronos(new_tree, lambda = 0)

################## end tree setup #######################


euctree = read.tree("plots/final_tree_uniqueID.tre")
euctree = new_tree
eucdat = read.csv("outputs/all_eucs_PCA_dataframe.csv", header = T) %>% 
                    dplyr::select(UniqueID, Drought2020Scale.mean,AP,AICorrected,PDryQ,
                                  PrecDeficit,height_max_m,Habit_num,
                                  RegenStrategysimplified_num,LeafAreaEst_cm2,Aridclass_num, Section,Genus) %>% 
                    na.omit()
eucdat = eucdat[euctree$tip.label, ]

######## Reconstructing ancestral state for discrete variables #########


habit = eucdat$Habit_num
names(habit) = row.names(eucdat)

ERreconstruction = ace(habit, euctree, type = "discrete", model = "ER")
SYMreconstruction = ace(habit, euctree, type = "discrete", model = "SYM")
ARDreconstruction = ace(habit, euctree, type = "discrete", model = "ARD")

