librarian::shelf(caper, ape, geiger, phytools)

euc.tree = read.tree("plots/final_tree_CCAUniqueID_tips_dichotomous_ultra-extended.tre")
plot(euc.tree)
euc1 = euc.tree[[1]]
euc.traits = read.csv("data/UniqueID_tree_tips_data_for_ASR.csv", header = T)
glimpse(euc.traits)
pruned.tree = drop.tip(euc.tree,euc.tree$tip.label[-match(euc.traits$CCAtiplabel)])

traits.tree=comparative.data(phy=pruned.tree, data=eye.traits, names.col=Binomial, vcv=T, na.omit=F)

#With Primates
combined.mod1=pgls(ActivityCycle~residualAD, data=traits.tree, lambda = "ML")
summary(combined.mod1)

combined.mod2=pgls(ActivityCycle~residualVA, data=traits.tree, lambda = "ML")
summary(combined.mod2) #significant! 

combined.mod3=pgls(TrophicLevel~residualAD, data=traits.tree, lambda = "ML")
summary(combined.mod3)

combined.mod4=pgls(TrophicLevel~residualVA, data=traits.tree, lambda = "ML")
summary(combined.mod4)

#Without Primates
combined.prim1=pgls(xActivityCycle~xresidAD, data=traits.tree, lambda = "ML")
summary(combined.prim1)

combined.prim2=pgls(xActivityCycle~xresidVA, data=traits.tree, lambda = "ML")
summary(combined.prim2) #still significant after removing primates

combined.prim3=pgls(xTrophicLevel~xresidAD, data=traits.tree, lambda = "ML")
summary(combined.prim3)

combined.prim4=pgls(xTrophicLevel~xresidVA, data=traits.tree, lambda = "ML")
summary(combined.prim4)

##make plots, Use 1 column
acuity=read.csv("/Users/JosephOuta/Desktop/acuity.csv", header=T, row.names = 1)
acuity=as.matrix(acuity)[,1]
acuity

acuityphy1 = contMap(pruned.tree, acuity, outline=FALSE)
plot.acuity=plot(acuityphy1, outline = FALSE)

##make plots, Use 1 column
size=read.csv("/Users/JosephOuta/Desktop/eyesize.csv", header=T, row.names = 1)
size=as.matrix(size)[,1]
size

sizephy1 = contMap(pruned.tree, size, outline=FALSE)
plot.size=plot(sizephy1, outline = FALSE)

#second plot
plumageM=read.csv("/Users/michaelsheehan/Downloads/plumageM.csv", header=T, row.names = 1)
plumageM=as.matrix(plumageM)[,1]
plumageM

obj2=contMap(bird1, plumageM, outline=FALSE)
plot.plumageM=plot(obj2, direction="leftwards", outline = FALSE)