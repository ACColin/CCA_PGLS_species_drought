## Quick mapping of drought incidence on original Thornhill for AES2021

# Keeping the original Thornhill topology (and taxa name) and mapping the mean of drought incidence at the tip.


DI <- read.csv("../data/CCAraw_Thornhill_Traits.csv", header = T) # getting the drought data in
DI %>% 
  group_by(UniqueID) %>% 
  summarise(average = mean(Drought2020Scale))

Thornhill.taxa <- read.csv("../data/drought.incidence.csv", header = T) # name to match with thornhill

all <- left_join(Thornhill.taxa, DI, by = "UniqueID") %>% # matching name and drought incidence data based on unique ID
  glimpse()
view(all)

info <- all %>% # averaging values per species to add data on the tree
  group_by(thornhill_bionomial_tip_match) %>% 
  summarise(drought.incidence = mean(Drought2020Scale))

#write.csv(info,"../outputs/info.csv")
# getting the tip labels from thornhill phylogeny (to further get the data at the tips)
thornhill_labels <- read.csv("../data/tree_labels.csv")
tree <- read.tree("../data/Eucalypts_ML2_dated_r8s.tre")
info <- read.csv("../outputs/info.csv")
# plotting the data
tree$tip.label<-thornhill_labels$name4[match(tree$tip.label, thornhill_labels$original)] # changing tip label to match victoria species name

ggtree(tree)

dotTree(tree,x,length=10,ftype="i", type = 'fan')

p <- ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=drought.incidence)) +
  scale_color_gradientn(colours=c("blue","red")) +
  geom_tiplab(aes(color=drought.incidence), hjust = -.1)
print(p)
ggsave("../plots/drought.incidence.mapped.colours.pdf", width = 50, height = 120, limitsize = FALSE)

p1 <- ggtree(tree, aes(color=drought.incidence), layout = 'circular', ladderize = F, continuous = T, size = 2) + scale_color_gradientn(colours=c("red", "orange", "green","cyan","blue")) +
  geom_tiplab(aes(color=drought.incidence), hjust = -.1) +
  xlim(0, 1,2) +
  theme(legend.position = x(.05, .85))