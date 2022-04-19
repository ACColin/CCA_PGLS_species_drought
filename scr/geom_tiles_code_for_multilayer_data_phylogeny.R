library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)
library(tidyverse)

# load data from TDbook, including tree_hmptree, 
# df_tippoint (the abundance and types of microbes),
# df_ring_heatmap (the abundance of microbes at different body sites),
# and df_barplot_attr (the abundance of microbes of greatest prevalence)
tree <- tree_hmptree
dat1 <- df_tippoint
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# adjust the order
dat2$Sites <- factor(dat2$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)","Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
dat3$Sites <- factor(dat3$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)", "Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
# extract the clade label information. Because some nodes of tree are annotated to genera,
# which can be displayed with high light using ggtree.
nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) + # plots the tree
  geom_hilight(data=nodedf, mapping=aes(node=node), # highlights clades, here in grey
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf, # prints the name of the clades from the dataframe that are highlighted
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=1.4,
                fontface="italic"
  )

p <- p %<+% dat1 + geom_star( # ggstar:geom_star gives a geometric shape to the tips like ggplot2:geom_point
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p <- p + new_scale_fill() + # ggnewscale:new_scale
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
  ) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )





  ######      ###          #####
 ##    Trying to make my own ##
######      ###          #####
tree <- read.tree("plots/final_tree_uniqueID.tre")
info <- read.csv("data/DroughtImpact_SummaryBy_Completetraits.csv")

dat <- info %>% 
  select(UniqueID, height_max_m) %>% 
  as.data.frame() %>% 
  glimpse()
dat2 <- info %>% 
  select(UniqueID, AICorrected) %>% 
  as.data.frame() %>% 
  glimpse()
dat3 <- info %>% 
  select(UniqueID, Drought2020Scale.mean) %>% 
  as.data.frame() %>% 
  glimpse()

## layer 1
p <- ggtree(tree, layout="fan", size=0.5, open.angle=10) +
  geom_fruit(data=dat, geom=geom_tile, mapping=(aes(y=UniqueID, fill=height_max_m)),color = "white", offset = 0.04, size = 0.2) #+
  # Color scale for maximum height:
  scale_color_viridis_c(option = "D")
print(p)
ggsave("plots/test_geomfruit_param.pdf", width = 50, height = 120, limitsize = FALSE)

p <- ggtree(tree, layout="fan", open.angle = 15, size = 0.2) +
  geom_tippoint(dat=dat, aes(colour=height_max_m)) +
  geom_tiplab(aes(colour=height_max_m),
              align = T, linetype = 3, size = 1, linesize = 0.2, show.legend = F
              ) +
  scale_color_manual(guide=guide_legend(keywidth=0.5, keyheight=0.5, order=2, override.aes=list(size=2,alpha=1)))

## layer 2
p1 <- p +
  new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile, mapping=aes(y=UniqueID, fill=AICorrected)) +
  # Color scale applied to geoms added after new_scale_fill()
  scale_color_viridis_c(option = "D")
print(p1)

# layer 3
p2 <- p1 +
  new_scale_fill() +
  geom_fruit(data=dat3, geom=geom_tile, mapping=aes(y=UniqueID, fill=Drought2020Scale.mean)) +
  scale_color_gradientn(colours=c("blue","red"))
print(p2)
ggsave("plots/test_geomfruit_param_LAYERS.pdf", width = 50, height = 120, limitsize = FALSE)
