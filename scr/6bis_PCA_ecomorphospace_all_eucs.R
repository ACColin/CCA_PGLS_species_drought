

################# Different PCA methods #################

########### PCA with cowplot
library(tidyverse)
library(cowplot)

t=read.csv("outputs/all_eucs_PCA_dataframe.csv", header = T, check.names=F) %>% 
  dplyr::select(Drought2020Scale.mean,AP,AICorrected,PDryQ,
                PrecDeficit,height_max_m,Habit_num,
                RegenStrategysimplified_num,LeafAreaEst_cm2,Aridclass_num) %>% 
  na.omit()
p=prcomp(t)
pct=paste0(colnames(p$x)," (",sprintf("%.1f",p$sdev/sum(p$sdev)*100),"%)")
p2=as.data.frame(p$x)
p2$k=factor(cutree(hclust(dist(t)),k=12))
load=p$rotation

plots=lapply(seq(1,7,2),function(i){
  x=sym(paste0("PC",i))
  y=sym(paste0("PC",i+1))
  
  mult=min(max(p2[,i])/max(load[,i]),max(p2[,i+1])/max(load[,i+1]))
  colors=hcl(head(seq(15,375,length=length(unique(p2$k))+1),-1),120,50)
  
  ggplot(p2,aes(!!x,!!y))+
    geom_segment(data=load,aes(x=0,y=0,xend=mult*!!x,yend=mult*!!y),arrow=arrow(length=unit(.3,"lines")),color="gray60",size=.4)+
    annotate("text",x=(mult*load[,i]),y=(mult*load[,i+1]),label=rownames(load),size=2.5,vjust=ifelse(load[,i+1]>0,-.5,1.4))+
    geom_polygon(data=p2%>%group_by(k)%>%slice(chull(!!x,!!y)),aes(color=k,fill=k),size=.3,alpha=.2)+
    geom_point(aes(color=k),size=.6)+
    #geom_text(aes(label=rownames(t),color=k),size=2.5,vjust=-.6)+
    # ggrepel::geom_text_repel(aes(label=rownames(t),color=k),max.overlaps=Inf,force=5,size=2.2,min.segment.length=.1,segment.size=.2)+
    labs(x=pct[i],y=pct[i+1])+
    scale_x_continuous(breaks=seq(-100,100,20),expand=expansion(mult=.06))+
    scale_y_continuous(breaks=seq(-100,100,20),expand=expansion(mult=.06))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    theme(aspect.ratio=1,
          axis.text=element_text(color="black",size=6),
          axis.text.x=element_text(margin=margin(.2,0,0,0,"cm")),
          axis.text.y=element_text(angle=90,vjust=1,hjust=.5,margin=margin(0,.2,0,0,"cm")),
          axis.ticks=element_line(size=.3,color="gray60"),
          axis.ticks.length=unit(-.13,"cm"),
          axis.title=element_text(color="black",size=8),
          legend.position="none",
          panel.background=element_rect(fill="white"),
          panel.border=element_rect(color="gray60",fill=NA,size=.4),
          panel.grid=element_blank())
})

plot_grid(plotlist=plots)
ggsave("a.png",height=12,width=12)



########### Script from Victoria
pcaDF = read.csv("outputs/all_eucs_PCA_dataframe.csv", header = T)
PCADF3 <-  pcaDF %>% dplyr::select(AICorrected, UniqueID, Drought2020Scale.mean, Taxonsimplified_DN,
                                   Aridclass_num, RegenStrategysimplified_num, Habit_num, height_max_m,
                                   Wooddensity.g.cm3..mean,
                                   d13C.mean, totalN.mean)

PCADF3 <- PCADF3 %>% rename(Height=height_max_m, Regen=RegenStrategysimplified_num,
                            WD=Wooddensity.g.cm3..mean,
                            dC13=d13C.mean, Nleaf=totalN.mean)
which(colnames(PCADF3)=="Height")

#only take continuous var
pcaCCA2<- PCA(PCADF3[-c(2:11)], scale.unit=TRUE,  graph=FALSE)
And here different plots, the last bit is to add the ellipses
plot

#positively correlated variables are grouped together.

#Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).

#The distance between variables and the origin measures the quality of the variables on the factor map. Variables that are away from the origin are well represented on the factor map.

fviz_pca_var(pcaCCA2, col.var = "black")

#plot corr with dimensions

corrplot(var2$cos2, is.corr=FALSE)

#plot with individuals (uniqueid) plot by average of drought impact but not variables included

fviz_pca_ind(pcaCCA2.1, col.ind =PCADF_notreg $Drought2020Scale.mean, 
             gradient.cols = c("blue", "red"),
             legend.title = "Drought Impact",  geom="point", label ="var")

fviz_pca_ind(pcaCCA2.1,  geom.ind = "point", col.ind=PCADF_notreg $Drought2020Scale.mean,
             gradient.cols = c("blue", "red"),
             legend.title = "Drought Impact", label ="var")

##this includes variables and individuals-uniqueid coloured by drought impact 2020

fviz_pca_biplot(pcaCCA2.1,  label= "var",  col.var="black", labelsize = 5,
                col.ind =PCADF_notreg $Drought2020Scale.mean,
                gradient.cols = c("blue", "red"), pointsize = 2.5, legend.title = "Drought Impact", title = NULL) + theme(text = element_text(size = 15),
 axis.title = element_text(size = 15),axis.text = element_text(size = 15))

ggsave("PCA_aitraitsnotreg_22-8-22.jpg", width = 40, height = 20, units=c("cm"), dpi=400)



#to add ellipse but only possible with categorical

fviz_pca_biplot(pcaCCA2.1, label="var",  col.var="black", labelsize = 5, habillage=PCADF3$DroughtClass,
                addEllipses=TRUE, ellipse.level=0.95, pointsize = 2.5, legend.title = "Drought Impact", title = NULL) + theme(text = element_text(size = 15),
     axis.title = element_text(size = 15),
     axis.text = element_text(size = 15))

ggsave("PCA_aitraitsEllipsesDrought_4-8-22.jpg", width = 40, height = 20, units=c("cm"), dpi=400)