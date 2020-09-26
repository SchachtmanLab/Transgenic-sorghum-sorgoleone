BiocManager::install("vegan")
BiocManager::install("ggplot2")
BiocManager::install("ape")
library(vegan)
library(ggplot2)
library(ape)
setwd("~/Desktop/transgenic_sorghum")

#CAP analysis
sorg.map <- read.delim("map_2016_2015_2014.txt") #change this
row.names(sorg.map)<-sorg.map$SampleID
sorg.map$Type <- factor(sorg.map$Type)
sorg.map$Year <- factor(sorg.map$Year)
sorg.map$Plantstage <- factor(sorg.map$Plantstage)
sorg.map$Genotype <- factor(sorg.map$Genotype)
gh.bray <- read.table("bray_curtis_dm.txt", header = T, row.names = 1)
gh.bray <- gh.bray[match(row.names(sorg.map), row.names(gh.bray)), match(row.names(sorg.map), colnames(gh.bray))]
bray.cap.whole <- capscale(as.dist(gh.bray) ~   Genotype + Condition(Plantstage),  data = sorg.map, add = T) #change

#PERMONOVA analysis to test whether sorgoleone can significantly affect the microbial communities structure
anova(bray.cap.whole, by = "terms")
summary(bray.cap.whole)
percent_explained <- bray.cap.whole$CCA$eig / sum(bray.cap.whole$CCA$eig) * 100

#use adnois() to obtain the explained variance by genotype
adonis(gh.bray ~Genotype, data = sorg.map, strata=paste(sorg.map$Plantstage), add=T)


#

setwd("~/Desktop/transgenic_sorghum")
setwd("~/Desktop/peng_wang_06052018/github")
AllPCA <-read.csv("test_mu_plots_2016.csv",row.names=1,check.names=FALSE) 

AllPCA$Plantstage  <- factor(AllPCA$Plantstage, levels = c("7 leaf stage","Boot","Grain fill"))
AllPCA$Type  <- factor(AllPCA$Type, levels = c("Endosphere","Rhizosphere","Soil"))

color_flag <- c("#2c7fb8","#fdbb84") ##fee8c8
p <- ggplot(AllPCA,aes(CAP1,MDS1,color=Genotype, shape=Genotype, stroke = 1))
p + geom_point(data=AllPCA,aes(CAP1,MDS1,color=Genotype,fill=Genotype, stroke = 1),size=5)+
  scale_colour_manual(values=c("black","black","white","white","white","white","white","white","grey","blue","white","white","white","white","grey","white","white"))+ #DIMBOA high
  scale_shape_manual(values=c(22,21,17, 19,15, 15, 15, 17,16,15, 1,17,1, 15, 2,16,17,18,3,4,0,1,5,6,15, 16, 17, 18,20))+
  scale_fill_manual(values=color_flag) +
  facet_grid(Plantstage  ~ Type, scales = "free") + #####try again make all axis free
  
  theme_bw() +
  #theme(panel.spacing = unit(0.00, "lines")) +
  #scale_linetype_manual(values = c("SBR" = "dotted", "SSR"="dashed"), name="Segment")+
  theme(strip.text.x = element_text(size = 15, colour = "black", angle = 0)) +
  theme(strip.text.y = element_text(size = 16, colour = "black", angle = 270)) +
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_rect(fill="grey", colour="white",size=0.5)) +
  #theme(axis.text.x = element_text(angle = 180, hjust = 1, size=14,color="black",face="bold")) +
  theme(axis.text.x = element_text(size=15,color="black",face="bold")) +
  theme(axis.text.y = element_text(size=15,color="black",face="bold")) +
  theme(axis.title.x = element_text(face = "bold",size = rel(1.4))) + 
  #theme(axis.title.x = element_text(size = rel(2), angle = 90)) + 
  theme(axis.title.y = element_text(face = "bold",size = rel(1.4))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.text = element_text(face = "bold",size = 15)) +
  theme(legend.title  = element_text(face = "bold",size = 15)) +
  theme(legend.position = "bottom") +
  guides(shape=guide_legend(override.aes = list(size=3, stroke=1))) +
  geom_vline(xintercept = 0, alpha = 0.7, linetype="dotted", 
             color = "black", size=0.5) +
  geom_hline(yintercept = 0, alpha = 0.7, linetype="dotted", 
             color = "black", size=0.5) +
  theme(panel.spacing = unit(0.6, "lines")) +
  #stat_ellipse(type = "euclid")+
  labs(x = "CAP1", y="MDS1") +
  #labs(fill="Genotype")   +
  theme(aspect.ratio=0.95) 

