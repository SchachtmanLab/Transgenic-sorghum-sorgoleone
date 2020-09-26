setwd("~/Desktop/Transgenic_sorghum/Network/family_rel/AUG_RNAi_core_soil")
setwd("~/Desktop/Transgenic_sorghum/Network/family_rel/RNAi_all")
setwd("~/Desktop/Transgenic_sorghum/Network/family_rel/pvals_RNAi_RHZ")
setwd("~/Desktop/Transgenic_sorghum/Network/family_rel")
setwd("~/Desktop/Transgenic_sorghum/Network/family_near_soil_Jul_RNAi")
setwd("~/Desktop/Transgenic_sorghum/Network/family_near_RHZ_Jul_WT")
setwd("~/Desktop/Transgenic_sorghum/Network/family_near_RHZ_Jul_RNAi")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/core_DIMBOA")
setwd("~/Desktop/Kaiju_tax/Family_PHYLUM/statistical analysis/only_diff_for_network")

##move the column sequence in the data frame
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}


library(reshape2)
library(dplyr)
d <- read.csv("spear_WT_core_Aug.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_family_RNAi_all.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_family_WT_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_family_RNAi_JUL_soil.out.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("otu_wt_aug_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_AUG_RNAi_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_AUG_RNAi_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_JUL_WT_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_JUL_RNAi_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_SEP_WT_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_SEP_RNAi_soil.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_JUL_WT_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_JUL_RNAi_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/AUG_RHZ_OTU/WT/WT_new")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/core_DIMBOA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/core_DIMBOA_GABA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/pvals_Core_DIMBOA_GABA_GO_DEGs")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock_2nd/analysis/cor_Core_DIMBOA_absolu")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock_2nd/analysis/cor_Core_GABA_absolu")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/DIMBOA-network")

#2020MAY#make average for different dataframe (1) need to make the index for each df.
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/GABA-network/Netowrk_Havelock_GABA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/LOW_GABA_for_DIMBOA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/High_GABA_for_DIMBOA")

d10 <- read.csv("core_table_80_M_GABA_10_ana.csv", header=TRUE,row.names=1, sep=",")
d37 <- read.csv("core_table_80_M_GABA_37_GABA_ana.csv", header=TRUE,row.names=1, sep=",")
d109 <- read.csv("core_table_80_M_GABA_109_GABA_ana.csv", header=TRUE,row.names=1, sep=",")

#for low GABA low DIMBOA
d6 <- read.csv("P006_low_gaba_low_dimboa.csv", header=TRUE,row.names=1, sep=",")
d84 <- read.csv("P084_low_gaba_low_dimboa.csv", header=TRUE,row.names=1, sep=",")

###DIMBOA_FUNGI_BACTERIA

d <- read.csv("cor_bac_fungi_DIMBOA_family_ana.csv", header=TRUE,row.names=1, sep=",")


##make average for different data frames
library(data.table)
head(d10[,c(1:14)],10)
mean_3s=rbindlist(list(d10,d37,d109))[,lapply(.SD,mean), list(OUT_ID)] ## two different column names:Lat Lon
write.csv(mean_3s, file = "cor_OTU_mean_10_37_109_M_GABA.csv")
dt = rbindlist(list(d10,d37,d109))
dt[,print(.SD), list(OUT_ID)]  #Lon, Lat


mean_2s=rbindlist(list(d6,d84))[,lapply(.SD,mean), list(OUT_ID)] ## two different column names:Lat Lon
write.csv(mean_2s, file = "cor_OTU_mean_6_84_LGABA_LDIMBOA.csv")
dt = rbindlist(list(d10,d37,d109))
dt[,print(.SD), list(OUT_ID)]  #Lon, Lat

###############################



d <- read.csv("cor_OTU_AUG_WT_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_AUG_RNAi_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_OTU_SEP_RNAi_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_core_DIMBOA.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_core_DIMBOA_GABA.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")


d <- read.csv("cor_core_DIMBOA_GABA.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_Core_GABA_DIMBOA_DEGs_from_GO.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_Core_DIMBOA_absolu.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_Core_GABA_abosu.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")

#############################
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/pvals_Medium_DIMBOA_Root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_Medium_DIMBOA_Root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_Low_DIMBOA_Root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_High_GABA_ROOT")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_Medium_GABA_Root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_low_GABA_root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Hi_DIMBOA_rhz")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Medium_DIMBOA_rhz")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Low_DIMBOA_rhz")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Low_GABA_rhz")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Hi_GABA_rhz")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/Medium_GABA_rhz")
library(reshape2)
#2018Havelock reanalysis
d <- read.delim("core_Medium_GABA.txt")
#d <- read.csv("High_DIMBOA_Root.csv", header=TRUE,row.names=1, sep=",")
#p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")
p <- read.delim("pvals.two_sided.txt")
x=melt(p)
write.csv(x, file="pvals.two_sided_1.csv")

x=melt(d)
write.csv(x, file="High_DIMBOA_Root_1.csv")

s <- read.csv("High_DIMBOA_Root_1.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided_1.csv", header=TRUE,row.names=1, sep=",")

#sp1 <- cbind(s,p, by=c("variable"))
sp1 <- cbind(s,p, by=c("OTU.ID"))
colnames(sp1)=c("tax1","tax2","r","id1","id2","p","by")
sp1 <- sp1 %>%  
  select("tax1","tax2","r","p")
write.csv(sp1, file="spear_cor_Core_High_DIMBOA_Root_2.csv")
#############################
#############################
#no p value
d <- read.csv("cor_OTU_mean_10_37_109_M_GABA.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("cor_core_table_80_L_GABA_B73_ana.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("cor_H_DIMBOA_all.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("PB73_low_gaba_high_dimboa.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_P018.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("cor_P100.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("cor_bac_fungi_DIMBOA_family_ana.csv", header=TRUE,row.names=1, sep=",")



###
###




x=melt(p)
write.csv(x, file="pvals.two_sided_1.csv")

x=melt(d)
write.csv(x, file="High_DIMBOA_Root_1.csv")

write.csv(x, file="spear_cor_Core_GABA_mean_10_37_109.csv")
write.csv(x, file="spear_cor_Core_GABA_Low_B73.csv")
write.csv(x, file="spear_cor_Core_LGABA_dimboa_H_B73.csv")
write.csv(x, file="spear_cor_Core_HGABA_Hdimboa_100.csv")

write.csv(x, file="spear_cor_Core_DIMBOA_fungi_bac.csv")
#myvars <- names(mydata) %in% c("v1", "v2", "v3")
#newdata <- mydata[!myvars]
#########################
#########################As example
##OTU-WT-AUG, DIMBOA, DIMBOA_GABA_DEGs
s <- read.csv("spear_cor_Core_GABA_absolu.csv", header=TRUE,row.names=1, sep=",")
nrow(s)
s1<-s %>% slice(1:1000000)
s2<-s %>% slice(1000001:2000000)
s3<-s %>% slice(2000001:3000000)
s4<-s %>% slice(3000001:4000000)
s5<-s %>% slice(4000001:4355569)

s6<-s %>% slice(5000001:6000000)
s7<-s %>% slice(6000001:7000000)
s8<-s %>% slice(7000001:8000000)
s9<-s %>% slice(8000001:9000000)
s10<-s %>% slice(9000001:10000000)
s11<-s %>% slice(10000001:11000000)
s12<-s %>% slice(11000001:11519236)


write.csv(s1, file="spear_cor_Core_GABA_absolu_1.csv")                                            
write.csv(s2, file="spear_cor_Core_GABA_absolu_2.csv")                                            
write.csv(s3, file="spear_cor_Core_GABA_absolu_3.csv")  
write.csv(s4, file="spear_cor_Core_GABA_absolu_4.csv")  
write.csv(s5, file="spear_cor_Core_GABA_absolu_5.csv")


write.csv(s6, file="spear_cor_core_DIMBOA_GABA_DEG_6.csv")                                            
write.csv(s7, file="spear_cor_core_DIMBOA_GABA_DEG_7.csv")                                            
write.csv(s8, file="spear_cor_core_DIMBOA_GABA_DEG_8.csv")  
write.csv(s9, file="spear_cor_core_DIMBOA_GABA_DEG_9.csv")  
write.csv(s10, file="spear_cor_core_DIMBOA_GABA_DEG_10.csv")  
write.csv(s11, file="spear_cor_core_DIMBOA_GABA_DEG_11.csv")                                            
write.csv(s12, file="spear_cor_core_DIMBOA_GABA_DEG_12.csv")                                            

##OTU-RNAI-AUG
s <- read.csv("spear_RNAI-soil-aug-otu.csv", header=TRUE,row.names=1, sep=",")
nrow(s)
s1<-s %>% slice(1:1000000)
s2<-s %>% slice(1000001:2000000)
s3<-s %>% slice(2000001:2259009)

write.csv(s1, file="spear_RNAI-soil-aug-otu_1.csv")                                            
write.csv(s2, file="spear_RNAI-soil-aug-otu_2.csv")                                            
write.csv(s3, file="spear_RNAI-soil-aug-otu_3.csv")  
#########################
#########################As example
##OTU-soil-WT-JUL, RNAi, and soil-WT-SEP, RNAi
s <- read.csv("spear_WT-RHZ-AUG-otu.csv", header=TRUE,row.names=1, sep=",")
nrow(s)
s1<-s %>% slice(1:1000000)
s2<-s %>% slice(1000001:2000000)
s3<-s %>% slice(2000001:2307361)
#s3<-s %>% slice(2000001:2105401) #jul-soil
s3<-s %>% slice(2000001:3000000)
s4<-s %>% slice(3000001:3279721)

write.csv(s1, file="spear_RNAi-RHZ-JUL-otu_1.csv")                                            
write.csv(s2, file="spear_RNAi-RHZ-JUL-otu_2.csv")                                            
write.csv(s3, file="spear_RNAi-RHZ-JUL-otu_3.csv") 
write.csv(s4, file="spear_RNAi-RHZ-JUL-otu_4.csv")  
##
##OTU-RHZ-WT-AUG, RNAi
s <- read.csv("spear_WT-RHZ-AUG-otu.csv", header=TRUE,row.names=1, sep=",")
nrow(s)
s1<-s %>% slice(1:1000000)
s2<-s %>% slice(1000001:1830609)

s2<-s %>% slice(1000001:2000000)
#s3<-s %>% slice(2000001:2105401) #jul-soil
s3<-s %>% slice(2000001:2307361)

write.csv(s1, file="spear_WT-RHZ-AUG-otu_1.csv")                                            
write.csv(s2, file="spear_WT-RHZ-AUG-otu_2.csv")                                            
write.csv(s3, file="spear_WT-RHZ-AUG-otu_3.csv") 


##OTU-RNAI-JUL
s <- read.csv("High_DIMBOA_Root.csv", header=TRUE,row.names=1, sep=",")
nrow(s)
s1<-s %>% slice(1:1000000)
s2<-s %>% slice(1000001:2000000)
s3<-s %>% slice(2000001:2259009)

write.csv(s1, file="spear_RNAI-soil-aug-otu_1.csv")                                            
write.csv(s2, file="spear_RNAI-soil-aug-otu_2.csv")                                            
write.csv(s3, file="spear_RNAI-soil-aug-otu_3.csv")  

#names
#n <- 3
#s4 <- data.frame(s3,i=rep(1:n,ea=NROW(s3)))
#write.csv(s4, file="spear_wt-soil-aug-otu_4.csv")  

name=d[,1]
#name <- select(d,starts_with("OTU_id"))
#name$ID <- seq.int(nrow(name))
#moveme(names(name), "ID first")
#name <- tibble::rowid_to_column(name, "ID")

write.csv(name, file="names.csv") 
name <- read.csv("names.csv", header=T,row.names=1, sep=",")
#name <- read.csv("names.csv", header=TRUE,row.names=1, sep=",")
n <- 200
names1 <- data.frame(name,i=rep(1:n,ea=NROW(name)))
write.csv(names1, file="names_n.csv")  

#JUL_SOIL_WT and RNAi and DIMBOA
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")
nrow(p)
p1<-p %>% slice(1:1000000)
p2<-p %>% slice(1000001:2000000)
p3<-p %>% slice(2000001:3000000)
p4<-p %>% slice(3000001:4000000)
p5<-p %>% slice(4000001:4355569)


p6<-p %>% slice(5000001:6000000)
p7<-p %>% slice(6000001:7000000)
p8<-p %>% slice(7000001:8000000)
p9<-p %>% slice(8000001:9000000)
p10<-p %>% slice(9000001:10000000)
p11<-p %>% slice(10000001:11000000)
p12<-p %>% slice(11000001:11519236)




write.csv(p1, file="pvals.two_sided_1.csv")                                            
write.csv(p2, file="pvals.two_sided_2.csv")  
write.csv(p3, file="pvals.two_sided_3.csv")    
write.csv(p4, file="pvals.two_sided_4.csv")   
write.csv(p5, file="pvals.two_sided_5.csv")  
write.csv(p6, file="pvals.two_sided_6.csv")                                            
write.csv(p7, file="pvals.two_sided_7.csv")  
write.csv(p8, file="pvals.two_sided_8.csv")    
write.csv(p9, file="pvals.two_sided_9.csv")   
write.csv(p10, file="pvals.two_sided_10.csv")  
write.csv(p11, file="pvals.two_sided_11.csv")                                            
write.csv(p12, file="pvals.two_sided_12.csv")  



head(s1, n=3)
tail(s3, n=10)
tail(s, n=10)
head(p1, n=3)
tail(s5, n=10)

##merge-1 #for july soil WT and RNAi. Merge the p value and R value
####this is for the Havelock 2018 data analysis
s <- read.csv("High_DIMBOA_Root_1.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided_1.csv", header=TRUE,row.names=1, sep=",")

#sp1 <- cbind(s,p, by=c("variable"))
sp1 <- cbind(s,p, by=c("OTU.ID"))
colnames(sp1)=c("tax1","tax2","r","id1","id2","p","by")
sp1 <- sp1 %>%  
  select("tax1","tax2","r","p")
write.csv(sp1, file="spear_cor_Core_High_DIMBOA_Root_2.csv")

sp1 <- cbind(s1,p1, by=c("variable"))
sp2 <- cbind(s2,p2, by=c("variable"))
sp3 <- cbind(s3,p3, by=c("variable"))
sp4 <- cbind(s4,p4, by=c("variable"))
sp5 <- cbind(s5,p5, by=c("variable"))
sp6 <- cbind(s6,p6, by=c("variable"))
sp7 <- cbind(s7,p7, by=c("variable"))
sp8 <- cbind(s8,p8, by=c("variable"))
sp9 <- cbind(s9,p9, by=c("variable"))
sp10 <- cbind(s10,p10, by=c("variable"))
sp11 <- cbind(s11,p11, by=c("variable"))
sp12 <- cbind(s12,p12, by=c("variable"))


write.csv(sp1, file="spear_cor_Core_High_DIMBOA_Root_1.csv")
write.csv(sp2, file="spear_cor_Core_GABA_absolu_sp2.csv")
write.csv(sp3, file="spear_cor_Core_GABA_absolu_sp3.csv")
write.csv(sp4, file="spear_cor_Core_GABA_absolu_sp4.csv")
write.csv(sp5, file="spear_cor_Core_GABA_absolu_sp5.csv")

write.csv(sp6, file="spear_cor_core_DIMBOA_GABA_DEG_sp6.csv")
write.csv(sp7, file="spear_cor_core_DIMBOA_GABA_DEG_sp7.csv")
write.csv(sp8, file="spear_cor_core_DIMBOA_GABA_DEG_sp8.csv")
write.csv(sp9, file="spear_cor_core_DIMBOA_GABA_DEG_sp9.csv")
write.csv(sp10, file="spear_cor_core_DIMBOA_GABA_DEG_sp10.csv")
write.csv(sp11, file="spear_cor_core_DIMBOA_GABA_DEG_sp11.csv")
write.csv(sp12, file="spear_cor_core_DIMBOA_GABA_DEG_sp12.csv")

##merge-2
#As example
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network/OTU_soil_WT_aug/merge")
m1 <- read.csv("spear_cor_Core_GABA_absolu_sp1.csv", header=TRUE,row.names=1, sep=",")
m2 <- read.csv("spear_cor_Core_GABA_absolu_sp2.csv", header=TRUE,row.names=1, sep=",")
m3 <- read.csv("spear_cor_Core_GABA_absolu_sp3.csv", header=TRUE,row.names=1, sep=",")
m4 <- read.csv("spear_cor_Core_GABA_absolu_sp4.csv", header=TRUE,row.names=1, sep=",")
m5 <- read.csv("spear_cor_Core_GABA_absolu_sp5.csv", header=TRUE,row.names=1, sep=",")


m3 <- read.csv("spear_WT-RHZ-AUG-otu_sp3.csv", header=TRUE,row.names=1, sep=",")
#m123 <- merge(m1,m2, m3, by=c("ID")) # merge by column

#JUL_SOIL_WT and RNAi, SEP_SOIL_WT, RNAi and DIMBOA for input in igraph
m1 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp1.csv", header=TRUE,row.names=1, sep=",")
m2 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp2.csv", header=TRUE,row.names=1, sep=",")
m3 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp3.csv", header=TRUE,row.names=1, sep=",")
m4 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp4.csv", header=TRUE,row.names=1, sep=",")
m5 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp5.csv", header=TRUE,row.names=1, sep=",")
m6 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp6.csv", header=TRUE,row.names=1, sep=",")
m7 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp7.csv", header=TRUE,row.names=1, sep=",")
m8 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp8.csv", header=TRUE,row.names=1, sep=",")
m9 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp9.csv", header=TRUE,row.names=1, sep=",")
m10 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp10.csv", header=TRUE,row.names=1, sep=",")
m11 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp11.csv", header=TRUE,row.names=1, sep=",")
m12 <- read.csv("spear_cor_core_DIMBOA_GABA_DEG_sp12.csv", header=TRUE,row.names=1, sep=",")

#merge by rows
m_all <- rbind(m1, m2, m3, m4, m5) 
nrow(m1)
nrow(m2)
nrow(m3)
nrow(m4)
nrow(m5)
nrow(m_all)
write.csv(m_all, file="spear_cor_Core_GABA_absolu_sp12345.csv")

m12345 <- read.csv("spear_cor_Core_GABA_absolu_sp12345.csv", header=TRUE,row.names=1, sep=",")
tail(m12345, n=5)



################

#use BiocManager::install() to install the following R packages
library(ggplot2)
library(corrr)
library(ggraph)
library(tidyr)
library(tidygraph)
library(igraph)
library(reshape)
###OTU-JUL-SOIL-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_JUL_WT")
d <- read.csv("spear_WT-soil-JUL-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_WT-soil-JUL-otu_sp123_degree_posi_0.6.csv", header=TRUE,row.names=1, sep=",")
#stacknodes <- read.csv("spear_WT-soil-JUL-otu_sp123_degree.csv", header=TRUE,row.names=1, sep=",")
###OTU-JUL-SOIL-RNAi
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_JUL_RNAi")
d <- read.csv("spear_RNAi-soil-JUL-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAi-soil-JUL-otu_sp123_degree.csv", header=TRUE,row.names=1, sep=",")
###OTU-AUG-SOIL-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_WT_aug/merge")
d <- read.csv("spear_wt-soil-aug-otu_123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_wt-soil-aug-otu_123_degree_posi_0.6.csv", header=TRUE,row.names=1, sep=",")
graph_cors <- d %>%
  filter(r > .60) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(abs(p) <= .01) %>%
  graph_from_data_frame(vertices = stacknodes, directed = FALSE)
###OTU-AUG-SOIL-RNAI
setwd("/Users/peng/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_RNAi_aug/merge")
d <- read.csv("spear_RNAI-soil-aug-otu_123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAI-soil-aug-otu_123_degree.csv", header=TRUE,row.names=1, sep=",")
graph_cors <- d %>%
  filter(abs(r) > .80) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(abs(p) <= .01) %>%
  graph_from_data_frame(vertices = stacknodes, directed = FALSE)
### OTU-SEP-SOIL-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_WT_sep")
d <- read.csv("spear_WT-soil-SEP-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_WT-soil-SEP-otu_sp123_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU-SEP-SOIL-RNAi
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_soil/OTU_soil_RNAi_sep")
d <- read.csv("spear_RNAi-soil-SEP-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAi-soil-SEP-otu_sp123_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU-JUL-RHZ-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/JUL_RHZ_OTU/WT")
d <- read.csv("spear_WT-RHZ-JUL-otu_sp1234.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_WT-RHZ-JUL-otu_sp1234_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU-JUL-RHZ-RNAi
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/JUL_RHZ_OTU/RNAi")
d <- read.csv("spear_RNAi-RHZ-JUL-otu_sp1234.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAi-RHZ-JUL-otu_sp1234_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU-AUG-RHZ-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/AUG_RHZ_OTU/WT/WT_new")
d <- read.csv("spear_WT-RHZ-AUG-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_WT-RHZ-AUG-otu_sp123_0.8_degree.csv", header=TRUE,row.names=1, sep=",") #spear_WT-RHZ-AUG-otu_sp123_degree.csv
#OTU-AUG-RHZ-RNAi
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/AUG_RHZ_OTU/RNAi")
d <- read.csv("spear_RNAi-RHZ-AUG-otu_sp123.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAi-RHZ-AUG-otu_sp123_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU-SEP-RHZ-WT
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/SEP_RHZ_OTU/WT")
d <- read.csv("spear_WT-RHZ-SEP-otu_sp12.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_WT-RHZ-SEP-otu_sp12_degree.csv", header=TRUE,row.names=1, sep=",")
#OTU_SEP-RHZ-RNAi
setwd("~/Desktop/Transgenic_sorghum/Network/OTU_network_RHZ/SEP_RHZ_OTU/RNAi")
d <- read.csv("spear_RNAi-RHZ-SEP-otu_sp12.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_RNAi-RHZ-SEP-otu_sp12_degree.csv", header=TRUE,row.names=1, sep=",")

###OTU-DIMBOA_OTU
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/core_DIMBOA")
d <- read.csv("spear_cor_core_DIMBOA_sp12345.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_DIMBOA_sp12345_degree_0.7.csv", header=TRUE,row.names=1, sep=",")

setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/core_DIMBOA_GABA")
d <- read.csv("spear_cor_core_DIMBOA_GABA_sp12345.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_DIMBOA_GABA_sp12345_degree_0.7.csv", header=TRUE,row.names=1, sep=",")

##for the all DIMBAO and GABA and DEGs
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock/pvals_Core_DIMBOA_GABA_GO_DEGs")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/Network/Network_2018_Havelock_2nd/analysis/cor_Core_GABA_absolu")
d <- read.csv("spear_cor_Core_GABA_absolu_sp12345.csv", header=TRUE,row.names=1, sep=",")

##2020 redo
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/GABA-network/Netowrk_Havelock_GABA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/DIMBOA-network")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/LOW_GABA_for_DIMBOA")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_ANALYSIS_09302019/ALL_ROOT_RHZ/DESEQ2/Network/High_GABA_for_DIMBOA")
#fungi and bacteria
setwd("~/Desktop/Kaiju_tax/Family_PHYLUM/statistical analysis/only_diff_for_network")


### 2020 August redo redo
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_High_DIMBOA_Root")
d <- read.csv("spear_cor_Core_High_DIMBOA_Root_2.csv", header=TRUE,row.names=1, sep=",")



d <- read.csv("spear_cor_Core_GABA_mean_10_37_109.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("spear_cor_Core_GABA_Low_B73.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("spear_cor_Core_LGABA_dimboa_H_B73.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("spear_cor_Core_LGABA_dimboa_L_6_84.csv", header=TRUE,row.names=1, sep=",")
d <- read.csv("fungi_bac_family_speaman_p_value.csv", header=TRUE,row.names=1, sep=",")





d <- read.csv("spear_cor_Core_dimboa_H.csv", header=TRUE,row.names=1, sep=",")
d<- read.csv("spear_cor_Core_dimboa_L.csv", header=TRUE,row.names=1, sep=",")

d <- read.csv("spear_cor_Core_HGABA_Hdimboa_18.csv", header=TRUE,row.names=1, sep=",")
d<- read.csv("spear_cor_Core_HGABA_Hdimboa_100.csv", header=TRUE,row.names=1, sep=",")


d <- read.csv("fungi_bac_family_speaman_p_value.csv", header=TRUE,row.names=1, sep=",")
###for only the p value
d[d == 1] <- NA
d=na.omit(d)
nrow(d)

#old with p value
##stacknodes <- read.csv("spear_DIMBOA_GABA_DEG_all_degree_0.7_p_less_equal_0.01.csv", header=TRUE,row.names=1, sep=",")
###no p value
stacknodes <- read.csv("spear_cor_Core_DIMBOA_fungi_bac_0.9.csv", header=TRUE,row.names=1, sep=",")
stacknodes <- read.csv("spear_cor_Core_DIMBOA_fungi_bac_0.7_001.csv", header=TRUE,row.names=1, sep=",")


#spear_cor_Core_GABA_mean_10_37_109_degree_0.7
#spear_cor_Core_GABA_Low_degree_0.7
#spear_cor_Core_DIMBOA_H_degree_0.7.csv
#spear_cor_Core_LGABA_dimboa_H_B73_degree_0.7

#GABA-2020
#output the degree--ONLY NEED TO CALCULATE THE CONNECTION 
graph_cors <- d %>%  
  filter(abs(r) > .5) %>%  # filter() out any correlations with an absolute value less than some threshold.
  #filter(abs(p) == 0) %>%
  filter(abs(p) <= 0.001) %>%
  graph_from_data_frame(directed = FALSE)
#filter(abs(value) > .9) %>% #no pvalue
#filter(p <= .01) %>%   ##this is for DIMBOA and GABA

x= degree(graph_cors)
write.csv(x, file="spear_cor_High_DIMBOA_Root.csv")

write.csv(x, file="spear_cor_Core_DIMBOA_fungi_bac_0.6_001.csv")


##transgenic sorghum
#output the degree--ONLY NEED TO CALCULATE THE CONNECTION 
graph_cors <- d %>%  
  filter(r > .8) %>%  # filter() out any correlations with an absolute value less than some threshold.
  #filter(abs(p) <= .01) %>%
  filter(p <= .01) %>%   ##this is for DIMBOA and GABA
  graph_from_data_frame(directed = FALSE)
x= degree(graph_cors)
write.csv(x, file="spear_WT-RHZ-AUG-otu_sp123_0.8_posi.csv")


###fungi bacteria
setwd("~/Desktop/Kaiju_tax/Family_PHYLUM/statistical analysis/only_diff_for_network")
d <- read.csv("fungi_bac_family_speaman_p_value.csv", header=TRUE,row.names=1, sep=",")
d[d == 1] <- NA
d=na.omit(d)
nrow(d)

###input the nodes
stacknodes <- read.csv("spear_cor_High_DIMBOA_Root.csv", header=TRUE,row.names=1, sep=",")
stacknodes=subset(stacknodes, tax %in% c("Atheliaceae",	"Ceratobasidiaceae",	"Botryobasidiaceae",	"Nectriaceae",	"Gloniaceae",	"Aspergillaceae",	"Neocallimastigaceae",	"Herpotrichiellaceae",	"Mortierellaceae",	"Trichocomaceae",	"Myxotrichaceae",	"Xylariaceae",	"Melampsoraceae",	"Chaetomiaceae",	"Magnaporthaceae",	"Chlamydiaceae",	"Chitinophagaceae",	"Nitrospiraceae",	"Polyangiaceae",	"Flavobacteriaceae",	"Geodermatophilaceae",	"Bacillaceae",	"Archangiaceae",	"Cytophagaceae",	"Ktedonobacteraceae",	"Kofleriaceae",	"Thermomonosporaceae",	"Planctomycetaceae",	"Xanthomonadaceae",	"Gemmatimonadaceae"))
stacknodes=select(stacknodes,tax,connection)
#input the data for DIMBOA_bac_fungi
graph_cors <- d %>%
  filter(abs(r) > .7) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(abs(p) <= 0.001) %>%
  subset( tax2 %in% c("Atheliaceae",	"Ceratobasidiaceae",	"Botryobasidiaceae",	"Nectriaceae",	"Gloniaceae",	"Aspergillaceae",	"Neocallimastigaceae",	"Herpotrichiellaceae",	"Mortierellaceae",	"Trichocomaceae",	"Myxotrichaceae",	"Xylariaceae",	"Melampsoraceae",	"Chaetomiaceae",	"Magnaporthaceae",	"Chlamydiaceae",	"Chitinophagaceae",	"Nitrospiraceae",	"Polyangiaceae",	"Flavobacteriaceae",	"Geodermatophilaceae",	"Bacillaceae",	"Archangiaceae",	"Cytophagaceae",	"Ktedonobacteraceae",	"Kofleriaceae",	"Thermomonosporaceae",	"Planctomycetaceae",	"Xanthomonadaceae",	"Gemmatimonadaceae"))%>%  #replace the Alteromonadaceae to Thermomonosporaceae
  subset( tax1 %in% c("Atheliaceae",	"Ceratobasidiaceae",	"Botryobasidiaceae",	"Nectriaceae",	"Gloniaceae",	"Aspergillaceae",	"Neocallimastigaceae",	"Herpotrichiellaceae",	"Mortierellaceae",	"Trichocomaceae",	"Myxotrichaceae",	"Xylariaceae",	"Melampsoraceae",	"Chaetomiaceae",	"Magnaporthaceae",	"Chlamydiaceae",	"Chitinophagaceae",	"Nitrospiraceae",	"Polyangiaceae",	"Flavobacteriaceae",	"Geodermatophilaceae",	"Bacillaceae",	"Archangiaceae",	"Cytophagaceae",	"Ktedonobacteraceae",	"Kofleriaceae",	"Thermomonosporaceae",	"Planctomycetaceae",	"Xanthomonadaceae",	"Gemmatimonadaceae"))%>%
  graph_from_data_frame(vertices = stacknodes, directed = F)

#input the data for TE
graph_cors <- d %>%
  filter(abs(r) > .8) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(abs(p) <= .01) %>%
  graph_from_data_frame(vertices = stacknodes, directed = FALSE)

#input the data for DIMBOA GABA DEGs
graph_cors <- d %>%
  filter(abs(r) > .7) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(p <= .01) %>%
  graph_from_data_frame(vertices = stacknodes, directed = FALSE)


#input the data for DIMBOA GABA 2020 August root
### 2020 August redo redo
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_High_DIMBOA_Root")
setwd("~/Desktop/Havelock_2018/RHZ_2018_corn_analysis/NEW_NEW_ANALYSIS_07312020/Sfigure10/Network/pvals_Medium_DIMBOA_Root")

d <- read.csv("spear_cor_Core_High_DIMBOA_Root_2.csv", header=TRUE,row.names=1, sep=",")
graph_cors <- d %>%  
  filter(abs(r) > .6) %>%  # filter() out any correlations with an absolute value less than some threshold.
  #filter(abs(p) <= .01) %>%
  filter(p <= .001) %>%   ##this is for DIMBOA and GABA
  graph_from_data_frame(directed = FALSE)
x= degree(graph_cors)
write.csv(x, file="spear_nodes.csv")
stacknodes <- read.csv("spear_nodes.csv", header=TRUE,row.names=1, sep=",")

graph_cors <- d %>%
  filter(abs(r) > .6) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(p <= .001) %>%
  graph_from_data_frame(vertices = rownames(stacknodes), directed = FALSE)


#calculate the modularity
#https://kateto.net/networks-r-igraph
#High modularity for a partitioning reflects dense connections within communities and sparse connections across communities
WC <-edge.betweenness.community(graph_cors)  #jul-s-wt
WC2 <-edge.betweenness.community(graph_cors)


WC_100 <-edge.betweenness.community(graph_cors2)  #jul-s-wt
WC1 <-edge.betweenness.community(graph_cors) #aug-s-wt
WC2 <-edge.betweenness.community(graph_cors) #sep-s-wt
WC3 <-edge.betweenness.community(graph_cors) #jul-s-rnai
WC4 <-edge.betweenness.community(graph_cors) #aug-s-rnai
WC5 <-edge.betweenness.community(graph_cors) #sep-s-rnai
#-rhz
WC6 <-edge.betweenness.community(graph_cors)  #jul-r-wt
WC7 <-edge.betweenness.community(graph_cors) #aug-r-wt
WC8 <-edge.betweenness.community(graph_cors) #sep-r-wt
WC9 <-edge.betweenness.community(graph_cors) #jul-r-rnai
WC10 <-edge.betweenness.community(graph_cors) #aug-r-rnai
WC11 <-edge.betweenness.community(graph_cors) #sep-r-rnai
#
#2020-AUgust-root
modularity(WC)
modularity(WC2)

#hub
set.seed(1)
plot(graph_cors, directed=F,vertex.label=NA,layout=layout.kamada.kawai,edge.curved=0) #layout.fruchterman.reingold

##temporary
#node-betweeness
g.b <- betweenness(graph_cors, v = V(graph_cors), directed = F, weights = NULL,
                   nobigint = TRUE, normalized = FALSE)
melt(g.b)
bet <- melt(g.b)
write.csv(bet, file="betweenness-pos-100.csv")
#hist(g.b, breaks = 20)
table(g.b)
#node degree
g.outd<- degree(graph_cors, mode = c("out"))
#hist(g.outd, breaks = 20)
table(g.outd)
melt(g.outd)
degree <- melt(g.outd)
write.csv(degree, file="degree-pos_100.csv")
#######



#input the data
#graph_cors <- d %>%
# filter(abs(r) > .80) %>%  # filter() out any correlations with an absolute value less than some threshold.
#filter(abs(p) <= .01) %>%
#graph_from_data_frame(vertices = stacknodes, directed = FALSE)
#calculate the modularity
#WC <-edge.betweenness.community(graph_cors)  
#modularity(WC)
#Find the hubs
hub <- hub_score(graph_cors)$vector
hub_value <- sort(hub, decreasing = T)
hub_top<- names(sort(hub, decreasing = T)) #
hub_top[1:1]
hub_value_export <- melt(hub_value)
write.csv(hub_value_export, file="hub_value_export.csv")
#eaverage path
average.path.length(graph_cors, directed=F, unconnected=T)
mean_distance(graph_cors, directed=F, unconnected=T)
#diamter
diameter(graph_cors,directed = TRUE, unconnected = TRUE, weights = NULL)
get_diameter(graph_cors) 
#clustering coefficient
transitivity(graph_cors, type = "average")
#total links and nodes
gsize(graph_cors)
gorder(graph_cors)
#remove mudules less than 5 links
number <- cluster_louvain(graph_cors)
table(number$membership)

#to extract the modules
xxx=as.data.frame(number$membership)
DIMBOA_bac_fungi <- cbind(xxx,stacknodes)
write.csv(DIMBOA_bac_fungi, "file=module.csv")



Small = which(table(number$membership) < 5)
Keep = V(graph_cors)[!(number$membership %in% Small)]
## Get subgraph & plot
graph_cors_2 = induced_subgraph(graph_cors, Keep)
gsize(graph_cors_2)
gorder(graph_cors_2)
#R2 of law
fit_power_law(graph_cors)
#node-betweeness
g.b <- betweenness(graph_cors, v = V(graph_cors), directed = F, weights = NULL,
                   nobigint = TRUE, normalized = FALSE)
melt(g.b)
hist(g.b, breaks = 20)
table(g.b)
#node degree
g.outd<- degree(graph_cors, mode = c("out"))
#hist(g.outd, breaks = 20)
table(g.outd)

#random network (100 random runs)
g.big.er = erdos.renyi.game(104,2048, type=c("gnm"), directed =F, loops=F) # WIld type august
fit_power_law(g.big.er)
average.path.length(g.big.er, directed=F, unconnected=T)
transitivity(g.big.er, type = "average")
WC <-edge.betweenness.community(g.big.er)  
modularity(WC)
par(mfrow=c(10,10),oma = c(0,0,0,0), mar = c(0,0,0,0))
lst <- list() #create a list to store your results
for(i in 1:100) #this statement does the repetition (looping)
{
  g.big.er = erdos.renyi.game(57,248, type=c("gnm"), directed =F, loops=F)
  x <- average.path.length(g.big.er, directed=F, unconnected=T) ##x, y are vectors 
  y <- transitivity(g.big.er, type = "average")
  fit_power_law(g.big.er)
  z <- diameter(g.big.er,directed = TRUE, unconnected = TRUE, weights = NULL)
  print(x)  
  print(y) 
  print(z) 
  lst[i] <- x
  lst[i] <- y
  lst[i] <- z
}
str(lst)
output<- lst
write.csv(output, file = "output_100_times_random_network.csv")




#for (year in c(2010,2011,2012,2013,2014,2015)){
# print(paste("The year is", year))
#}

##


#Make the plots-I
set.seed(1)
E(graph_cors)$color <-rgb(0,0,0, alpha=0.2)
ego <- names(which.max(degree(graph_cors)))
ego_2<- names(sort(degree(graph_cors), decreasing = T)) #
V(graph_cors)[V(graph_cors)!=ego]$color='grey'
V(graph_cors)[ego]$color = 'red'
#V(graph_cors)[hub_top[1:1]]$color = c('red')
#V(graph_cors)[hub_top[1:5]]$color = c('red','blue','dark green','black','purple')
V(graph_cors)$size = degree(graph_cors)*0.2
par(mfrow=c(1,1),oma = c(0,4,0,0), mar = c(0,0,4,0))

set.seed(1)
plot(graph_cors, directed=F, vertex.label=NA,vertex.size=4, vertex.color=membership(WC),edge.width=0.5,layout=layout.graphopt(graph_cors,niter=300,spring.length=2,spring.constant=1.5),edge.curved=0) #layout.mds#layout.graphopt#layout.kamada.kawai #vertex.label=V(graph_cors)$name
##plot(graph_cors, directed=F, vertex.label=NA,vertex.color=membership(WC1),edge.width=0.5,layout=layout.reingold.tilford(graph_cors,circular=T),edge.curved=0) #layout.mds#layout.graphopt#layout.kamada.kawai #vertex.label=V(graph_cors)$name
#try edges
V(graph_cors)[V(graph_cors)!=ego]$color=membership(WC11)
V(graph_cors)[ego]$color = membership(WC11)
V(graph_cors)$color <- membership(WC11)
V(graph_cors)$width <- 0
#V(graph_cors)$frame.color <- "white"
set.seed(123)
plot(graph_cors, directed=F, vertex.label=NA, vertex.size=2,vertex.color=membership(WC),edge.color=membership(WC), edge.width=2,layout=layout.graphopt(graph_cors,niter=300,spring.length=2,spring.constant=1.5),edge.curved=0) #layout.mds#layout.graphopt#layout.kamada.kawai #vertex.label=V(graph_cors)$name



#Make the plots-II
set.seed(1)
set.seed(2)
E(graph_cors)$color <-rgb(0,0,0, alpha=0.2)
ego <- names(which.max(degree(graph_cors)))
ego_2<- names(sort(degree(graph_cors), decreasing = T)) #
V(graph_cors)[V(graph_cors)!=ego]$color='grey'
V(graph_cors)[ego]$color = 'red'
V(graph_cors)[hub_top[1:1]]$color = c('red')
#V(graph_cors)[hub_top[1:5]]$color = c('red','blue','dark green','black','purple')

#V(graph_cors)$size = degree(graph_cors)*0.1 #High GABA use 0.1, rest use 0.2

V(graph_cors)$size = degree(graph_cors)#for ROOT 
par(mfrow=c(1,1),oma = c(3,2,0,0),
    mar = c(0,0,0,0))
set.seed(1) # Low GABA
set.seed(2) #medium GABA
set.seed(3) #low DIMBOA
set.seed(4) #high DIMBOA
plot(graph_cors, directed=F,vertex.label=NA,vertex.size=4, vertex.color=membership(WC),layout=layout_with_fr,edge.width=1,edge.curved=0)
#re scale
l <- layout_with_fr(graph_cors)
l <- layout_with_kk(graph_cors)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(graph_cors, directed=F,vertex.label=NA,vertex.size=4, vertex.color=membership(WC),layout=layout_with_fr,edge.width=1,edge.curved=0,rescale=T, layout=l*1.5)

###TEST
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 

# Remove layouts that do not apply to our graph.

layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]



par(mfrow=c(4,4), mar=c(1,1,1,1))

for (layout in layouts) {
  
  print(layout)
  
  l <- do.call(layout, list(graph_cors)) 
  
  plot(graph_cors, edge.arrow.mode=0, layout=l, main=layout,vertex.label=NA,vertex.color=membership(WC),edge.curved=0) } #vertex.color=membership(WC),directed=F,
###TEST




###OKOK
#h-d-root
E(graph_cors)$color <-rgb(0,0,0, alpha=0.2)
ego <- names(which.max(degree(graph_cors)))
ego_2<- names(sort(degree(graph_cors), decreasing = T)) #
V(graph_cors)[V(graph_cors)!=ego]$color='grey'
V(graph_cors)[ego]$color = 'red'
V(graph_cors)[hub_top[1:1]]$color = c('red')

V(graph_cors)$size = degree(graph_cors) *0.5#for ROOT 
par(mfrow=c(1,1),oma = c(3,2,0,0),
    mar = c(0,0,0,0))

set.seed(2)
plot(graph_cors, directed=F,vertex.label=NA,vertex.color=membership(WC2),layout=layout.fruchterman.reingold,edge.width=1,edge.curved=0)

##grey network
plot(graph_cors, directed=F,vertex.label=NA,vertex.color="grey",vertex.size=1,layout=layout.fruchterman.reingold,edge.width=1,edge.curved=0)
###highlight the nodes
inc.edges <-incident(graph_cors,V(graph_cors)[name=="Otu147"], mode="all")
ecol <-rep("grey",ecount(graph_cors))
ecol[inc.edges] <- "blue"
ew <-rep(1,ecount(graph_cors))
ew[inc.edges] <- 2
vcol <-rep("grey",vcount(graph_cors))
vcol[V(graph_cors)$name=="Otu147"] <- "red"
#E(graph_cors)$weight <- edge.betweenness(graph_cors)
set.seed(123456)
plot(graph_cors, directed=F, vertex.size=2,edge.width=ew, vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0,vertex.label=NA,layout=layout.fruchterman.reingold,edge.curved=0) #edge.width=E(graph_cors)$width, 

#set.seed(123)
#plot(graph_cors, directed=F,vertex.label=NA,layout=layout.fruchterman.reingold,edge.curved=0)

#3D
coords <- layout_with_fr(graph_cors,dim=3)
rglplot(graph_cors, directed=F,vertex.label=NA,edge.curved=0,layout=coords)
graph_cors



##power of law
fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}
