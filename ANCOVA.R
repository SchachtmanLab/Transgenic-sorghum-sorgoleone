
##network stat analysis
setwd("~/Desktop/R/ANCOVA")
dat <- read.csv("Degree_2016_TE_OTU_cooc.csv", header=TRUE,row.names=1, sep=",")
dat <- read.csv("Closeness_2016_TE_OTU_cooc.csv", header=TRUE,row.names=1, sep=",")
dat$Date <- factor(dat$Date)


Network_result <- aov(Degree_norm_log ~ Sample_type*Date*Genotype,data=dat)
print(summary(Network_result))

Network_result <- aov(Size~Sample_type*Time*Genotype,data=dat)
print(summary(Network_result))

Network_result <- aov(Size~Sample_type*Genotype,data=dat)
print(summary(Network_result))

Network_result <- aov(Connectivity~Sample_type*Genotype,data=dat)
print(summary(Network_result))

Network_result <- aov(Connectivity~Genotype,data=dat)
print(summary(Network_result))

Network_result <- aov(Connectivity~Genotype*Time,data=dat)
print(summary(Network_result))

###genotype and time interaction
Network_result <- aov(Clustering_coefficient~Genotype*Time,data=dat)
print(summary(Network_result))

Network_result <- aov(Closeness_centrality~Genotype*Time,data=dat)
print(summary(Network_result))

