setwd("/Users/peng/Desktop/github")
setwd("~/Desktop/github")
library(reshape2)
library(dplyr)
##############################################
#Data wrangling for the output from the SparCC
##############################################
d <- read.csv("cor_OTU_AUG_WT_RHZ.csv", header=TRUE,row.names=1, sep=",")
p <- read.csv("pvals.two_sided.csv", header=TRUE,row.names=1, sep=",")
s=melt(d)
p=melt(p)
sp1 <- cbind(s,p, by=c("variable"))
head(sp1)
name=names(d)
name= data.frame(mm = name)
nrow(name)
names1 <- data.frame(name,i=rep(1:nrow(name),ea=NROW(name)))
sp2 = cbind(sp1,names1)
head(sp2)
colnames(sp2)=c("tax1", "r","tax1_1","p","by","tax2","i")
sp <- sp2 %>%  
  select("tax1","tax2","r","p")
write.csv(sp, file="spear_cor_Core_AUG_RHZ.csv")


#################
#Network analysis
#################
#use BiocManager::install() to install the following R packages
library(ggplot2)
library(ggraph)
library(tidyr)
library(tidygraph)
library(igraph)
library(reshape)
library(dplyr)
library(tibble)
###OTU-AUG-RHZ-WT as an example
d <- read.csv("spear_cor_Core_AUG_RHZ.csv", header=TRUE,row.names=1, sep=",")

#output the degree--ONLY NEED TO CALCULATE THE CONNECTION 
graph_cors <- d %>%  
  filter(r > 0.8) %>%  # filter() out any correlations with an absolute value less than certain threshold.
  filter(p <= .01) %>%  
  graph_from_data_frame(directed = FALSE)

#calculate the degree for each OTUs after the filteration
stacknodes1 = degree(graph_cors)
write.csv(stacknodes1, file="spear_WT-RHZ-AUG-otu_0.8_posi.csv")

stacknodes <- data.frame(tax = stacknodes1) %>%
  rownames_to_column()

colnames(stacknodes) <-  c("tax","connection")

#only keep the nodes' degree whose p value and r value p<=0.01 and r>0.8
graph_cors <- d %>%
  filter(r > .8) %>%  # filter() out any correlations with an absolute value less than some threshold.
  filter(p <= .01) %>%
  graph_from_data_frame(vertices = stacknodes, directed = FALSE)


#calculate the modularity
#https://kateto.net/networks-r-igraph
#High modularity for a partitioning reflects dense connections within communities and sparse connections across communities
WC <-edge.betweenness.community(graph_cors) 
modularity(WC)

#Node-betweeness
g.b <- betweenness(graph_cors, v = V(graph_cors), directed = F, weights = NULL,
                   nobigint = TRUE, normalized = FALSE)
melt(g.b)
bet <- melt(g.b)
write.csv(bet, file="betweenness-pos.csv")
#hist(g.b, breaks = 20)
table(g.b)

#Node degree
g.outd<- degree(graph_cors, mode = c("out"))
#hist(g.outd, breaks = 20)
table(g.outd)
melt(g.outd)
degree <- melt(g.outd)
write.csv(degree, file="degree-pos.csv")

#Find the hubs
hub <- hub_score(graph_cors)$vector
hub_value <- sort(hub, decreasing = T)
hub_top<- names(sort(hub, decreasing = T)) #
hub_top[1:1] #to see the first hub OTU
hub_value_export <- melt(hub_value)
write.csv(hub_value_export, file="hub_value_export.csv")


#average path
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

##remove mudules less than 5 links
number <- cluster_louvain(graph_cors)
table(number$membership)
#to extract the modules
xxx=as.data.frame(number$membership)
RHZ_AUG_WT <- cbind(xxx,stacknodes)
write.csv(RHZ_AUG_WT, "file_module.csv")
Small = which(table(number$membership) < 5)
Keep = V(graph_cors)[!(number$membership %in% Small)]
#Get subgraph & plot
graph_cors_2 = induced_subgraph(graph_cors, Keep)
gsize(graph_cors_2)
gorder(graph_cors_2)


#R2 of law
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

fit_power_law(graph_cors)


#a loop for calculating random network with 100 random runs
g.big.er = erdos.renyi.game(910,3576, type=c("gnm"), directed =F, loops=F) # WIld type august
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


#Make the plots-I # same node size
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

#Make the plots-II # highlight edges
set.seed(1)
##plot(graph_cors, directed=F, vertex.label=NA,vertex.color=membership(WC1),edge.width=0.5,layout=layout.reingold.tilford(graph_cors,circular=T),edge.curved=0) #layout.mds#layout.graphopt#layout.kamada.kawai #vertex.label=V(graph_cors)$name
#try edges
V(graph_cors)[V(graph_cors)!=ego]$color=membership(WC)
V(graph_cors)[ego]$color = membership(WC)
V(graph_cors)$color <- membership(WC)
V(graph_cors)$width <- 0
#V(graph_cors)$frame.color <- "white"
plot(graph_cors, directed=F, vertex.label=NA, vertex.size=2,vertex.color=membership(WC),edge.color=membership(WC), edge.width=2,layout=layout.graphopt(graph_cors,niter=300,spring.length=2,spring.constant=1.5),edge.curved=0) #layout.mds#layout.graphopt#layout.kamada.kawai #vertex.label=V(graph_cors)$name



#Make the plots-III highlight the hub OTUs
set.seed(1)
E(graph_cors)$color <-rgb(0,0,0, alpha=0.2)
ego <- names(which.max(degree(graph_cors)))
ego_2<- names(sort(degree(graph_cors), decreasing = T)) #
V(graph_cors)[V(graph_cors)!=ego]$color='grey'
V(graph_cors)[ego]$color = 'red'
V(graph_cors)[hub_top[1:1]]$color = c('red')  #hight the hub
#V(graph_cors)[hub_top[1:5]]$color = c('red','blue','dark green','black','purple')
V(graph_cors)$size = degree(graph_cors)#
par(mfrow=c(1,1),oma = c(3,2,0,0),
    mar = c(0,0,0,0))
plot(graph_cors, directed=F,vertex.label=NA,vertex.size=4,layout=layout_with_fr,edge.width=1,edge.curved=0)


#Make the plots-IV  
#####from website, source: https://kateto.net/netscix2016.html
###
l <- layout_with_fr(graph_cors)
l <- layout_with_kk(graph_cors)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(graph_cors, directed=F,vertex.label=NA,vertex.size=4, vertex.color=membership(WC),layout=layout_with_fr,edge.width=1,edge.curved=0,rescale=T, layout=l*1.5)
###to test the different layout methods
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(4,4), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <- do.call(layout, list(graph_cors)) 
  plot(graph_cors, edge.arrow.mode=0, layout=l, main=layout,vertex.label=NA,vertex.color=membership(WC),edge.curved=0) } #vertex.color=membership(WC),directed=F,


####Make the plots-V 
##just grey network
plot(graph_cors, directed=F,vertex.label=NA,vertex.color="grey",vertex.size=1,layout=layout.fruchterman.reingold,edge.width=1,edge.curved=0)
###highlight the edges
inc.edges <-incident(graph_cors,V(graph_cors)[name=="Otu147"], mode="all")
ecol <-rep("grey",ecount(graph_cors))
ecol[inc.edges] <- "blue"
ew <-rep(1,ecount(graph_cors))
ew[inc.edges] <- 2
vcol <-rep("grey",vcount(graph_cors))
vcol[V(graph_cors)$name=="Otu147"] <- "red"
#E(graph_cors)$weight <- edge.betweenness(graph_cors)
plot(graph_cors, directed=F, vertex.size=2,edge.width=ew, vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0,vertex.label=NA,layout=layout.fruchterman.reingold,edge.curved=0) #edge.width=E(graph_cors)$width, 

###Make the plots-VI, actually used in the paper
E(graph_cors)$color <-rgb(0,0,0, alpha=0.2)
ego <- names(which.max(degree(graph_cors)))
ego_2<- names(sort(degree(graph_cors), decreasing = T)) #
V(graph_cors)[V(graph_cors)!=ego]$color='grey'
#V(graph_cors)[ego]$color = 'red'
#V(graph_cors)[hub_top[1:1]]$color = c('red')
V(graph_cors)$size = degree(graph_cors) *0.2  # can change to scale the node size
par(mfrow=c(1,1),oma = c(3,2,0,0),
    mar = c(0,0,0,0))
plot(graph_cors, directed=F,vertex.label=NA,vertex.color=membership(WC),layout=layout.fruchterman.reingold,edge.width=1,edge.curved=0)


#Make the plots-VII - 3D
coords <- layout_with_fr(graph_cors,dim=3)
rglplot(graph_cors, directed=F,vertex.label=NA,vertex.color=membership(WC),edge.curved=0,layout=coords)





