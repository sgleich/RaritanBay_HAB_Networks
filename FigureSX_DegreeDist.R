### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 15, 2022 ###
### Figure S? - Degree distribution plots ###

# Libraries 
library(tidyverse)
library(reshape2)
library(igraph)
library(ggplot2)
library(ggpubr)

# Choose network
net <- netwk6

# Calculate information about network to parameterize random networks
nodez <-length(V(net)) # Number of nodes in network
edgez <- length(E(net)) # Number of edges in network
k <- nodez*(nodez-1)/2 # Number of possible edges in the network
edgeProb <- edgez/k

# Make 100 random networks and calculate degree distribution
outDeg <- NULL
outStats <- NULL
for (i in 1:100){
  randNet<- erdos.renyi.game(nodez,edgeProb)
  randDeg <- degree(randNet)
  outDeg <- cbind(outDeg,randDeg)
  meanDeg <- mean(randDeg)
  meanClust <- transitivity(randNet, type="global")
  meanPath <- mean_distance(randNet)
  vec <- c(meanDeg,meanClust,meanPath)
  outStats <- cbind(outStats,vec)
}

rownames(outStats) <- c("meanDeg","meanClust","meanPath")

outDeg <- as.data.frame(outDeg)
for (i in 1:length(outDeg)){
  colnames(outDeg)[i] <- paste(colnames(outDeg)[i],i,sep="_")
  colnames(outStats)[i] <- paste(colnames(outStats)[i],i,sep="_")
}

outDegMelt <- melt(outDeg)
outDegFrac <- outDegMelt %>%group_by(variable,value)%>%tally()
outDegTotal <- outDegMelt %>% group_by(variable)%>%summarize(tot=sum(value))
outDegFrac <- left_join(outDegFrac,outDegTotal)
outDegFrac$prob <- outDegFrac$n/outDegFrac$tot

# Now take the degree distribution vector of your network and make dataframe with it (vector is called 'deg' here)
deg <- degree(net)
sumz <- sum(deg)
myDeg <- data.frame(variable=rep("my_deg",length(deg)),value=deg,tot=rep(sumz,length(deg)))
myDeg <- myDeg%>%group_by(variable,value,tot)%>%tally()
myDeg$prob <- myDeg$n/myDeg$tot
myDeg <- myDeg[c(1:2,4,3,5)]

outDegFrac <- rbind(myDeg,outDegFrac)

degDist6 <- ggplot(outDegFrac,aes(x=value, y=prob))+geom_line(aes(group=variable),size=0.05)+theme_classic()+xlab("Degree")+ylab("Probability")+geom_line(data = outDegFrac[outDegFrac$variable == 'my_deg',], aes(group=variable, size=1), size=0.8,linetype="dashed",color="red")+ggtitle("Degree Distribution\nSite 6")

degDist1
degDist2
degDist3
degDist4
degDist5
degDist6

ggarrange(degDist1,degDist2,degDist3,degDist4,degDist5,degDist6,ncol=3,nrow=2)
ggsave("DegDist.pdf",width=12,height=7)

# Network stats and random network stats
sort(deg)
length(V(net))
length(E(net))
mean(deg)
transitivity(net, type="global")
mean_distance(net)
outStats <- as.data.frame(t(outStats))
mean(outStats$meanDeg)
sd(outStats$meanDeg)
mean(outStats$meanClust)
sd(outStats$meanClust)
mean(outStats$meanPath)
sd(outStats$meanPath)
