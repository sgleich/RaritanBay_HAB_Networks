### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 15, 2022 ###
### Figure 4 - Networks at each site with HAB species color coded ###

# Libraries
library(mgcv)
library(tidyverse)
library(reshape2)
library(igraph)
library(huge)
library(pulsar)
library(NetGAM)
library(missForest)
library(psych)


### Figure 3 - Networks at each site ###
# Load data
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfZoo$TOTAL <- NULL
dfZoo$X.1 <- NULL
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)

environ <- read.csv(file="Environ.csv",header=TRUE)
environ$variable <- rownames(environ)
environ <- environ[environ$YEAR != 2022,]
# environ6 <- environ[environ$SITE == 6,]
environ <- environ[environ$DEPTH == 1,]
environ$CONDUCTIVITY..mS.cm. <- NULL
environ$DO....SAT.<- NULL
environ$variable <- NULL

colz <- colsplit(environ$DATE,"/",c("month","day","year"))
environ$X <- paste(colz$month,environ$YEAR,environ$SITE,sep="")
environ$DATE <- NULL
environ$YEAR <- NULL
environ$SEASON <- NULL
environ$DEPTH <- NULL
environ$SITE <- NULL

colnames(environ) <- c("TEMP","pH","SAL","DO","Secchi","NO3","NH4","SRP","Si","X")
environ <- mutate_all(environ, function(x) as.numeric(as.character(x)))

environImp <- missForest::missForest(environ)
environImp <- environImp$ximp

dfTotal$SEAS <- NULL
dfTotal$DATE <- NULL
dfTotal$SEASON <- NULL
dfTotal$DEPTH <- NULL
dfTotal$SITE <- as.numeric(dfTotal$SITE)

# Make network
makeNet <- function(df,site){
  dfSite <- subset(df,SITE==site)
  dfSite$X <- as.character(dfSite$X)
  monthz <- substr(dfSite$X,1,2)
  monthz <- ifelse(monthz>12,substr(monthz,1,1),monthz)
  dfSite$MOY <- as.numeric(monthz)
  MCount <- c((dfSite$MOY - 3) + 12*(dfSite$YEAR - 1))
  dfSite$MOY <- NULL
  Xz <- dfSite$X
  dfSite <- dfSite[,3:ncol(dfSite)]
  dfSite <- as.data.frame(t(dfSite))
  dfSite <- subset(dfSite,rowSums(dfSite)!=0)
  dfSite <- as.data.frame(t(dfSite))
  dfSite$YEAR <- NULL
  
  # NetGAM
  gamOut <- netGAM.df(dfSite,as.numeric(monthz),MCount,clrt=FALSE)
  gamOut$X <- as.numeric(Xz)
  gamOut <- left_join(gamOut,environ)
  gamOut$X <- NULL
  
  # Glasso
  gamOut <- huge.npn(gamOut)
  lams  <- getLamPath(getMaxCov(as.matrix(gamOut)), .01, len=30)
  hugeargs <- list(lambda=lams, verbose=FALSE)
  netOut <- batch.pulsar(gamOut, fun=huge::huge, fargs=hugeargs,rep.num=50, criterion='stars')
  fit <- refit(netOut)$refit
  fit.fin <- fit$stars
  fit.fin <- as.matrix(fit.fin)
  colnames(fit.fin)<- colnames(gamOut)
  rownames(fit.fin)<- colnames(gamOut)
  cor <- corr.test(gamOut,method="spearman")
  cor <- cor$r
  fitCor <- fit.fin*cor
  return(fitCor)}

net1 <- makeNet(dfTotal,1)
net2 <- makeNet(dfTotal,2)
net3 <- makeNet(dfTotal,3)
net4 <- makeNet(dfTotal,4)
net5 <- makeNet(dfTotal,5)
net6 <- makeNet(dfTotal,6)

# Edgelists
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

edge1 <- net1  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)
edge2 <- net2  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)
edge3 <- net3  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)
edge4 <- net4  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)
edge5 <- net5  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)
edge6 <- net6  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(edge = value)%>%filter(edge!=0)

mat1 <- edge1
mat1$edge <- NULL
mat1 <- as.matrix(mat1)

mat2 <- edge2
mat2$edge <- NULL
mat2 <- as.matrix(mat2)

mat3 <- edge3
mat3$edge <- NULL
mat3 <- as.matrix(mat3)

mat4 <- edge4
mat4$edge <- NULL
mat4 <- as.matrix(mat4)

mat5 <- edge5
mat5$edge <- NULL
mat5 <- as.matrix(mat5)

mat6 <- edge6
mat6$edge <- NULL
mat6 <- as.matrix(mat6)

netwk1 <- graph_from_edgelist(mat1,directed=FALSE)
netwk2 <- graph_from_edgelist(mat2,directed=FALSE)
netwk3 <- graph_from_edgelist(mat3,directed=FALSE)
netwk4 <- graph_from_edgelist(mat4,directed=FALSE)
netwk5 <- graph_from_edgelist(mat5,directed=FALSE)
netwk6 <- graph_from_edgelist(mat6,directed=FALSE)

commEdges <- intersection(netwk1,netwk2,netwk3,netwk4,netwk5,netwk6)
E(commEdges)

# Color by HABZ
habz <- c("PSEUDO","SNOW","SCRIP","PMIC","PMIN","PLIKE","OSCIL","KARLO","HTRI","HROT","HAKASH","DINOP","ASAN","ALEX")
environ<- c("TEMP","pH","SAL","DO","Secchi","NO3","NH4","SRP" ,"Si")

# Network 1
# Node colors
V(netwk1)$col<- "grey"
V(netwk1)$col <- ifelse(V(netwk1)$name=="PSEUDO","forestgreen",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="SNOW","dodgerblue",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="SCRIP","darkgoldenrod3",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="PMIC","orchid",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="PMIN","hotpink1",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="PLIKE","coral1",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="OSCIL","cyan",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="KARLO","darkolivegreen1",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="HTRI","indianred",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="HROT","khaki1",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="HAKASH","aquamarine",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="DINOP","cornflowerblue",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="ASAN","firebrick1",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name=="ALEX","hotpink4",V(netwk1)$col)
V(netwk1)$col <- ifelse(V(netwk1)$name %in% environ,"grey20",V(netwk1)$col)

pdf("../../Site1_Netwk.pdf",width=6,height=4)
plot(netwk1,vertex.label=NA,vertex.size=7,vertex.color=V(netwk1)$col,layout=layout_with_fr(netwk1),main="Site 1")
dev.off()

# Network 2
# Node colors
V(netwk2)$col <- "grey"
V(netwk2)$col <- ifelse(V(netwk2)$name=="PSEUDO","forestgreen",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="SNOW","dodgerblue",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="SCRIP","darkgoldenrod3",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="PMIC","orchid",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="PMIN","hotpink1",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="PLIKE","coral1",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="OSCIL","cyan",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="KARLO","darkolivegreen1",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="HTRI","indianred",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="HROT","khaki1",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="HAKASH","aquamarine",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="DINOP","cornflowerblue",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="ASAN","firebrick1",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name=="ALEX","hotpink4",V(netwk2)$col)
V(netwk2)$col <- ifelse(V(netwk2)$name %in% environ,"grey20",V(netwk2)$col)

pdf("../../Site2_Netwk.pdf",width=6,height=4)
plot(netwk2,vertex.label=NA,vertex.size=7,vertex.color=V(netwk2)$col,layout=layout_with_fr(netwk2),main="Site 2")
dev.off()

# Network 3
# Nodes
V(netwk3)$col <- "grey"
V(netwk3)$col <- ifelse(V(netwk3)$name=="PSEUDO","forestgreen",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="SNOW","dodgerblue",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="SCRIP","darkgoldenrod3",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="PMIC","orchid",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="PMIN","hotpink1",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="PLIKE","coral1",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="OSCIL","cyan",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="KARLO","darkolivegreen1",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="HTRI","indianred",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="HROT","khaki1",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="HAKASH","aquamarine",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="DINOP","cornflowerblue",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="ASAN","firebrick1",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name=="ALEX","hotpink4",V(netwk3)$col)
V(netwk3)$col <- ifelse(V(netwk3)$name %in% environ,"grey20",V(netwk3)$col)

pdf("../../Site3_Netwk.pdf",width=6,height=4)
plot(netwk3,vertex.label=NA,vertex.size=7,vertex.color=V(netwk3)$col,layout=layout_with_fr(netwk3),main="Site 3")
dev.off()

# Network 4
# Nodes
V(netwk4)$col <- "grey"
V(netwk4)$col <- ifelse(V(netwk4)$name=="PSEUDO","forestgreen",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="SNOW","dodgerblue",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="SCRIP","darkgoldenrod3",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="PMIC","orchid",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="PMIN","hotpink1",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="PLIKE","coral1",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="OSCIL","cyan",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="KARLO","darkolivegreen1",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="HTRI","indianred",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="HROT","khaki1",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="HAKASH","aquamarine",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="DINOP","cornflowerblue",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="ASAN","firebrick1",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name=="ALEX","hotpink4",V(netwk4)$col)
V(netwk4)$col <- ifelse(V(netwk4)$name %in% environ,"grey20",V(netwk4)$col)

pdf("../../Site4_Netwk.pdf",width=6,height=4)
plot(netwk4,vertex.label=NA,vertex.size=7,vertex.color=V(netwk4)$col,layout=layout_with_fr(netwk4),main="Site 4")
dev.off()

# Network 5
# Nodes
V(netwk5)$col <- "grey"
V(netwk5)$col <- ifelse(V(netwk5)$name=="PSEUDO","forestgreen",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="SNOW","dodgerblue",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="SCRIP","darkgoldenrod3",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="PMIC","orchid",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="PMIN","hotpink1",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="PLIKE","coral1",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="OSCIL","cyan",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="KARLO","darkolivegreen1",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="HTRI","indianred",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="HROT","khaki1",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="HAKASH","aquamarine",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="DINOP","cornflowerblue",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="ASAN","firebrick1",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name=="ALEX","hotpink4",V(netwk5)$col)
V(netwk5)$col <- ifelse(V(netwk5)$name %in% environ,"grey20",V(netwk5)$col)

pdf("../../Site5_Netwk",width=6,height=4)
plot(netwk5,vertex.label=NA,vertex.size=7,vertex.color=V(netwk5)$col,layout=layout_with_fr(netwk5),main="Site 5")
dev.off()

# Network 6
V(netwk6)$col <- "grey"
V(netwk6)$col <- ifelse(V(netwk6)$name=="PSEUDO","forestgreen",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="SNOW","dodgerblue",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="SCRIP","darkgoldenrod3",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="PMIC","orchid",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="PMIN","hotpink1",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="PLIKE","coral1",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="OSCIL","cyan",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="KARLO","darkolivegreen1",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="HTRI","indianred",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="HROT","khaki1",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="HAKASH","aquamarine",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="DINOP","cornflowerblue",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="ASAN","firebrick1",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name=="ALEX","hotpink4",V(netwk6)$col)
V(netwk6)$col <- ifelse(V(netwk6)$name %in% environ,"grey20",V(netwk6)$col)

pdf("../../Site6_Netwk.pdf",width=6,height=4)
plot(netwk6,vertex.label=NA,vertex.size=7,vertex.color=V(netwk6)$col,layout=layout_with_fr(netwk6),main="Site 6")
dev.off()
