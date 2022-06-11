# RaritanBay_HAB_Networks
## By: Megan B. Rothenberger, Samantha J. Gleich, and Evan Flint
### Last Updated: June 11, 2022
Network analysis of the Rarity Bay Estuary phytoplankton community composition time-series dataset. See Rothenberger et al. (2014; https://doi.org/10.1007/s12237-013-9714-0) and  Rothenberger and Calomeni (2016; https://doi.org/10.1016/j.jembe.2016.03.015) for more information.

## Load packages
```
library(formattable)
library(htmltools)
library(webshot)  
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggarrange)
library(randomcoloR)
library(NetGAM)
library(mgcv)
library(igraph)
library(huge)
library(batchtools)
library(pulsar)
library(psych)
```

## Load in data and wrangle
```
### Load in data and wrangle ###
# Load data
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfZoo$TOTAL <- NULL
dfZoo$X.1 <- NULL
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)

# Month info
dfTotal$X <- as.character(dfTotal$X)
monthz <- substr(dfTotal$X,1,2)
monthz <- ifelse(monthz>12,substr(monthz,1,1),monthz)
dfTotal$month <- monthz
dfTotal$SEAS <- NULL
dfTotal$YEAR <- NULL
```
## Make taxa barplots for each site
```
### Make taxa barplots for each site ###
dfMelt <- melt(dfTotal,id.vars=c("SITE","month","X"))
dfCat <- read.csv("RB_Code.csv",header=TRUE)
colnames(dfCat) <- c("variable","Group")
dfMelt <- left_join(dfMelt,dfCat)

# Wrangle + Summarize
dfSumz <- dfMelt%>%group_by(Group,month,SITE,X)%>%summarize(s=sum(value))
ochro <- c("Tribophyte","Raphidophyte","Chrysophyte")
dfSumz$Group <- as.character(dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group%in%ochro, "Other Ochrophytes",dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group=="Euglenoid", "Euglenozoa",dfSumz$Group)
dfMeanz <- dfSumz%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()

# Set up plot
colrs <- randomcoloR::distinctColorPalette(length(unique(dfMeanz$Group)))
dfMeanz$month <- factor(dfMeanz$month,levels=c("4","5","6","7","8","9","10","11"))
dfMeanz$month <- as.numeric(as.character(dfMeanz$month))
dfMeanz <- subset(dfMeanz,!is.na(Group))

sitePlot <- function(df,ste,title){
  subs <- dfMeanz %>% filter(SITE==ste)%>%as.data.frame()
  subs$m <- ifelse(is.na(subs$m),0,subs$m)
  p1 <- ggplot(subs,aes(x=month, y=m,fill=Group))+geom_area(position="fill",color="#525252")+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic()+xlab("Month")+ylab("Relative Abundance")+scale_x_continuous(breaks=c(4,5,6,7,8,9,10,11),labels=c("April","May","June","July","August","September","October","November"))+ggtitle(paste(title))+theme(axis.text.x = element_text(angle = 45,hjust=1))
  return(p1)}

ste1 <- sitePlot(dfMeanz,1,"Site 1")
ste1
ste2 <- sitePlot(dfMeanz,2,"Site 2")
ste2
ste3 <- sitePlot(dfMeanz,3,"Site 3")
ste3
ste4 <- sitePlot(dfMeanz,4,"Site 4")
ste4
ste5 <- sitePlot(dfMeanz,5,"Site 5")
ste5
ste6 <- sitePlot(dfMeanz,6,"Site 6")
ste6

# ggarrange(ste1,ste2,ste3,ste4,ste5,ste6,ncol=3,nrow=2,common.legend = TRUE)
# ggsave("../../Figs/Figure2.pdf",width=12,height=8)
```
## Make taxa barplots pre-post Hurricane Sandy
```
# Load data
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfZoo$TOTAL <- NULL
dfZoo$X.1 <- NULL
dfZoo$TOTAL <- NULL
dfZoo$X.1 <- NULL
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)

# Month info
dfTotal$X <- as.character(dfTotal$X)
monthz <- substr(dfTotal$X,1,2)
monthz <- ifelse(monthz>12,substr(monthz,1,1),monthz)
dfTotal$month <- monthz
dfTotal$SEAS <- NULL

# Melt data
dfMelt <- melt(dfTotal,id.vars=c("SITE","month","X","YEAR"))
dfCat <- read.csv("RB_Code.csv",header=TRUE)
colnames(dfCat) <- c("variable","Group")
dfMelt <- left_join(dfMelt,dfCat)

# Wrangle + Summarize
dfSumz <- dfMelt%>%group_by(Group,month,SITE,X,YEAR)%>%summarize(s=sum(value))
ochro <- c("Tribophyte","Raphidophyte","Chrysophyte")
dfSumz$Group <- as.character(dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group%in%ochro, "Other Ochrophytes",dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group=="Euglenoid", "Euglenozoa",dfSumz$Group)

dfMeanz1 <- dfSumz%>%filter(YEAR==1|YEAR==2|YEAR==3)%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()
dfMeanz2 <- dfSumz%>%filter(YEAR==4)%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()
dfMeanz3 <- dfSumz%>%filter(YEAR==5|YEAR==6)%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()
dfMeanz4 <- dfSumz%>%filter(YEAR>6)%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()

# Set up plot
# colrs <- randomcoloR::distinctColorPalette(length(unique(dfMeanz$Group)))

plotYear <- function(dfMeanz,title){
  dfMeanz$month <- factor(dfMeanz$month,levels=c("4","5","6","7","8","9","10","11"))
  dfMeanz$month <- as.numeric(as.character(dfMeanz$month))
  dfMeanz <- subset(dfMeanz,!is.na(Group))
  siteMeanz <- subset(dfMeanz,SITE==6)
  siteMeanz$m <- ifelse(is.na(siteMeanz$m),0,siteMeanz$m)
  p1 <- ggplot(siteMeanz,aes(x=month, y=m,fill=Group))+geom_area(position="fill",color="#525252")+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic()+xlab("Month")+ylab("Relative Abundance")+scale_x_continuous(breaks=c(4,5,6,7,8,9,10,11),labels=c("April","May","June","July","August","September","October","November"))+ggtitle(paste(title))+theme(axis.text.x = element_text(angle = 45,hjust=1))
  return(p1)}

y1 <- plotYear(dfMeanz1,"Before Sandy\n2010-2012\nSite 6")
y2 <- plotYear(dfMeanz2,"Directly After Sandy\n2013\nSite 6")
y3 <- plotYear(dfMeanz3,"2-3 Years After Sandy\n2014-2015\nSite 6")
y4 <- plotYear(dfMeanz4,"4-10 Years After Sandy\n2016-2021\nSite 6")

# ggarrange(y1,y2,y3,y4,ncol=2,nrow=2,common.legend = TRUE)
# ggsave("../../Figs/Figure4.pdf",width=8,height=4)
```
## Make a network at each site
```
# Load data
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfZoo$TOTAL <- NULL
dfZoo$X.1 <- NULL
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)

# Make network function
makeNet <- function(df,site){
  dfSite <- subset(df,SITE==site)
  dfSite$X <- as.character(dfSite$X)
  monthz <- substr(dfSite$X,1,2)
  monthz <- ifelse(monthz>12,substr(monthz,1,1),monthz)
  dfSite$MOY <- as.numeric(monthz)
  MCount <- c((dfSite$MOY - 3) + 12*(dfSite$YEAR - 1))
  dfSite$MOY <- NULL
  dfSite <- dfSite[,5:ncol(dfSite)]
  dfSite <- as.data.frame(t(dfSite))
  dfSite <- subset(dfSite,rowSums(dfSite)!=0)
  dfSite <- as.data.frame(t(dfSite))
  
  # NetGAM - remove temporal signal
  gamOut <- netGAM.df(dfSite,as.numeric(monthz),MCount,clrt=FALSE)
  
  # Glasso - graphical lasso network analysis
  gamOut <- huge.npn(gamOut)
  lams  <- getLamPath(getMaxCov(as.matrix(gamOut)), .01, len=30)
  hugeargs <- list(lambda=lams, verbose=FALSE)
  netOut <- batch.pulsar(gamOut, fun=huge::huge, fargs=hugeargs,rep.num=50, criterion='stars')
  fit <- refit(netOut)$refit
  fit.fin <- fit$stars
  fit.fin <- as.matrix(fit.fin)
  colnames(fit.fin)<- colnames(gamOut)
  rownames(fit.fin)<- colnames(gamOut)
  
  # SCC - spearman correlation coefficient to use as edge weights
  cor <- corr.test(gamOut,method="spearman")
  cor <- cor$r
  fitCor <- fit.fin*cor
  return(fitCor)}

# Make networks for each site
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

# Prepare data for igraph
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

# Make igraph network objects
netwk1 <- graph_from_edgelist(mat1,directed=FALSE)
netwk2 <- graph_from_edgelist(mat2,directed=FALSE)
netwk3 <- graph_from_edgelist(mat3,directed=FALSE)
netwk4 <- graph_from_edgelist(mat4,directed=FALSE)
netwk5 <- graph_from_edgelist(mat5,directed=FALSE)
netwk6 <- graph_from_edgelist(mat6,directed=FALSE)
