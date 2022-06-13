# RaritanBay_HAB_Networks
## By: Megan B. Rothenberger, Samantha J. Gleich, and Evan Flint
### Last Updated: June 12, 2022
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
  
  # NetGAM
  gamOut <- netGAM.df(dfSite,as.numeric(monthz),MCount,clrt=FALSE)
  gamOut$X <- as.numeric(Xz)
  gamOut <- left_join(gamOut,environ)
  gamOut$X <- NULL
  gamOut$YEAR <- NULL
  
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
```
## Make network image plots for each site
```
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

pdf("Site1_Netwk.pdf",width=6,height=4)
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

pdf("Site2_Netwk.pdf",width=6,height=4)
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

pdf("Site3_Netwk.pdf",width=6,height=4)
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

pdf("Site4_Netwk.pdf",width=6,height=4)
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

pdf("Site5_Netwk",width=6,height=4)
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

pdf("Site6_Netwk.pdf",width=6,height=4)
plot(netwk6,vertex.label=NA,vertex.size=7,vertex.color=V(netwk6)$col,layout=layout_with_fr(netwk6),main="Site 6")
dev.off()
```
## Calculate network statistics
```
mean(degree(netwk1))
transitivity(netwk1, type="global")
mean_distance(netwk1)
```
## Calculate
```
edgez <- edge6
edgez <- as.data.frame(edgez)
sum(edgez$V1%in%environ|edgez$V2%in%environ)
write.csv(edgez,"Site6_Edgelist_Env.csv")

edgezHab <- subset(edgez,Var1%in%habz|Var2%in%habz)
fin <- NULL
fin2 <- NULL
for(species in 1:14){
  sub <- subset(edgezHab,Var1==habz[species]|Var2==habz[species])
  namez <- V(netwk6)$name
  namez <- namez[-c(109:117)]
  ind <- which(namez==habz[species])
  namez <- namez[-c(ind)]
  s1 <- sum(sub$Var1 %in% namez)
  s2 <- sum(sub$Var2 %in% namez)
  s3 <- sum(sub$Var1 %in% environ)
  s4 <- sum(sub$Var2 %in% environ)
  s5 <- sum(sub$edge>0)
  s6 <- sum(sub$edge<0)
  out <- c(s1+s2,s3+s4)
  fin <- cbind(fin,out)
  out2 <- c(s5,s6)
  fin2 <- cbind(fin2,out2)
}
fin <- as.data.frame(t(fin))
rownames(fin) <- habz
colnames(fin)<- c("Biotic","Abiotic")

fin2 <- as.data.frame(t(fin2))
rownames(fin2) <- habz
colnames(fin2)<- c("Positive","Negative")

fin$Species <-rownames(fin) 
finMelt <- melt(fin,id.vars=c("Species"))

bioabio <- finMelt%>%ggplot(aes(x=Species,y=value,group=Species,fill=variable))+geom_bar(stat="identity",color="black")+scale_fill_manual(name="Association",values=c("gray","gray20"),label=c("Biotic","Abiotic"))+theme_classic(base_size=12)+scale_x_discrete(breaks=c("ALEX","ASAN","DINOP","HAKASH","HROT","HTRI","KARLO","OSCIL","PLIKE","PMIC","PMIN","PSEUDO","SCRIP","SNOW"),labels=c(expression(italic("Alexandrium spp.")),expression(italic("Akashiwo sanguinea")),expression(italic("Dinophysis acuminata")),expression(italic("Heterosigma akashiwo")),expression(italic("Heterocapsa rotundata")),expression(italic("Heterocapsa triquetra")),expression(italic("Karlodinium spp.")),expression(italic("Oscillatoria spp.")),expression(italic("Pfiestieria-like")),expression(italic("Prorocentrum micans")),expression(italic("Prorocentrum minimum")),expression(italic("Pseudo-nitzschia spp.")),expression(italic("Scrippsiella trochoidea")),expression(italic("Snowella spp."))))+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylab("Number of Associations")+xlab("HAB Species")
# ggsave("Biotic_Abiotic_Site6.pdf",width=8,height=4)

fin2$Species <-rownames(fin2) 
finMelt2 <- melt(fin2,id.vars=c("Species"))

posneg <- finMelt2%>%ggplot(aes(x=Species,y=value,group=Species,fill=variable))+geom_bar(stat="identity",color="black")+scale_fill_manual(name="Association",values=c("dodgerblue","indianred"))+theme_classic(base_size=12)+scale_x_discrete(breaks=c("ALEX","ASAN","DINOP","HAKASH","HROT","HTRI","KARLO","OSCIL","PLIKE","PMIC","PMIN","PSEUDO","SCRIP","SNOW"),labels=c(expression(italic("Alexandrium spp.")),expression(italic("Akashiwo sanguinea")),expression(italic("Dinophysis acuminata")),expression(italic("Heterosigma akashiwo")),expression(italic("Heterocapsa rotundata")),expression(italic("Heterocapsa triquetra")),expression(italic("Karlodinium spp.")),expression(italic("Oscillatoria spp.")),expression(italic("Pfiestieria-like")),expression(italic("Prorocentrum micans")),expression(italic("Prorocentrum minimum")),expression(italic("Pseudo-nitzschia spp.")),expression(italic("Scrippsiella trochoidea")),expression(italic("Snowella spp."))))+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylab("Number of Associations")+xlab("HAB Species")+ylim(0,15)

ggarrange(bioabio,posneg,ncol=1,nrow=2)
ggsave("BiAbio_PosNeg.pdf",width=8,height=12)
```
## Network Stats
```
#### Deg dist + network stats ###
net <- read.csv(file.choose(),header=TRUE,row.names=NULL)
net <- graph_from_edgelist(as.matrix(net),directed=FALSE)

net <- netwk6

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
```
