# RaritanBay_HAB_Networks
## By: Megan B. Rothenberger, Samantha J. Gleich, and Evan Flint
### Last Updated: June 11, 2022
Network analysis of the Rarity Bay Estuary phytoplankton community composition time-series dataset. See Rothenberger et al. (2014; https://doi.org/10.1007/s12237-013-9714-0) and  Rothenberger and Calomeni (2016; https://doi.org/10.1016/j.jembe.2016.03.015) for more information.

## Load in data and wrangle
```
### Load in data and wrangle ###
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Load data
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)
dfTotal$X <- as.character(dfTotal$X)

# Get month vector
try <- substr(dfTotal$X,1,2)
try <- ifelse(try>12,substr(try,1,1),try)
dfTotal$month <- try
dfTotal$SEAS <- NULL
dfTotal$YEAR <- NULL

# Melt dataframe and add taxonomic group info (i.e. chlorophyte, other ochrophyte,cryptomonad, cyanobacteria, diatom, dinoflagellate, euglenozoa, haptophyte,mesozooplankton, microzooplankton)
dfMelt <- melt(dfTotal,id.vars=c("SITE","month","X"))
dfCat <- read.csv("RB_Code.csv",header=TRUE)
colnames(dfCat) <- c("variable","Group")
dfMelt <- left_join(dfMelt,dfCat)
dfMelt <- subset(dfMelt,!is.na(Group))
```
## Make taxa barplots for each site
```
### Make taxa barplots ###
# Sum the number of each taxa in each sample
dfSumz <- dfMelt%>%group_by(Group,month,SITE,X)%>%summarize(s=sum(value))

# Change names of some groups
ochro <- c("Tribophyte","Raphidophyte","Chrysophyte")
dfSumz$Group <- as.character(dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group%in%ochro, "Other Ochrophytes",dfSumz$Group)
dfSumz$Group <- ifelse(dfSumz$Group=="Euglenoid", "Euglenozoa",dfSumz$Group)

# Calculate the mean abundance of each taxonomic group at each month and site
dfMeanz <- dfSumz%>%group_by(month,SITE,Group)%>%summarize(m=mean(s))%>%as.data.frame()
dfMeanz <- subset(dfMeanz, !is.na(Group))

# Generate a random color palette
colrs <- randomcoloR::distinctColorPalette(length(unique(dfMeanz$Group)))

# Set up x-axis information
dfMeanz$month <- factor(dfMeanz$month,levels=c("4","5","6","7","8","9","10","11"))
dfMeanz$month <- as.numeric(as.character(dfMeanz$month))

# Plot (change SITE for different sites)
p <- dfMeanz %>% filter(SITE==6)%>%ggplot(aes(x=month, y=m,fill=Group))+geom_area(position="fill",color="#525252")+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic()+xlab("Month")+ylab("Relative Abundance")+scale_x_continuous(breaks=c(4,5,6,7,8,9,10,11),labels=c("April","May","June","July","August","September","October","November"))+ggtitle("Site 6")+theme(axis.text.x = element_text(angle = 45,hjust=1))


# ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2,common.legend = TRUE)
# ggsave("../Fig2_Prelim.pdf",width=10,height=6)
```
## Make taxa barplots pre-post Hurricane Sandy
```
# Load data
dfZoo <- read.csv("RB_Zoo.csv",header=TRUE)
dfPhyto <- read.csv("RB_Phyto.csv",header=TRUE)
dfZoo2 <- mutate_all(dfZoo, function(x) as.numeric(as.character(x)))
dfPhyto2 <- mutate_all(dfPhyto, function(x) as.numeric(as.character(x)))

# Combine phytoplankton and zooplankton data
dfTotal <- left_join(dfZoo2,dfPhyto2)
dfTotal$X <- as.character(dfTotal$X)

# Get month vector
try <- substr(dfTotal$X,1,2)
try <- ifelse(try>12,substr(try,1,1),try)
dfTotal$month <- try
dfTotal$SEAS <- NULL

# Subset site 6
dfTotal6 <- subset(dfTotal,SITE==6)
dfTotal6$X <- NULL
dfTotal6$SITE <- NULL

# Wrangle data
dfMelt <- melt(dfTotal6,id.vars=c("YEAR","month"))
dfCat <- read.csv("RB_Code.csv",header=TRUE)
colnames(dfCat) <- c("variable","Group")
dfMelt <- left_join(dfMelt,dfCat)
dfMelt <- subset(dfMelt,!is.na(Group))

# Find total abundance of each taxonomic group for each month and year combination
dfSumz <- dfMelt%>%group_by(Group,month,YEAR)%>%summarize(s=sum(value))

# Set up the pre-post Sandy groupings
dfSumz$period <- ifelse(dfSumz$YEAR==1 |dfSumz$YEAR==2|dfSumz$YEAR==3,"Before Sandy",NA)
dfSumz$period <- ifelse(dfSumz$YEAR==4,"Directly After Sandy",dfSumz$period)
dfSumz$period <- ifelse(dfSumz$YEAR==5 |dfSumz$YEAR==6,"2-3 Years After Sandy",dfSumz$period)
dfSumz$period <- ifelse(dfSumz$YEAR>6,"4-10 Years After Sandy",dfSumz$period)
dfSumz$Group <- as.character(dfSumz$Group)

# Edit taxonomic groups names
dfSumz$Group2 <- dfSumz$Group
dfSumz$Group2 <- as.character(dfSumz$Group)
dfSumz$Group2 <- ifelse(dfSumz$Group2=="Euglenoid", "Euglenozoa",dfSumz$Group2)
dfSumz <- subset(dfSumz,!is.na(Group2))
ochro <- c("Tribophyte","Raphidophyte","Chrysophyte")
dfSumz$Group2 <- as.character(dfSumz$Group2)
dfSumz$Group2 <- ifelse(dfSumz$Group2%in%ochro, "Other Ochrophytes",dfSumz$Group2)
dfSumz$Group2 <- ifelse(dfSumz$Group2=="Euglenoid", "Euglenozoa",dfSumz$Group2)

# Find mean of each group for each month during each pre-post Sandy period
dfMeanz <- dfSumz%>%group_by(month,period,Group2)%>%summarize(m=mean(s))%>%as.data.frame()
dfMeanz$month <- as.numeric(as.character(dfMeanz$month))

# Change period to look at plots for other pre-post Sandy periods
p<- dfMeanz %>% filter(period=="Directly After Sandy")%>%ggplot(aes(x=month, y=m,fill=Group2))+geom_area(position="fill",color="#525252")+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic()+xlab("Month")+ylab("Relative Abundance")+scale_x_continuous(limits=c(4,11),breaks=c(4,5,6,7,8,9,10,11),labels=c("April","May","June","July","August","September","October","November"))+ggtitle("Directly After Sandy\n2013\nSite 6")+theme(axis.text.x = element_text(angle = 45,hjust=1))


# ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,common.legend = TRUE)
# ggsave("../../Fig3_Prelim.pdf",width=10,height=7)
```
## Next
