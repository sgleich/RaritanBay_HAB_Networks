### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 15, 2022 ###
### Figure 3 - Taxonomic breakdown by Sandy ###

# Libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(randomcoloR)

# Load
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

# Melt
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
colrs <- randomcoloR::distinctColorPalette(length(unique(dfMeanz1$Group)))

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

ggarrange(y1,y2,y3,y4,ncol=2,nrow=2,common.legend = TRUE)
# ggsave("Figure3.pdf",width=8,height=4)