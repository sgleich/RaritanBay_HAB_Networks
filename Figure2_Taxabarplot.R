### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 15, 2022 ###
### Figure 2 - Taxonomic breakdown by season ###

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

# Melt
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

ggarrange(ste1,ste2,ste3,ste4,ste5,ste6,ncol=3,nrow=2,common.legend = TRUE)
# ggsave("Figure2.pdf",width=12,height=8)