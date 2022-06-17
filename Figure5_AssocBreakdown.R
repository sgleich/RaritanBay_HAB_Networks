### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 16, 2022 ###
### Figure 5 - Site 6 detailed analysis ###

# Libraries
library(igraph)
library(patchwork)
library(tidyverse)
library(ggplot2)

# Continue from where you left off with Figure4_Networks.R
edgezHab <- subset(edge6,Var1%in%habz|Var2%in%habz)
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

# Biotic vs. Abiotic edges in Site 6 network
bioabio <- finMelt%>%ggplot(aes(x=Species,y=value,group=Species,fill=variable))+geom_bar(stat="identity",color="black")+scale_fill_manual(name="Association",values=c("gray","gray20"),label=c("Biotic","Abiotic"))+theme_classic(base_size=12)+scale_x_discrete(breaks=c("ALEX","ASAN","DINOP","HAKASH","HROT","HTRI","KARLO","OSCIL","PLIKE","PMIC","PMIN","PSEUDO","SCRIP","SNOW"),labels=c(expression(italic("Alexandrium spp.")),expression(italic("Akashiwo sanguinea")),expression(italic("Dinophysis acuminata")),expression(italic("Heterosigma akashiwo")),expression(italic("Heterocapsa rotundata")),expression(italic("Heterocapsa triquetra")),expression(italic("Karlodinium spp.")),expression(italic("Oscillatoria spp.")),expression(italic("Pfiestieria-like")),expression(italic("Prorocentrum micans")),expression(italic("Prorocentrum minimum")),expression(italic("Pseudo-nitzschia spp.")),expression(italic("Scrippsiella trochoidea")),expression(italic("Snowella spp."))))+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylab("Number of Associations")+xlab("HAB Species")+ylim(0,15)
# ggsave("Biotic_Abiotic_Site6.pdf",width=8,height=4)

fin2$Species <-rownames(fin2) 
finMelt2 <- melt(fin2,id.vars=c("Species"))

# Positive vs. negative edges in Site 6 network
posneg <- finMelt2%>%ggplot(aes(x=Species,y=value,group=Species,fill=variable))+geom_bar(stat="identity",color="black")+scale_fill_manual(name="Association",values=c("dodgerblue","indianred"))+theme_classic(base_size=12)+scale_x_discrete(breaks=c("ALEX","ASAN","DINOP","HAKASH","HROT","HTRI","KARLO","OSCIL","PLIKE","PMIC","PMIN","PSEUDO","SCRIP","SNOW"),labels=c(expression(italic("Alexandrium spp.")),expression(italic("Akashiwo sanguinea")),expression(italic("Dinophysis acuminata")),expression(italic("Heterosigma akashiwo")),expression(italic("Heterocapsa rotundata")),expression(italic("Heterocapsa triquetra")),expression(italic("Karlodinium spp.")),expression(italic("Oscillatoria spp.")),expression(italic("Pfiestieria-like")),expression(italic("Prorocentrum micans")),expression(italic("Prorocentrum minimum")),expression(italic("Pseudo-nitzschia spp.")),expression(italic("Scrippsiella trochoidea")),expression(italic("Snowella spp."))))+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylab("Number of Associations")+xlab("HAB Species")+ylim(0,15)

# Categorize HAB associations by taxonomic breakdown
dfCat <- read.csv("RB_Code.csv",header=TRUE)
colnames(dfCat) <- c("variable","Group")
ochro <- c("Tribophyte","Raphidophyte","Chrysophyte")
dfCat$Group <- ifelse(dfCat$Group%in%ochro, "Other Ochrophytes",dfCat$Group)
dfCat$Group <- ifelse(dfCat$Group=="Euglenoid", "Euglenozoa",dfCat$Group)
dfCat2 <- data.frame(variable=c("TEMP","pH","SAL","DO","Secchi","NO3","NH4","SRP","Si"),Group=c(rep("Environmental variable")))
dfCat <- rbind(dfCat,dfCat2)

dfFin <- NULL
for (num in c(1,3:14)){
  dfSub <- subset(edgezHab, Var1==habz[num]|Var2==habz[num])
  for(row in 1:nrow(dfSub)){
    if(dfSub$Var2[row]!=habz[num]){
      dfSub$Var1[row] <- dfSub$Var2[row]
      dfSub$Var2[row] <- habz[num]}
  }
  colnames(dfSub) <- c("variable","Var2","edge")
  dfSub2 <- left_join(dfSub,dfCat)
  dfSum <- dfSub2%>%group_by(Group)%>%tally()%>%as.data.frame()
  dfSum$HAB <- habz[num]
  dfFin <- rbind(dfFin,dfSum)
}

dfFin[nrow(dfFin)+1,] <- c("Diatom",0,"SNOW")
dfFin$n <- as.numeric(dfFin$n)

# Add environmental variable color into color palette
# colrs <- randomcoloR::distinctColorPalette(length(unique(dfFin$Group)))
colrs2 <- colrs
colrs2[6] <- "gray20"
colrs2[7:12] <- colrs[6:11]

# HAB edges associated with different taxonomic groups 
groups <- dfFin%>%ggplot(aes(x=HAB,y=n,group=HAB,fill=Group))+geom_bar(stat="identity",color="black")+theme_classic(base_size=12)+scale_x_discrete(breaks=c("ALEX","ASAN","DINOP","HAKASH","HROT","HTRI","KARLO","OSCIL","PLIKE","PMIC","PMIN","PSEUDO","SCRIP","SNOW"),labels=c(expression(italic("Alexandrium spp.")),expression(italic("Akashiwo sanguinea")),expression(italic("Dinophysis acuminata")),expression(italic("Heterosigma akashiwo")),expression(italic("Heterocapsa rotundata")),expression(italic("Heterocapsa triquetra")),expression(italic("Karlodinium spp.")),expression(italic("Oscillatoria spp.")),expression(italic("Pfiestieria-like")),expression(italic("Prorocentrum micans")),expression(italic("Prorocentrum minimum")),expression(italic("Pseudo-nitzschia spp.")),expression(italic("Scrippsiella trochoidea")),expression(italic("Snowella spp."))))+theme(axis.text.x = element_text(angle = 45,hjust=1))+ylab("Number of Associations")+xlab("HAB Species")+ylim(0,15)+scale_fill_manual(name="Taxonomic Group",values=c(colrs2))

# Combine 3 plots
bioabio+posneg+groups
# ggsave("Figure5.pdf",width=16,height=5)

