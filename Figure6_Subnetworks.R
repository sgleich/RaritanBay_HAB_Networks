### Raritan Bay ENA ###
### By: Samantha Gleich ###
### Last Updated: June 17, 2022 ###
### Figure 6 - HAB subnetworks at Site 6 ###

# Libraries
library(igraph)
library(tidyverse)

# Start with this site 6 network and then separate it into HAB subnetworks
plot(netwk6,vertex.label=NA,vertex.size=7,vertex.color=V(netwk6)$col,layout=layout_with_fr(netwk6),main="Site 6")

# Heterosigma akashiwo
hakas <- subset(edge6,Var1=="HAKASH"|Var2=="HAKASH")
d <- data.frame(hakas$Var1,hakas$Var2)
hakashG <- graph_from_edgelist(as.matrix(d),directed=FALSE)
E(hakashG)$weight <- hakas$edge
E(hakashG)$colz <- ifelse(E(hakashG)$weight>0,"deepskyblue4","indianred")
V(hakashG)$col <- ifelse(V(hakashG)$name=="HAKASH","aquamarine","grey")
V(hakashG)$col <- ifelse(V(hakashG)$name=="SCRIP","darkgoldenrod3",V(hakashG)$col)
V(hakashG)$col <- ifelse(V(hakashG)$name=="PMIC","orchid",V(hakashG)$col)
V(hakashG)$name <- ifelse(V(hakashG)$name=="HAKASH","Heterosigma akashiwo",V(hakashG)$name)
V(hakashG)$name <- ifelse(V(hakashG)$name=="CHLAMY","Chlamydomonas spp.",V(hakashG)$name)
V(hakashG)$name <- ifelse(V(hakashG)$name=="CYCLO","Cyclotella spp.",V(hakashG)$name)
V(hakashG)$name <- ifelse(V(hakashG)$name=="PHAC","Phacus spp.",V(hakashG)$name)
V(hakashG)$name <- ifelse(V(hakashG)$name=="PMIC","Prorocentrum micans",V(hakashG)$name)
V(hakashG)$name <- ifelse(V(hakashG)$name=="SCRIP","Scrippsiella trochoidea",V(hakashG)$name)

pdf("HAKASH.pdf",width=4,height=4)
plot(hakashG,vertex.label=V(hakashG)$name,vertex.label.cex=0.6,vertex.label.dist=3,vertex.label.degree=c(120,100,84,90,84,89),vertex.label.color="black",vertex.size=9,vertex.label.font=c(4,4,4,4,4,4),vertex.color=V(hakashG)$col,main=expression(italic("Heterosigma akashiwo")),layout=layout_in_circle(hakashG),edge.color=E(hakashG)$colz,edge.width=(abs(E(hakashG)$weight)*6))
dev.off()

# Oscillatoria spp.
oscil <- subset(edge6,Var1=="OSCIL"|Var2=="OSCIL")
d <- data.frame(oscil$Var1,oscil$Var2)
oscilG <- graph_from_edgelist(as.matrix(d),directed=FALSE)
E(oscilG)$weight <- oscil$edge
E(oscilG)$colz <- ifelse(E(oscilG)$weight>0,"deepskyblue4","indianred")
V(oscilG)$col <- ifelse(V(oscilG)$name=="OSCIL","cyan","grey")
V(oscilG)$name <- ifelse(V(oscilG)$name=="OSCIL","Oscillatoria spp.",V(oscilG)$name)
V(oscilG)$name <- ifelse(V(oscilG)$name=="AULA","Aulacoseira spp.",V(oscilG)$name)
V(oscilG)$name <- ifelse(V(oscilG)$name=="CHMIN","Chaetoceros minimus",V(oscilG)$name)
V(oscilG)$name <- ifelse(V(oscilG)$name=="GYROS","Gyrosigma spp.",V(oscilG)$name)
V(oscilG)$name <- ifelse(V(oscilG)$name=="ODON","Odontella spp.",V(oscilG)$name)

pdf("OSCIL.pdf",width=4,height=4)
plot(oscilG,vertex.label=V(oscilG)$name,vertex.label.cex=0.6,vertex.label.dist=3,vertex.label.degree=c(120,100,84,90,84),vertex.label.color="black",vertex.size=9,vertex.label.font=c(4,4,4,4,4),vertex.color=V(oscilG)$col,main=expression(italic("Oscillatoria spp.")),layout=layout_in_circle(oscilG),edge.color=E(oscilG)$colz,edge.width=(abs(E(oscilG)$weight)*6))
dev.off()

# Akashiwo sanguinea
asan <- subset(edge6,Var1=="ASAN"|Var2=="ASAN")
d <- data.frame(asan$Var1,asan$Var2)
asanG <- graph_from_edgelist(as.matrix(d),directed=FALSE)
E(asanG)$weight <- asan$edge
E(asanG)$colz <- ifelse(E(asanG)$weight>0,"deepskyblue4","indianred")

V(asanG)$col <- ifelse(V(asanG)$name=="ASAN","firebrick1","grey")
V(asanG)$col <- ifelse(V(asanG)$name=="ALEX","hotpink4",V(asanG)$col)
V(asanG)$col <- ifelse(V(asanG)$name=="HROT","khaki1",V(asanG)$col)
V(asanG)$col <- ifelse(V(asanG)$name=="Secchi","gray20",V(asanG)$col)

V(asanG)$name <- ifelse(V(asanG)$name=="ASAN","Akashiwo sanguinea",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="BOSMI","Bosmina spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="EUPLO","Euplotes spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="HEMICYC","Hemicyclops spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="PLEOP","Pleopsis spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="STROM"," Strombidium spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="XANZOEA","Unknown Arthropod",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="ALEX","Alexandrium spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="AMPH","Amphidinium spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="GYROD","Gyrodinium spp.",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="HROT","Heterocapsa rotundata",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="MOON","Unknown Dinoflagellate",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="OXMAR","Oxhyrris marina",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="CYLCLO","Cylindrotheca closterium",V(asanG)$name)
V(asanG)$name <- ifelse(V(asanG)$name=="Secchi","Secchi Depth",V(asanG)$name)


pdf("ASAN.pdf",width=4,height=4)
plot(asanG,vertex.label=V(asanG)$name,vertex.label.cex=0.6,vertex.label.dist=3,vertex.label.degree=c(120,100,100,100,92,79,79,79,86,46,40,90,88.5,88.4,88.6),vertex.label.color="black",vertex.size=9,vertex.label.font=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4),vertex.color=V(asanG)$col,main=expression(italic("Akashiwo sanguinea")),layout=layout_in_circle(asanG),edge.color=E(asanG)$colz,edge.width=(abs(E(asanG)$weight)*6))
dev.off()

# Dinophysis acuminata 
dinop<- subset(edge6,Var1=="DINOP"|Var2=="DINOP")
d <- data.frame(dinop$Var1,dinop$Var2)
dinopG <- graph_from_edgelist(as.matrix(d),directed=FALSE)
E(dinopG)$weight <- dinop$edge
E(dinopG)$colz <- ifelse(E(dinopG)$weight>0,"deepskyblue4","indianred")

V(dinopG)$col <- ifelse(V(dinopG)$name=="DINOP","cornflowerblue","grey")
V(dinopG)$col <- ifelse(V(dinopG)$name=="SRP","gray20",V(dinopG)$col)

V(dinopG)$name <- ifelse(V(dinopG)$name=="DINOP","Dinophysis acuminata",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="ASTAD","Oikopleura spp.",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="BRACH","Brachionus calyciflorus",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="COPPART","Unknown Arthropod",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="LABO","Laboea",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="PODON","Pleopsis polyphemoides",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="STROB","Strobilidium spp.",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="TEMOR","Temora longicornis",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="AULA","Aulacoseira spp.",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="PHAC","Phacus spp.",V(dinopG)$name)
V(dinopG)$name <- ifelse(V(dinopG)$name=="SRP","SRP",V(dinopG)$name)


pdf("DINOP.pdf",width=4,height=4)
plot(dinopG,vertex.label=V(dinopG)$name,vertex.label.cex=0.6,vertex.label.dist=3,vertex.label.degree=c(120,75,100,80,91,78,78,90,78,46,38.9),vertex.label.color="black",vertex.size=9,vertex.label.font=c(4,4,4,4,4,4,4,4,4,4,4),vertex.color=V(dinopG)$col,main=expression(italic("Dinophysis acuminata")),layout=layout_in_circle(dinopG),edge.color=E(dinopG)$colz,edge.width=(abs(E(dinopG)$weight)*6))
dev.off()


