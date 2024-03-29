# RaritanBay_HAB_Networks
## By: Megan B. Rothenberger, Samantha J. Gleich, and Evan Flint
### Last Updated: June 20, 2023
Network analysis of the Raritan Bay Estuary phytoplankton community composition time-series dataset. See Rothenberger et al. (2014; https://doi.org/10.1007/s12237-013-9714-0) and  Rothenberger and Calomeni (2016; https://doi.org/10.1016/j.jembe.2016.03.015) for more information about the time-series dataset and previous work that has been conducted in this system. 

See Rothenberger et al. (2023; https://doi.org/10.1016/j.hal.2023.102411) for the manuscript associated with this analysis. 

## Required packages
To run the R scripts in this repository you will need to have the following packages installed: 
- `mgcv`
- `tidyverse`
- `compositions`
- `stats`
- `psych`
- `pulsar`
- `batchtools`
- `huge`
- `igraph`
- `reshape2`
- `randomcoloR`
- `ggpubr`
- `ggplot2`
- `NetGAM`
- `missForest`
- `patchwork`

## Make taxa barplots for each site
\
![](static/Figure2.png)\
**The R script "Figure2_Taxabarplot.R" will demonstrate how to wrangle and plot data to produce a plot like this one.** 


## Make taxa barplots pre-post Hurricane Sandy
\
![](static/Figure3.png)\
**The R script "Figure3_Taxabarplot.R" will demonstrate how to wrangle and plot data to produce a plot like this one.** 


## Make a network at each site
\
![](static/Figure4.png)\
**The R script "Figure4_Networks.R" will demonstrate how to produce networks like the those depicted here.**

## Look at the Site 6 HAB Associations in more detail (biotic vs. abiotic edges, positive vs. negative edges, and taxonomic breakdown of edges) 
\
![](static/Figure5.png)\
**The R script "Figure5_AssocBreakdown.R" will demonstrate how to wrangle and plot data to produce a plots like these.**

## Look at the variables that are associated with specific HAB species of interest.
\
![](static/Figure6.png)\
**The R script "Figure6_Subnetworks.R" will demonstrate how make subnetworks like these.**

