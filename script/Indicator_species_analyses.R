## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate whether any fungal OTUs were significantly associated with 
## high or low elevation sites and the biotic and abiotic characteristics associated 
## with each elevation group (i.e. differences in soil chemistry, climate, and plant
## community).

#=========================================================================================
# EM indicator species based on OTU count
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------
otu.data <- read.csv(paste0(dat.dir, 'SCM_EM_otu_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#-----------------------------------------------------------------------------------------
# Indicator species analysis
#-----------------------------------------------------------------------------------------
#<< Overall EM, with singletons included >>-----------------------------------------------
#--isolate OTUs with great than 4 ocurrences
em.indsp <- otu.data[6:length(otu.data)]
em.indsp <- em.indsp[, colSums(em.indsp) > 4]
#--pull out elevation group data for grouping
Pl_comm = otu.data$plant.comm
#--Association between species patterns and high/low elevation sites
indsp.em <- multipatt(em.indsp, Pl_comm, func = "IndVal") # uses default func = "IndVal"
#--print summary
summary(indsp.em)
#--Components A and B. 
summary(indsp.em, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.em, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.em$sign 

#=========================================================================================
# EM indicator species based on abundance data
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------
abund.data <- read.csv(paste0(dat.dir, 'SCM_EM_root_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#-----------------------------------------------------------------------------------------
# Indicator analysis
#-----------------------------------------------------------------------------------------
#--isolate OTUs with great than 4 ocurrences
abund.indsp <- abund.data[6:length(abund.data)]
abund.indsp <- abund.indsp[, colSums(abund.indsp) > 4]
#--pull out elevation group data for grouping
El_group = abund.data$plant.comm
#--Association between species patterns and high/low elevation sites
indsp.abund <- multipatt(abund.indsp, El_group) # uses default func = "IndVal"
#--print summary
summary(indsp.abund)
#--Components A and B. 
summary(indsp.abund, indvalcomp = TRUE) 
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.abund, alpha = 1)
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.abund$sign 

#=========================================================================================
# FE indicator species based on OTU count
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------
fe.data <- read.csv(paste0(dat.dir, 'SCM_FE_site_x_species_matrix.csv'),
                       as.is = T, header = T)
#-----------------------------------------------------------------------------------------
# Indicator analysis
#-----------------------------------------------------------------------------------------
#<< Overall FE, with singletons included >>-----------------------------------------------
#--isolate OTUs with greater with 4 ocurrences
fe.indsp <- fe.data[2:41]
fe.indsp <- fe.indsp[, colSums(fe.indsp) > 4]
#--pull out elevation group data for grouping
El_group = fe.data$plant.comm
#--Association between species patterns and high/low elevation sites
indsp.fe <- multipatt(fe.indsp, El_group) # uses default func = "IndVal"
#--print summary
summary(indsp.fe)
#--Components A and B. 
summary(indsp.fe, indvalcomp = TRUE)
# A = the probability that the surveyed site belongs to the target site group given the 
# fact that the species has been found.
# B = the probability of finding the species in sites belonging to the site group.

#--shows all species with lower indval, not just sig.
summary(indsp.fe, alpha = 1) 
#--shows all spp. even those with highest indval that do not show up in other summaries
indsp.em$sign
