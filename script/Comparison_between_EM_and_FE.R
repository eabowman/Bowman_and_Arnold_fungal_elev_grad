## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the relationship of EM and FE abundance and diversity to partially
## assess if there are any correlations between the two fungal communties.

#=========================================================================================
# Relationship of EM and FE abundance and diversity
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------

#--For abundance analyses
em.stat <- read.csv(paste0(dat.dir, 'SCM_EM_per_tree_stats.csv'), as.is = T, header = T)
abund.data <- read.csv(paste0(dat.dir, 'SCM_EM_root_based_site_x_species_matrix.csv'),
                       as.is = T, header = T)
fe.stat <- read.csv(paste0(dat.dir, 'SCM_FE_per_tree_stats.csv'), as.is = T, header = T)

#--For diversity analyses
em.div <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_landscape.csv'), as.is = T,
                   header = T)
fe.div <- read.csv(paste0(dat.dir, 'SCM_FE_div_pooled_by_landscape.csv'), as.is = T,
                   header = T)

#--results table
results.co <- data.frame(test = c('abundance', 'diversity'))
#-----------------------------------------------------------------------------------------
# FE abundance as a function of EM abundance
#-----------------------------------------------------------------------------------------
# << Clean up EM data files >> -----------------------------------------------------------
#--Uncomment to remove those trees with 0 abundance
#em.stat <- em.stat[which(!em.stat$abundance == 0.00), ]

#--Remove non-EM otus
root.out <- c("otu58", "otu125", "otu153", "otu44", "otu51", "otu68", "otu122", "otu135")
em.abun <- abund.data[which(!colnames(abund.data) %in% root.out)]

# << Clean up FE data files >> -----------------------------------------------------------
#--isolate fe abundance data
fe.abun <- fe.stat$abundance

# << Make data frame with EM and FE abundance >> -----------------------------------------
#--add LB012 to tree_no and make dataframe
tree_no <- fe.stat$tree_number
tree_no <- c(tree_no, 'LB012')
abund.co <- data.frame (tree_number = tree_no,
                        em.abund = NA,
                        fe.abund = NA)
#--add em abundance data
for (i in em.stat$tree_number) {
  abund.co[abund.co$tree_number==i, 'em.abund'] <-
    em.stat[em.stat$tree_number==i, 'abundance']
}
#--add em abundance data
for (i in fe.stat$tree_number) {
  abund.co[abund.co$tree_number==i, 'fe.abund'] <-
    fe.stat[fe.stat$tree_number==i, 'abundance']
}
#--repair em abundance for LB022 and LB036
abund.co[abund.co$tree_number %in% c('LB022','LB036','LB061','LB028'),
         'em.abund'] <- 0.0000001
#--repair fe abundance data
abund.co[abund.co$tree_number == 'LB012', 'fe.abund'] <- 0.000001
abund.co[abund.co$fe.abund == 0, 'fe.abund'] <- 0.00000001
#--LB018, LB049, LB050, LB006, LB012 all NA, because no needles collected
#abund.co[abund.co$tree_number %in% c('LB006','LB018','LB049','LB050','LB012'),
#         'fe.abund'] <- NA

# << Linear regression >> ----------------------------------------------------------------
lm.ab.co <- lm(logit(abund.co$fe.abund) ~ log (abund.co$em.abund))
abun.co <- summary(lm.ab.co)
abun.co

#--add results to results table
results.co[results.co$test == 'abundance', 'df.1'] <- abun.co$fstatistic[2]
results.co[results.co$test == 'abundance', 'df.2'] <- abun.co$fstatistic[3]
results.co[results.co$test == 'abundance', 'R.sq'] <- abun.co$r.squared[1]
results.co[results.co$test == 'abundance', 'f.stat'] <- abun.co$fstatistic[1]
results.co[results.co$test == 'abundance', 'p'] <- abun.co$coefficients[2,4]

# << Plot >> -----------------------------------------------------------------------------
jpeg(filename = paste0(fig.dir, 'AppendixS13.jpeg'), width = 1400, height = 700)
#--Combine into one figure
par(mfrow = c(1,2), "mar"=c(6, 5, 5, 6))
#--plot abundance as a function of elevation
plot(abund.co$em.abund, abund.co$fe.abund, pch = 19,
      ylab = "FE abundance", xlab = "EM abundance", cex.lab = 1.5, cex.axis = 1.5)
#abline (m <- lm(abund.co$fe.abund ~ abund.co$em.abund))

#-----------------------------------------------------------------------------------------
# FE diversity as a function of EM diversity
#-----------------------------------------------------------------------------------------
# << Clean up data >> --------------------------------------------------------------------
#--Make data frame with both EM and FE diversity
div.co <- data.frame(elevation = em.div$Elevation,
                     topography = em.div$topography,
                     em.fa = em.div$fishers.alpha,
                     fe.fa = NA)
#--add FE data
for (e in div.co$elevation) {
  div.co[div.co$elevation==e, 'fe.fa'] <-
    fe.div[fe.div$Elevation==e, 'fishers.alpha']
}

#--Uncomment to remove outliers
#div.co <- div.co[-c(12,17),]

# << Linear regression >> ----------------------------------------------------------------
lm.div.co <- lm(log(div.co$fe.fa) ~ log(div.co$em.fa))
di.co <- summary(lm.div.co)
di.co

#--add results to results table
results.co[results.co$test == 'diversity', 'df.1'] <- di.co$fstatistic[2]
results.co[results.co$test == 'diversity', 'df.2'] <- di.co$fstatistic[3]
results.co[results.co$test == 'diversity', 'R.sq'] <- di.co$r.squared[1]
results.co[results.co$test == 'diversity', 'f.stat'] <- di.co$fstatistic[1]
results.co[results.co$test == 'diversity', 'p'] <- di.co$coefficients[2,4]

#--write out results table
write.csv(results.co, paste0(res.dir, 'Comparison_of_EM_and_FE.csv'), row.names = F)

# << Plot >> -----------------------------------------------------------------------------
#--plot abundance as a function of elevation
plot (div.co$em.fa, div.co$fe.fa, pch = 19,
      ylab = "FE Fisher's alpha", xlab = "EM Fisher's alpha",
      cex.lab = 1.5, cex.axis = 1.5)
#abline (m <- lm(abund.co$fe.abund ~ abund.co$em.abund))
#--close graph
dev.off()
