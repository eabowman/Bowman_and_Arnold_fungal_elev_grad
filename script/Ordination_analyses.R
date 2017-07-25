## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate fungal community patterns.

#=========================================================================================
# Permanova and Adonis
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data 
#-----------------------------------------------------------------------------------------

#--EM site x species matrix based on OTU count
otu.data <- read.csv(paste0(dat.dir,'SCM_EM_otu_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#--EM site x species matrix based on OTU count
abund.data <- read.csv(paste0(dat.dir,'SCM_EM_root_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#--FE site x species matrix based on OTU count
fe.data <- read.csv(paste0(dat.dir,'SCM_FE_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#--Environmental data for PERMANOVA analyses
env.data <- read.csv(paste0(dat.dir,'SCM_environmental_data.csv'))

#--make soil data numeric
env <- c("elevation","ph.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm", "k.ppm","zn.ppm",
         "fe.ppm","mn.ppm","cu.ppm","ni.ppm", "po4.p.ppm","so4.s.ppm","b.ppm","esp",
         "cec.meq.100g","northing","easting")
for (i in env) {
  env.data[[i]] <- as.numeric(env.data[[i]])
}

#--Permanova data-----
#--data for jaccard based dissimilarity matrix
jac.perm <- read.csv(paste0(dat.dir,
                           'SCM_EM_otu_site_x_species_clustered_bytopoandelev.csv'),
         as.is = T, header = T)

#--data for Morisita-horn based dissimilarity index
mor.perm <- read.csv(paste0(dat.dir,
                            'SCM_EM_root_site_x_species_clustered_bytopoandelev.csv'),
                     as.is = T, header = T)

#-----------------------------------------------------------------------------------------
# Outlier trees and non-EM OTUs for possible removal later
#-----------------------------------------------------------------------------------------
#--Outliers to remove for EM fungi
em.out <- c("LB025", "LB005", "LB034", "LB038", "LB040", "LB050")
#--Outliers to remove for FE fungi
fe.out <- c("LB043", "LB044")
#--Remove non-EM otus
root.out <- c("otu58", "otu125", "otu153", "otu44", "otu51", "otu68", "otu122", "otu135")

#-----------------------------------------------------------------------------------------
# Create result tables
#-----------------------------------------------------------------------------------------
#--Table of Anosim results
anosim.res <- c("root.j", "root.m", "em.j","em.m", "fe.j", "fe.m")
anosim.res <- data.frame(anosim.res)
#--Table of PERMANOVA results
permanova.res <- data.frame(c("em.j","em.m"))
colnames(permanova.res) <- "test"

#=========================================================================================
# Root fungi-Jaccard
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean up data
#-----------------------------------------------------------------------------------------
#<< Jaccard index, OTU data >>------------------------------------------------------------
#--all root fungi in site x species matrix based on otu count
nmds.rooto <- otu.data
#--Comment to include outliers 
nmds.rooto <- nmds.rooto[!(nmds.rooto$tree_number %in% em.out),]

#-----------------------------------------------------------------------------------------
# NMDS plot and Adonis
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.rooto[5:length(nmds.rooto)]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(!(colSums(comm.matrix) == 0))]

# << to remove singletons uncomment whole section, if not skip this section >> -----------
#--remove singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--if removed singletons then need to check row sums, first identify which rows are being
#--removed so same rows can be removed from supporing files for later analyses
#rows.removed <- which(rowSums(comm.matrix) == 0)
#comm.matrix <- comm.matrix[-rows.removed,]
#--remove same rows from nmds.rooto file for later analyses
#nmds.rooto <- nmds.rooto[-rows.removed,]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       try = 100, trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="root.j"), "stress.nmds"] <-
  jaccard.otu$stress

#--Plot NMDS of root community based on Jaccard index and OTU abundance
#plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5,
#     cex.axis = 1.5, yaxt = "n")
#axis(2, at = seq(-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
#color.vec <- data.frame(color = rep(NA, length(rownames(comm.matrix))),
#                                    p.group = nmds.rooto$group)
#color.vec <- sapply(color.vec$p.group, function (x) 
#  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
#ordipointlabel(morisita.abund, display = "sites")
#points(jaccard.otu, display = "sites", cex = 2, pch = 20,
#       col = color.vec,
#       bg = color.vec)
#legend("topright", legend = levels(otu.env.data$group), bty = "n",
#       col = color.vec, pch = 21, pt.bg = color.vec, cex = 2)
#ordihull(jaccard.otu, groups = nmds.rooto$group)

#<< BetaDisper: a multivariate analogue of Levene's test for homogeneity of variance >>---
betadisper <- betadisper(comm.dist.jaccard, group = nmds.rooto$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
root.j.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="root.j"), "F.betadisper"] <-
  root.j.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="root.j"), "df.betadisper.1"] <-
  root.j.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="root.j"), "df.betadisper.2"] <-
  root.j.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="root.j"), "p.betadisper"] <-
  root.j.homo$`Pr(>F)`[1]

#<< ANOSIM >> ----------------------------------------------------------------------------
# << ANOSIM: parameters for anosim not met; non-parametric analysis used instead >> ------
#--PERMANOVA: non-parametric
root.jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ group,
                            data = nmds.rooto)
#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="root.j"), "permanova.f"] <-
  root.jaccard.adonis$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res =="root.j"), "permanova.r2"] <-
  root.jaccard.adonis$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res =="root.j"), "permanova.p"] <-
  root.jaccard.adonis$aov.tab$`Pr(>F)`[1]

#=========================================================================================
# Root fungi-Morisita
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean up data
#-----------------------------------------------------------------------------------------
nmds.roota <- abund.data
#--Comment to include outliers
#--remove rows from otu.table and comm.matrix; outliers are LB005, LB034, LB038, LB040
nmds.roota <- nmds.roota[!(nmds.roota$tree_number %in% em.out),]
#-----------------------------------------------------------------------------------------
# NMDS plot and Adonis
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.roota[5:length(nmds.roota)]
#--Uncomment remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]

# << to remove singletons uncomment whole section, if not skip this section >> -----------
#--remove singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]

#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)
#--NMDS analysis
morisita.abund <- metaMDS(comm.dist.horn, dist = "bray", permutations = 999, try = 50, 
                          trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="root.m"), "stress.nmds"] <-
  morisita.abund$stress

#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'Figure5.jpeg'), width = 1200, height = 450,
     quality = 100)
par(mfrow = c(1,3), "mar"=c(6, 5, 5, 3))

#--Plot NMDS of root community based on morisita index and OTU abundance
plot(morisita.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis(2, at = seq(-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep (NA, length(rownames(comm.matrix))),
                                     p.group = nmds.roota$group)
color.vec <- sapply(color.vec$p.group, function(x) 
  if(x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
#ordipointlabel(morisita.abund, display = "sites")
points(morisita.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
#legend("topright", legend = levels(otu.env.data$group), bty = "n",
#       col = color.vec, pch = 21, pt.bg = color.vec, cex = 2)
ordihull(morisita.abund, groups = nmds.roota$group)

# <<BetaDisper: a multivariate analogue of Levene's test for homogeneity of variance >>---
betadisper <- betadisper(comm.dist.horn, group = nmds.roota$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
root.m.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="root.m"), "F.betadisper"] <-
  root.m.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="root.m"), "df.betadisper.1"] <-
  root.m.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="root.m"), "df.betadisper.2"] <-
  root.m.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="root.m"), "p.betadisper"] <-
  root.m.homo$`Pr(>F)`[1]

# << ANOSIM >>----------------------------------------------------------------------------
root.horn.anosim <- anosim(comm.matrix, grouping = nmds.roota$group,
                            distance = "morisita")
#--Add results to data frame
anosim.res [ which (anosim.res$anosim.res =="root.m"), "r"] <-
  root.horn.anosim$statistic
anosim.res [ which (anosim.res$anosim.res =="root.m"), "p"] <- root.horn.anosim$signif

#=========================================================================================
# Ectomycorrhizal fungi-Jaccard
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean up data
#-----------------------------------------------------------------------------------------
nmds.otu <- otu.data

#--Comment to include outliers
nmds.otu <- nmds.otu[!(nmds.otu$tree_number %in% em.out), ]
#--remove non-em fungi only
nmds.otu <- nmds.otu[,!colnames(nmds.otu) %in% root.out ]

#-----------------------------------------------------------------------------------------
# NMDS plot and Adonis
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.otu[5:length(nmds.otu)]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]

# << to remove singletons uncomment whole section, if not skip this section >> -----------
#--remove singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--if removed singletons then need to check row sums, first identify which rows are being
#--removed so same rows can be removed from supporing files for later analyses
#rows.removed <- which(rowSums(comm.matrix) == 0)
#comm.matrix <- comm.matrix[-rows.removed,]
#--remove same rows from nmds.rooto file for later analyses
#nmds.rooto <- nmds.rooto[-rows.removed,]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999, 
                       try = 100, trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="em.j"), "stress.nmds"] <-
  jaccard.otu$stress

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
#plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5, cex.axis = 1.5)
#color.vec <- data.frame(color = rep(NA, length(rownames(comm.matrix))),
#                        p.group = nmds.otu$group)
#color.vec <- sapply(color.vec$p.group, function(x) 
#  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
#points(jaccard.otu, display = "sites", cex = 2, pch = 20,
#       col = color.vec,
#       bg = color.vec)
#ordihull(jaccard.otu, groups = nmds.otu$group)

# << BetaDisper: a multivariate analogue of Levene's test for homogeneity of variance >>--
betadisper <- betadisper(comm.dist.jaccard, group = nmds.otu$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
em.j.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="em.j"), "F.betadisper"] <-
  em.j.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="em.j"), "df.betadisper.1"] <-
  em.j.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="em.j"), "df.betadisper.2"] <-
  em.j.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="em.j"), "p.betadisper"] <-
  em.j.homo$`Pr(>F)`[1]

# << ANOSIM: parameters for anosim not met; non-parametric analysis used instead >> ------
#--PERMANOVA: non-parametric
em.jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ group,
                            data = nmds.otu)
#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="em.j"), "permanova.f"] <-
  em.jaccard.adonis$aov.tab$F.Model[1]
anosim.res[which(anosim.res$anosim.res =="em.j"), "permanova.r2"] <-
  em.jaccard.adonis$aov.tab$R2[1]
anosim.res[which(anosim.res$anosim.res =="em.j"), "permanova.p"] <-
  em.jaccard.adonis$aov.tab$`Pr(>F)`[1]

#-----------------------------------------------------------------------------------------
# PERMANOVA
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- jac.perm[8:length(jac.perm)]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--adonis
em.jaccard.adonis <- adonis(formula = comm.dist.jaccard ~ Average.temperature.warm.quarter
               + Phosphate + Plant.community, data = jac.perm, permutations = 1000,
               strata = jac.perm$group)
#--add results to data.frame
#--Avg. temperature warm quarter f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.j"), "F.model.avg.temp.wm.quarter"] <- 
  em.jaccard.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "em.j"), "r2.avg.temp.wm.quarter"] <- 
  em.jaccard.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "em.j"), "p.avg.temp.wm.quarter"] <- 
  em.jaccard.adonis$aov.tab$`Pr(>F)`[1]
#--PO4.P f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.j"), "F.model.po4"] <- 
  em.jaccard.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "em.j"), "r2.po4"] <- 
  em.jaccard.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "em.j"), "p.po4"] <- 
  em.jaccard.adonis$aov.tab$`Pr(>F)`[2]
#--Plant community f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.j"), "F.model.pl.comm"] <- 
  em.jaccard.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "em.j"), "r2.pl.comm"] <- 
  em.jaccard.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "em.j"), "p.pl.comm"] <- 
  em.jaccard.adonis$aov.tab$`Pr(>F)`[3]

#--remove temporary data frames, vectors, etc...
rm(comm.matrix, jaccard.otu, comm.dist.jaccard, nmds.otu)

#=========================================================================================
# Ectomycorrhizal fungi-Morisita
#=========================================================================================
#--Remove non-em fungi only from dataframe
nmds.ab <- abund.data[!(abund.data$tree_number %in% root.out),]
#--Comment to include outliers and non-em fungi
#--remove rows from otu.table and comm.matrix; outliers are LB005, LB034, LB038, LB040
nmds.ab <- nmds.ab [!(nmds.ab$tree_number %in% em.out),]

#-----------------------------------------------------------------------------------------
# NMDS and ANOSIM
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.ab[5:length(nmds.ab)]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]

# << to remove singletons uncomment whole section, if not skip this section >> -----------
#--remove singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--if removed singletons then need to check row sums, first identify which rows are being
#--removed so same rows can be removed from supporing files for later analyses
#rows.removed <- which(rowSums(comm.matrix) == 0)
#comm.matrix <- comm.matrix[-rows.removed,]
#--remove same rows from nmds.rooto file for later analyses
#nmds.rooto <- nmds.rooto[-rows.removed,]

#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)

#--NMDS analysis
horn.abund <- metaMDS(comm.dist.horn, dist = "bray", permutations = 999, try = 50, 
                      trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="em.m"), "stress.nmds"] <-
  horn.abund$stress

#--Plot NMDS of EM community based on Morisita index and OTU abundance
plot(horn.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
#--color for groups
color.vec <- data.frame (color = rep (NA, length(rownames(comm.matrix))),
                         p.group = nmds.ab$group)
color.vec <- sapply (color.vec$p.group, function (x) 
  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
points(horn.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
ordihull(horn.abund, groups = nmds.ab$group)

#<< BetaDisper: a multivariate analogue of Levene's test for homogeneity of variance >>---
betadisper <- betadisper(comm.dist.horn, group = nmds.ab$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
em.m.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="em.m"), "F.betadisper"] <-
  em.m.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="em.m"), "df.betadijac.sper.1"] <-
  em.m.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="em.m"), "df.betadisper.2"] <-
  em.m.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="em.m"), "p.betadisper"] <-
  em.m.homo$`Pr(>F)`[1]

# << ANOSIM >>----------------------------------------------------------------------------
em.horn.anosim <- anosim(comm.matrix, grouping = nmds.ab$group, distance = "morisita")
#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="em.m"), "r"] <- em.horn.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="em.m"), "p"] <- em.horn.anosim$signif

#-----------------------------------------------------------------------------------------
# PERMANOVA
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- mor.perm[8:length(mor.perm)]

#--distance matrix using morisita index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)

#--adonis
em.horn.adonis <- adonis(formula = comm.dist.horn ~ Average.temperature.warm.quarter +
                         Plant.community + Phosphate, data = mor.perm,
                         permutations = 1000, strata = mor.perm$group)
#--add results to data.frame
#--Avg. temperature warm quarter f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.m"), "F.model.avg.temp.wm.quarter"] <- 
  em.horn.adonis$aov.tab$F.Model[1]
permanova.res[which(permanova.res$test == "em.m"), "r2.avg.temp.wm.quarter"] <- 
  em.horn.adonis$aov.tab$R2[1]
permanova.res[which(permanova.res$test == "em.m"), "p.avg.temp.wm.quarter"] <- 
  em.horn.adonis$aov.tab$`Pr(>F)`[1]
#--PO4.P f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.m"), "F.model.po4"] <- 
  em.horn.adonis$aov.tab$F.Model[2]
permanova.res[which(permanova.res$test == "em.m"), "r2.po4"] <- 
  em.horn.adonis$aov.tab$R2[2]
permanova.res[which(permanova.res$test == "em.m"), "p.po4"] <- 
  em.horn.adonis$aov.tab$`Pr(>F)`[2]
#--Plant community f.model, r2, p-value
permanova.res[which(permanova.res$test == "em.m"), "F.model.pl.comm"] <- 
  em.horn.adonis$aov.tab$F.Model[3]
permanova.res[which(permanova.res$test == "em.m"), "r2.pl.comm"] <- 
  em.horn.adonis$aov.tab$R2[3]
permanova.res[which(permanova.res$test == "em.m"), "p.pl.comm"] <- 
  em.horn.adonis$aov.tab$`Pr(>F)`[3]

#=========================================================================================
# Endophytic fungi-Jaccard
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean up data
#-----------------------------------------------------------------------------------------
nmds.fe <- fe.data
#--remove rows from otu.table and comm.matrix; outliers are LB043, LB044
#--Comment to include outliers
nmds.fe <- nmds.fe[!(nmds.fe$tree_number %in% fe.out),]

#-----------------------------------------------------------------------------------------
# NMDS and Anosim
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.fe[2:41]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]
#--comment to add singleton OTUs
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray", permutations = 999,
                       trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="fe.j"), "stress.nmds"] <-
  jaccard.otu$stress

#--Plot NMDS of FE community based on Jaccard index and OTU abundance
#plot(jaccard.otu, display = "sites", type = "n", cex.lab = 1.5, cex.axis = 1.5)
#--color for groups
#color.vec <- data.frame (color = rep (NA, length(rownames(comm.matrix))),
#                         p.group = nmds.fe$group)
#color.vec <- sapply (color.vec$p.group, function (x) 
#  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
#points(jaccard.otu, display = "sites", cex = 2, pch = 20,
#       col = color.vec,
#       bg = color.vec)
#ordihull(jaccard.otu, groups = nmds.fe$group)

# << BetaDisper: a multivariate analogue of Levene's test for homogeneity of variance >>--
betadisper <- betadisper(comm.dist.jaccard, group = nmds.fe$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
fe.j.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="fe.j"), "F.betadisper"] <-
  fe.j.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="fe.j"), "df.betadisper.1"] <-
  fe.j.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="fe.j"), "df.betadisper.2"] <-
  fe.j.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="fe.j"), "p.betadisper"] <-
  fe.j.homo$`Pr(>F)`[1]

# << ANOSIM >>----------------------------------------------------------------------------
fe.jaccard.anosim <- anosim(comm.matrix, grouping = nmds.fe$group,
                             distance = "jaccard")
#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="fe.j"), "r"] <-
  fe.jaccard.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="fe.j"), "p"] <- fe.jaccard.anosim$signif

#=========================================================================================
# Endophytic fungi-Morisita
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean up data
#-----------------------------------------------------------------------------------------
nmds.fe <- fe.data
#--remove rows from otu.table and comm.matrix; outliers are LB043, LB044
#--Commend to include outliers
nmds.fe <- nmds.fe[!(nmds.fe$tree_number %in% fe.out),]

#--make soil data numeric
env <- c("elevation","ph.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm", "k.ppm","zn.ppm",
         "fe.ppm","mn.ppm","cu.ppm","ni.ppm", "po4.p.ppm","so4.s.ppm","b.ppm","esp",
         "cec.meq.100g","northing","easting")
for (i in env) {
  nmds.fe[[i]] <- as.numeric (nmds.fe [[i]])
}

#-----------------------------------------------------------------------------------------
# NMDS and Anosim
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.fe[2:41]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]
#--comment to add singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--distance matrix using morisita index
comm.dist.horn <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
horn.abund <- metaMDS(comm.dist.horn, dist = "bray", permutations = 999, try = 50, 
                      trymax = 1000)
#--add stress of NMDS to results table
anosim.res[which(anosim.res$anosim.res =="fe.m"), "stress.nmds"] <-
  horn.abund$stress

#--Plot NMDS of FE community based on Morisita index and OTU abundance
plot(horn.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5)
#--add color groups
color.vec <- data.frame (color = rep (NA, length(rownames(comm.matrix))),
                         p.group = nmds.fe$group)
color.vec <- sapply (color.vec$p.group, function (x) 
  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
#color.vec <- as.data.frame (color.vec)
#ordipointlabel(morisita.abund, display = "sites")
points(horn.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
#legend("topright", legend = levels(otu.env.data$group), bty = "n",
#       col = color.vec, pch = 21, pt.bg = color.vec, cex = 2)
ordihull(horn.abund, groups = nmds.fe$group)

#--output figure
dev.off()

#<< BetaDisper:a multivariate analogue of Levene's test for homogeneity of variance >> ---
betadisper <- betadisper(comm.dist.horn, group = nmds.fe$group)
#--ANOVA to assess if the variances are different, the distances of group members to the
#--group centroid are subject to ANOVA
fe.m.homo <- anova(betadisper)
#--add results to table
anosim.res[which(anosim.res$anosim.res =="fe.m"), "F.betadisper"] <-
  fe.m.homo$`F value`[1]
anosim.res[which(anosim.res$anosim.res =="fe.m"), "df.betadisper.1"] <-
  fe.m.homo$Df[1]
anosim.res[which(anosim.res$anosim.res =="fe.m"), "df.betadisper.2"] <-
  fe.m.homo$Df[2]
anosim.res[which(anosim.res$anosim.res =="fe.m"), "p.betadisper"] <-
  fe.m.homo$`Pr(>F)`[1]

# << ANOSIM >> ---------------------------------------------------------------------------
fe.horn.anosim <- anosim(comm.matrix, grouping = nmds.fe$group, distance = "horn")

#--Add results to data frame
anosim.res[which(anosim.res$anosim.res =="fe.m"), "r"] <-
  fe.horn.anosim$statistic
anosim.res[which(anosim.res$anosim.res =="fe.m"), "p"] <- fe.horn.anosim$signif

#--export anosim table
write.csv(anosim.res, paste0(res.dir, "comm_comp_anosim_res.csv"),
          row.names = F)

#--output results
write.csv(permanova.res, paste0(res.dir, "comm_comp_permanova_res_Table3.csv"),
          row.names = F)
