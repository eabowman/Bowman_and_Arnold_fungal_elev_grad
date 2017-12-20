## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the soil and climate patterns as part of a study of these two
## fungal communities along an elevation gradient in the Santa Catalina Mountains

#=========================================================================================
# Analysis of soil community similarity
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frame
#-----------------------------------------------------------------------------------------

#--results table
soil.similarity <- data.frame(soil = c('all.sites', 'high', 'low'), df.1 = NA, df.2 = NA,
                              f.stat = NA, r.sq = NA, p = NA)

#--read in soil data
soil.data <- read.csv(paste0(dat.dir, 'SCM_soil.csv'), as.is = T, header = T)

#<< mean pairwise analysis of between group variation>> -------
soil.dist <- vegdist(soil.data[4:19], method = 'euclidean')

#--PERMANOVA comparing soil similarity across elevations
soil.all <- adonis(soil.dist ~ elevation, data = soil.data)

# << add results to results table >> -----------------------------------------------------

#--across all sites
soil.similarity[soil.similarity$soil == 'all.sites', 'df.1'] <- soil.all$aov.tab$Df[1]
soil.similarity[soil.similarity$soil == 'all.sites', 'df.2'] <- soil.all$aov.tab$Df[2]
soil.similarity[soil.similarity$soil == 'all.sites', 'f.stat'] <-
  soil.all$aov.tab$F.Model[1]
soil.similarity[soil.similarity$soil == 'all.sites', 'r.sq'] <- soil.all$aov.tab$R2[1]
soil.similarity[soil.similarity$soil == 'all.sites', 'p'] <- soil.all$aov.tab$`Pr(>F)`[1]

#--write out results table
write.csv(soil.similarity, paste0(res.dir,
         'Bowman_and_Arnold_soil_similarity_analysis.csv'),
          row.names = F)

#=========================================================================================
# Mantel test and Mantel correlogram
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Mantel test of correlation of environment by geographical distance
#-----------------------------------------------------------------------------------------

#--soil data distance matrix
env.dist <- vegdist(scale(soil.data[4:19]), method = "euclidean")

#--distance matrix of geographical data
geo.dist <- dist(soil.data[22:23])

#--mantel test
mantel <- mantel(env.dist, geo.dist, strata = soil.data$group)
print(mantel)

#--lay out of plots
par(mfrow = c(3,1))

#--Mantel correlogram
man.corr <- mantel.correlog(env.dist, geo.dist)
plot(man.corr)
title(main = 'All sites')

#-----------------------------------------------------------------------------------------
# Exploration of plant communities separately
# soil and geographical distance
#-----------------------------------------------------------------------------------------
#--partition soil.data into plant communities
op.soil <- soil.data[ which(soil.data$p.c == 'o.p'), ]
p.soil <- soil.data[ which(soil.data$p.c == 'p'), ]
pf.soil <- soil.data[ which(soil.data$p.c == 'p.f'), ]
pfmd.soil <- soil.data[ which(soil.data$p.c == 'p.f.md'), ]

#--soil data distance matrix
env.dist.op <- vegdist(scale(op.soil [4:19]), method = "euclidean")
env.dist.p <- vegdist(scale(p.soil [4:19]), method = "euclidean")
env.dist.pf <- vegdist(scale(pf.soil [4:19]), method = "euclidean")
env.dist.pfmd <- vegdist(scale(pfmd.soil [4:19]), method = "euclidean")

#--distance matrix of geographical data
geo.dist.op <- dist(op.soil[22:23])
geo.dist.p <- dist(p.soil[22:23])
geo.dist.pf <- dist(pf.soil[22:23])
geo.dist.pfmd <- dist(pfmd.soil[22:23])

#--mantel: oak pine
mantel.op <- mantel(env.dist.op, geo.dist.op)
print(mantel.op)

man.corr.op <- mantel.correlog(env.dist.op, geo.dist.op)
plot(man.corr.op)
title(main = 'Pine-oak')

#--mantel: Pine
mantel.p <- mantel(env.dist.p, geo.dist.p)
print(mantel.p)

man.corr.p <- mantel.correlog(env.dist.p, geo.dist.p)
plot(man.corr.p)
title(main = 'Pine')

#--mantel: Pine-Doug fir
mantel.pf <- mantel(env.dist.pf, geo.dist.pf)
print(mantel.pf)

man.corr.pf <- mantel.correlog(env.dist.pf, geo.dist.pf)
plot(man.corr.pf)
title(main = 'Pine-Doug fir')

#--mantel: Pine-Doug fir-mixed deciduous
mantel.pfmd <- mantel(env.dist.pfmd, geo.dist.pfmd)
print(mantel.pfmd)

man.corr.pfmd <- mantel.correlog(env.dist.pfmd, geo.dist.pfmd)
plot(man.corr.pfmd)
title(main = 'Pine-Doug fir-mixed deciduous')

#-----------------------------------------------------------------------------------------
# Parital mantel test of correlation of EM community by environment and geo. distance
#-----------------------------------------------------------------------------------------
#--data-----
#--data for jaccard based dissimilarity matrix
jac.perm <- read.csv(paste0(dat.dir,
                            'SCM_EM_otu_site_x_species_clustered_bytopoandelev.csv'),
                     as.is = T, header = T)

#--data for Morisita-horn based dissimilarity index
mor.perm <- read.csv(paste0(dat.dir,
                            'SCM_EM_root_site_x_species_clustered_bytopoandelev.csv'),
                     as.is = T, header = T)
#--FE site x species matrix based on OTU count for northing and easting data
fe.data <- read.csv(paste0(dat.dir,'SCM_FE_site_x_species_matrix.csv'),
                    as.is = T, header = T)

#--change plant community from categorical data to numeric data for distance matrices
jac.perm[jac.perm$Plant.community == 'Pinus-Quercus', 'Plant.community'] <- 1
jac.perm[jac.perm$Plant.community == 'Pinus', 'Plant.community'] <- 2
jac.perm[jac.perm$Plant.community == 'Pinus-Pseudotsuga', 'Plant.community'] <- 3
jac.perm[jac.perm$Plant.community == 'Pinus-Pseudotsuga-Mixed', 'Plant.community'] <- 4
mor.perm[mor.perm$Plant.community == 'Pinus-Quercus', 'Plant.community'] <- 1
mor.perm[mor.perm$Plant.community == 'Pinus', 'Plant.community'] <- 2
mor.perm[mor.perm$Plant.community == 'Pinus-Pseudotsuga', 'Plant.community'] <- 3
mor.perm[mor.perm$Plant.community == 'Pinus-Pseudotsuga-Mixed', 'Plant.community'] <- 4

#--add northing and easting to data frame
for(i in unique(jac.perm$Elevation)){
  jac.perm[jac.perm$Elevation == i, 'northing'] <-
    unique(fe.data[fe.data$elevation == i, 'northing'])
  jac.perm[jac.perm$Elevation == i, 'easting'] <-
    unique(fe.data[fe.data$elevation == i, 'easting'])
  mor.perm[mor.perm$Elevation == i, 'northing'] <-
    unique(fe.data[fe.data$elevation == i, 'northing'])
  mor.perm[mor.perm$Elevation == i, 'easting'] <-
    unique(fe.data[fe.data$elevation == i, 'easting'])
}

#--create distance matrices
env.dist <- vegdist(jac.perm[c('Plant.community','Phosphate',
                               'Average.temperature.warm.quarter')], method = 'euclidean')
geo.dist <- vegdist(jac.perm[c('northing','easting')], method = 'euclidean')
mor.comm.dist <- vegdist(mor.perm[8:130], method = 'morisita', binary = F)
jac.comm.dist <- vegdist(jac.perm[8:28], method = 'jaccard', binary = T)

#<< without elevational high and low groupings >> ----------------------------------------
#--Morisita based distance matrix
mantel.partial(mor.comm.dist, geo.dist, env.dist)

#--Jaccard based distance matrix
mantel.partial(jac.comm.dist, geo.dist, env.dist)

#<< with elevational high and low groupings >> ----------------------------------------
#--Morisita based distance matrix
mantel.partial(mor.comm.dist, geo.dist, env.dist, strata = mor.perm$group)

#--Jaccard based distance matrix
mantel.partial(jac.comm.dist, geo.dist, env.dist, strata = jac.perm$Group)

#-----------------------------------------------------------------------------------------
# Export distance matrices
#-----------------------------------------------------------------------------------------

#--Soil
soil.dist <- as.data.frame(t(combn(as.character(soil.data$tree_number),2)))
soil.dist$distance <- as.numeric(env.dist)
colnames(soil.dist) <- c('c1', 'c2', 'distance')
write.csv(soil.dist, paste0(res.dir, "Appendix S9.csv"))

#--geographical distance
dist.dist <- as.data.frame(t(combn(as.character (soil.data$tree_number),2)))
dist.dist$distance <- as.numeric(geo.dist)
colnames(dist.dist) <- c('c1', 'c2', 'distance')
write.csv(dist.dist, paste0(res.dir, "Appendix S10.csv"))

#=========================================================================================
# Soil as a function of elevation
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Linear regression and Anova/Kruskal test
#-----------------------------------------------------------------------------------------

#--anova and kruskal
soil.traits <- c("ph.su", "EC.ds.m")
soil.nonpara <- c('ca.ppm', 'mg.ppm', 'na.ppm', 'zn.ppm', 'cu.ppm', 'so4.s.ppm')

#--all kruskal except fe
tree.1 <- c('k.ppm', 'b.ppm', 'cec.meq.100g')
tree.37 <- 'po4.p.ppm'
tree.13 <- 'fe.ppm'
tree.20 <- 'mn.ppm'
tree.34 <- 'ni.ppm'

#<< ANOVA >>------------------------------------------------------------------------------
#--Make data frame for restuls
anova.soil <- data.frame(c(soil.traits, soil.nonpara,
                            tree.1, tree.34, tree.20, tree.13, tree.37))
colnames(anova.soil) <- "soil"
#--for loop to do anova and add results to data frame
for(i in soil.traits) {
  #--normally dist. data
  lm.i <- lm(soil.data[[i]] ~ soil.data$elevation)
  anova.i <- anova (lm.i)
  anova.soil[which (anova.soil$soil == i), "anova.f"] <- anova.i$`F value`[1]
  anova.soil[which (anova.soil$soil == i), "anova.p"] <- anova.i$`Pr(>F)`[1]
}

#--loop for non-parametric analyses
for(i in soil.nonpara) {
  #--log transformed and normally dist. data
  kruskal.i <- kruskal.test(soil.data[[i]] ~ soil.data$elevation)
  anova.soil[which (anova.soil$soil == i), "kruskal.f"] <- kruskal.i$statistic
  anova.soil[which (anova.soil$soil == i), "kruskal.p"] <- kruskal.i$p.value
}

#--loop for those with outliers
for(i in tree.1) {
  soil.data <- soil.data [which (!soil.data$tree_number == 'LB001'), ]
  #--log transformed and normally dist. data
  kruskal.i <- kruskal.test(soil.data[[i]] ~ soil.data$elevation)
  anova.soil [which (anova.soil$soil == i), "kruskal.f"] <- kruskal.i$statistic
  anova.soil [which (anova.soil$soil == i), "kruskal.p"] <- kruskal.i$p.value
}

#--for phosphate
po4.p <- soil.data [which (!soil.data$tree_number == 'LB037'), ]
#--log transformed and normally dist. data
kruskal.i <- kruskal.test(soil.data$po4.p.ppm ~ soil.data$elevation)
anova.soil[which(anova.soil$soil == 'po4.p.ppm'), "kruskal.f"] <- kruskal.i$statistic
anova.soil[which(anova.soil$soil == 'po4.p.ppm'), "kruskal.p"] <- kruskal.i$p.value
#--for iron
fe <- soil.data[which(!soil.data$tree_number == 'LB013'), ]
#--log transformed and normally dist. data
kruskal.i <- kruskal.test(soil.data$fe.ppm ~ soil.data$elevation)
anova.soil[which (anova.soil$soil == 'fe.ppm'), "kruskal.f"] <- kruskal.i$statistic
anova.soil[which (anova.soil$soil == 'fe.ppm'), "kruskal.p"] <- kruskal.i$p.value
#--for manganese
mn <- soil.data[which(!soil.data$tree_number == 'LB020'), ]
#--log transformed and normally dist. data
kruskal.i <- kruskal.test(soil.data$mn.ppm ~ soil.data$elevation)
anova.soil[which(anova.soil$soil == 'mn.ppm'), "kruskal.f"] <- kruskal.i$statistic
anova.soil[which(anova.soil$soil == 'mn.ppm'), "kruskal.p"] <- kruskal.i$p.value
#--for nickel
ni <- soil.data[which(!soil.data$tree_number == 'LB034'), ]
#--log transformed and normally dist. data
kruskal.i <- kruskal.test(soil.data$ni.ppm ~ soil.data$elevation)
anova.soil[which (anova.soil$soil == 'ni.ppm'), "kruskal.f"] <- kruskal.i$statistic
anova.soil[which (anova.soil$soil == 'ni.ppm'), "kruskal.p"] <- kruskal.i$p.value

#--write csv file of linear regressions
write.csv(anova.soil, paste0(res.dir, "Appendix S8.csv"),
          row.names = FALSE)

#<< Graph as function of elevation >>-----------------------------------------------------

#--Export as jpeg
jpeg(filename = paste0 (fig.dir,"Appendix S7.jpeg"),
     width = 1000, height = 800)
#--Combine into one figure
par(mfrow = c(4,4), mar = c(2,4,4,4), oma = c(4,4,4,4))

#--Phosphate
plot(po4.p$elevation, po4.p$po4.p.ppm, pch = 19,
      ylab = "Phosphate (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 50, by = 10), cex.axis = 1.5, las = 2)
abline(m <- lm (po4.p.ppm ~ elevation, data = po4.p))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Phosphate (ppm)", line = 4, cex = 1)

#--pH
plot(soil.data$elevation, soil.data$ph.su, pch = 19,
      ylab = "pH (SU)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 8, by = 1.0), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "pH (SU)", line = 4, cex = 1)

#--EC
plot(soil.data$elevation, soil.data$EC.ds.m, pch = 19,
      cex.axis = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0.00, 0.3, by = 0.05), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "El. conductivity (ds/m)", line = 4, cex = 1)

#--Calcium
plot(soil.data$elevation, soil.data$ca.ppm, pch = 19,
      ylab = "Calcium (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 4000, by = 1000), cex.axis = 1.5, las = 2)
abline(m <- lm(ca.ppm ~ elevation, data = soil.data))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Calcium (ppm)", line = 4, cex = 1)

#--Magnesium
plot(soil.data$elevation, soil.data$mg.ppm, pch = 19,
      ylab = "Magnesium (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 300, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(mg.ppm ~ elevation, data = soil.data))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Magnesium (ppm)", line = 3.5, cex = 1)

#--Sodium
plot(soil.data$elevation, soil.data$na.ppm, pch = 19,
      ylab = "Sodium (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 30, by = 5), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Sodium (ppm)", line = 4, cex = 1)

#--make soil data frame minus LB001 for potassium, boron, and cec
tree.1 <- soil.data[which(!soil.data$tree_number == 'LB001'), ]
#--Potassium
plot(tree.1$elevation, tree.1$k.ppm, pch = 19,
      ylab = "Potassium (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 300, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(k.ppm ~ elevation, data = tree.1))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Potassium (ppm)", line = 4, cex = 1)

#--Zinc
plot(soil.data$elevation, soil.data$zn.ppm, pch = 19,
      ylab = "Zinc (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 30, by = 5), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Zinc (ppm)", line = 4, cex = 1)

#--Iron
plot(fe$elevation, fe$fe.ppm, pch = 19,
      ylab = "Iron (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 200, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(fe.ppm ~ elevation, data = fe))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Iron (ppm)", line = 4, cex = 1)

#--Manganese
plot(mn$elevation, mn$mn.ppm, pch = 19,
      ylab = "Manganese (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 60, by = 10), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Manganese (ppm)", line = 4, cex = 1)

#--Copper
plot(soil.data$elevation, soil.data$cu.ppm, pch = 19,
      ylab = "Copper (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 10, by = 2), cex.axis = 1.5, las = 2)
abline(m <- lm(cu.ppm ~ elevation, data = soil.data))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Copper (ppm)", line = 4, cex = 1)

#--Nickel
plot(ni$elevation, ni$ni.ppm, pch = 19,
      ylab = "Nickel (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 0.7, by = 0.2), cex.axis = 1.5, las = 2)
abline(m <- lm(ni.ppm ~ elevation, data = ni))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Nickel (ppm)", line = 4, cex = 1)

#--Sulfate
plot(soil.data$elevation, soil.data$so4.s.ppm, pch = 19,
      ylab = "Sulfate (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 30, by = 5), cex.axis = 1.5, las = 2)
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Sulfate (ppm)", line = 4, cex = 1)

#--Boron
plot(tree.1$elevation, tree.1$b.ppm, pch = 19,
      ylab = "Boron (ppm)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 1.0, by = 0.2), cex.axis = 1.5, las = 2)
abline(m <- lm(b.ppm ~ elevation, data = tree.1))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "Boron (ppm)", line = 4, cex = 1)

#--CEC
plot(tree.1$elevation, tree.1$cec.meq.100g, pch = 19,
      ylab = "CEC (meq/100g)", xlab = "Elevation (m)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 16, by = 4), cex.axis = 1.5, las = 2)
abline(m <- lm(cec.meq.100g ~ elevation, data = tree.1))
mtext(side = 1, text = "Elevation (m)", line = 3, cex = 1)
mtext(side = 2, text = "CEC (meq/100g)", line = 4, cex = 1)

#--close
dev.off()

#=========================================================================================
# Correlationo of Phosphorus with other soil components positively correlated with elev.
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Linear regression and Anova/Kruskal test: Outliers removed
#-----------------------------------------------------------------------------------------
#--kruskal
soil.tran <- c('ca.ppm', 'mg.ppm', 'cu.ppm')
#--all kruskal except fe
tree.1 <- c('k.ppm', 'b.ppm', 'cec.meq.100g')
tree.13 <- 'fe.ppm'
tree.34 <- 'ni.ppm'

#--Remove LB037 from Phosphate
soil.data <- soil.data [-which (soil.data$tree_number == 'LB037'),]

#<< ANOVA >>------------------------------------------------------------------------------
#--Make data frame for results
anova.phos <- data.frame(c(soil.nonpara, tree.1, tree.34, tree.13))
colnames(anova.phos) <- "soil"
#--log transformed soil data
for(i in soil.tran) {
  lm.i <- lm(log(soil.data[[i]]) ~ log(soil.data$po4.p.ppm))
  anova.i <- anova(lm.i)
  anova.phos[which(anova.phos$soil == i), "anova.f"] <- anova.i$`F value`[1]
  anova.phos[which(anova.phos$soil == i), "anova.p"] <- anova.i$`Pr(>F)`[1]
}

#--loop for those with outliers
for(i in tree.1) {
  soil.data <- soil.data[which(!soil.data$tree_number == 'LB001'), ]
  #--log transformed data
  lm.i <- lm(log(soil.data[[i]]) ~ log(soil.data$po4.p.ppm))
  anova.i <- anova (lm.i)
  anova.phos[which(anova.phos$soil == i), "anova.f"] <- anova.i$`F value`[1]
  anova.phos[which(anova.phos$soil == i), "anova.p"] <- anova.i$`Pr(>F)`[1]
}

#--for iron
fe <- soil.data[which(!soil.data$tree_number == 'LB013'), ]
#--log transformed data
lm.fe <- lm(log(fe$fe.ppm) ~ log(fe$po4.p.ppm))
anova.fe <- anova(lm.fe)
anova.phos[which(anova.phos$soil == 'fe.ppm'), "anova.f"] <- anova.fe$`F value`[1]
anova.phos[which(anova.phos$soil == 'fe.ppm'), "anova.p"] <- anova.fe$`Pr(>F)`[1]

#--for nickel
ni <- soil.data[which(!soil.data$tree_number == 'LB034'), ]
#--log transformed and normally dist. data
kruskal.i <- kruskal.test(soil.data$ni.ppm ~ soil.data$po4.p.ppm)
anova.phos[which(anova.phos$soil == 'ni.ppm'), "kruskal.f"] <- kruskal.i$statistic
anova.phos[which(anova.phos$soil == 'ni.ppm'), "kruskal.p"] <- kruskal.i$p.value
#--log transformed data
lm.ni <- lm(log(ni$ni.ppm) ~ log(ni$po4.p.ppm))
anova.ni <- anova(lm.ni)
anova.phos[which(anova.phos$soil == 'ni.ppm'), "anova.f"] <- anova.ni$`F value`[1]
anova.phos[which(anova.phos$soil == 'ni.ppm'), "anova.p"] <- anova.ni$`Pr(>F)`[1]

#--write csv file of linear regressions
write.csv(anova.phos, paste0 (res.dir, "Appendix S14.csv"), row.names = FALSE)

#<< Graph as function of phosphorus >>----------------------------------------------------
#--Export as jpeg
jpeg(filename = paste0(fig.dir,"Appendix S13.jpeg"),
     width = 1000, height = 800)
#--Combine into one figure
par(mfrow = c(4,2), mar = c(2,8,4,4), oma = c(2,10,2,10))

#--Calcium
plot(soil.data$po4.p.ppm, soil.data$ca.ppm, pch = 19,
      ylab = "Calcium (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 4000, by = 1000), cex.axis = 1.5, las = 2)
abline(m <- lm(ca.ppm ~ po4.p.ppm, data = soil.data))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Calcium (ppm)", line = 4, cex = 1)

#--Magnesium
plot(soil.data$po4.p.ppm, soil.data$mg.ppm, pch = 19,
      ylab = "Magnesium (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 300, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(mg.ppm ~ po4.p.ppm, data = soil.data))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Magnesium (ppm)", line = 3.5, cex = 1)

#--make soil data frame minus LB001 for potassium, boron, and cec
tree.1 <- soil.data [which (!soil.data$tree_number == 'LB001'), ]
#--Potassium
plot(tree.1$po4.p.ppm, tree.1$k.ppm, pch = 19,
      ylab = "Potassium (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 300, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(k.ppm ~ po4.p.ppm, data = tree.1))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Potassium (ppm)", line = 4, cex = 1)

#--Iron
plot(fe$po4.p.ppm, fe$fe.ppm, pch = 19,
      ylab = "Iron (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 200, by = 50), cex.axis = 1.5, las = 2)
abline(m <- lm(fe.ppm ~ po4.p.ppm, data = fe))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Iron (ppm)", line = 4, cex = 1)

#--Copper
plot(soil.data$po4.p.ppm, soil.data$cu.ppm, pch = 19,
      ylab = "Copper (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = FALSE)
axis(2, at = seq (0, 10, by = 2), cex.axis = 1.5, las = 2)
abline(m <- lm(cu.ppm ~ po4.p.ppm, data = soil.data))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Copper (ppm)", line = 4, cex = 1)

#--Nickel
plot(ni$po4.p.ppm, ni$ni.ppm, pch = 19,
      ylab = "Nickel (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 0.7, by = 0.2), cex.axis = 1.5, las = 2)
abline(m <- lm(ni.ppm ~ po4.p.ppm, data = ni))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Nickel (ppm)", line = 4, cex = 1)

#--Boron
plot(tree.1$po4.p.ppm, tree.1$b.ppm, pch = 19,
      ylab = "Boron (ppm)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 1.0, by = 0.2), cex.axis = 1.5, las = 2)
abline(m <- lm(b.ppm ~ po4.p.ppm, data = tree.1))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "Boron (ppm)", line = 4, cex = 1)

#--CEC
plot(tree.1$po4.p.ppm, tree.1$cec.meq.100g, pch = 19,
      ylab = "CEC (meq/100g)", xlab = "Phosphate (ppm)",
      cex.axis = 1.5, cex.lab = 1.5, yaxt = "n", ann = F)
axis(2, at = seq (0, 16, by = 4), cex.axis = 1.5, las = 2)
abline(m <- lm(cec.meq.100g ~ po4.p.ppm, data = tree.1))
mtext(side = 1, text = "Phosphate (ppm)", line = 3, cex = 1)
mtext(side = 2, text = "CEC (meq/100g)", line = 4, cex = 1)

#--close
dev.off()

#=========================================================================================
# Correlation of climate data 
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Graphs
#-----------------------------------------------------------------------------------------
#<< Graph climate data as function of other climate traits >>-----------------------------
#--Export as jpeg
jpeg(filename = paste0 (fig.dir,"Appendix S8.jpeg"),
     width = 800, height = 800)
#--Combine into one figure
par(mfrow = c(2,1), mar = c(7,7,5,5), oma = c(0,8,0,8))

#--Avg. high temperature and precipitation
lm.prec <- lm(soil.data$annual.prec.cm ~ soil.data$average.temp.warm.quarter)
anova.prec <- anova(lm.prec)
# graph
plot(soil.data$average.temp.warm.quarter, soil.data$annual.prec.cm, pch = 19,
    ylab = "Avg. annual precipitation (cm)",
    xlab = "",
    cex.axis = 1.5, cex.lab = 1.5)
mtext(side = 1, text = "Avg. temperature in the \n warmest quarter (°C)",
      line = 4, cex = 1.5)
axis(2, at = seq (0, 4000, by = 1000), cex.axis = 1.5, las = 2)
abline(m <- lm(annual.prec.cm ~ average.temp.warm.quarter, data = soil.data))

#--Avg. high temperature and avg. snow
lm.cold <- lm(soil.data$average.temp.cold.quarter ~ soil.data$average.temp.warm.quarter)
anova.cold <- anova(lm.cold)
# graph
plot(soil.data$average.temp.warm.quarter, soil.data$average.temp.cold.quarter, pch = 19,
    ylab = "Avg. temperature in the \n coldest quarter (°C)",
    xlab = "",
    cex.axis = 1.5, cex.lab = 1.5)
mtext(side = 1, text = "Avg. temperature in the \n warmest quarter (°C)",
      line = 4, cex = 1.5)

axis(2, at = seq (0, 4000, by = 1000), cex.axis = 1.5, las = 2)
abline(m <- lm(average.temp.cold.quarter ~ average.temp.warm.quarter, data = soil.data))

#--close
dev.off()
