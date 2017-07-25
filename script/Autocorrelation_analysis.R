## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the soil and climate patterns as part of a study of these two
## fungal communities along an elevation gradient in the Santa Catalina Mountains

#=========================================================================================
# Analysis of autocorrelation of fungal communities
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frame
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
# EM data
#-----------------------------------------------------------------------------------------
#--make soil data numeric
env <- c("elevation","ph.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm", "k.ppm","zn.ppm",
         "fe.ppm","mn.ppm","cu.ppm","ni.ppm", "po4.p.ppm","so4.s.ppm","b.ppm","esp",
         "cec.meq.100g","northing","easting")
for (i in env) {
  nmds.rooto[[i]] <- as.numeric(nmds.rooto[[i]])
}

#<< Create distance matrices >> ----------------------------------------------------------
#--geographical distance using easting and northing columns
geo.dist <- vegdist(otu.data[178:179], method = 'euclidean')
#--low elevation
geo.dist.low <- vegdist(otu.data[otu.data$group == 'Low elev.',178:179],
                        method = 'euclidean')
#--high elevation
geo.dist.high <- vegdist(otu.data[otu.data$group == 'High elev.',178:179],
                         method = 'euclidean')

#--Jaccard index, OTU data
#--Comment to include outliers 
#otu.data <- otu.data[!(otu.data$tree_number %in% em.out),]
#--isolate otu data
comm.matrix <- otu.data[2:157]
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(!(colSums(comm.matrix) == 0))]
#--comment to add singletons
#comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--low elevation
#--isolate otu.data
comm.matrix.low <- otu.data[otu.data$group == 'Low elev.', 2:157]
#--Remove otu with no occurrences
comm.matrix.low <- comm.matrix.low[which(!(colSums(comm.matrix.low) == 0))]
#--distance matrix using jaccard index
comm.dist.jlow <- vegdist(comm.matrix.low, method = 'jaccard', binary = T)

#--high elevation
#--isolate otu.data
comm.matrix.high <- otu.data[otu.data$group == 'High elev.', 2:157]
#--Remove otu with no occurrences
comm.matrix.high <- comm.matrix.high[which(!(colSums(comm.matrix.high) == 0))]
#--distance matrix using jaccard index
comm.dist.jhigh <- vegdist(comm.matrix.high, method = 'jaccard', binary = T)

#--soil
soil.dist <- vegdist(otu.data[c(160:174)], method = 'euclidean')
#--low elevation
soil.dist.low <- vegdist(otu.data[otu.data$group == 'Low elev.',160:174])
#--high elevation
soil.dist.high <- vegdist(otu.data[otu.data$group == 'High elev.',160:174])

#<< EM x geographical distance >> --------------------------------------------------------
#--across all sites
geo.em <- mantel(geo.dist, comm.dist.jaccard)

#--only low elevation sites
geo.emlow <- mantel(geo.dist.low, comm.dist.jlow)

#--only high elevation sites
geo.emhigh <- mantel(geo.dist.high, comm.dist.jhigh)

#<< EM x soil >> -------------------------------------------------------------------------
#--across all sites
soil.em <- mantel(soil.dist, comm.dist.jaccard)

#--only low elevation sites
soil.emlow <- mantel(soil.dist.low, comm.dist.jlow)

#--only high elevation sites
soil.emhigh <- mantel(soil.dist.high, comm.dist.jhigh)

#<< geographical distance x soil >> ------------------------------------------------------
#--across all sites
soil.geo <- mantel(soil.dist, geo.dist)

#--only low elevation sites
soil.geolow <- mantel(soil.dist.low, geo.dist.low)

#--only high elevation sites
soil.geohigh <- mantel(soil.dist.high, geo.dist.high)

#<< parital mantel of geo dist x em with soil variation removed >> -----------------------
#--across all sites
em.geo.soil <- mantel.partial(comm.dist.jaccard, geo.dist, soil.dist)

#--only low elevation sites
lem.geo.soil <- mantel.partial(comm.dist.jlow, geo.dist.low, soil.dist.low)

#--only high elevation sites
hem.geo.soil <- mantel.partial(comm.dist.jhigh, geo.dist.high, soil.dist.high)

#=========================================================================================
# Analysis of distance across all sites and within elevation groups
#=========================================================================================
#--isolate elevation and coordinate data
distance <- data.frame(elevation = unique(otu.data$elevation),
                       northing = unique(otu.data$northing),
                       easting = unique(otu.data$easting))
distance$group <- ifelse(distance$elevation < 2300,'Low elev.','High elev.')

#--make distance matrices
all.dist <- vegdist(distance[2:3], method = 'euclidean')
high.dist <- vegdist(distance[distance$elevation > 2300, 2:3], method = 'euclidean')
low.dist <- vegdist(distance[distance$elevation < 2300, 2:3], method = 'euclidean')

#--PERMANOVA across all sites
all <- adonis(all.dist ~ elevation, data = distance)

#--PERMANOVA across High elev. sites
high <- adonis(high.dist ~ elevation, data = distance[which(distance$elevation > 2300),])

#--PERMANOVA across Low elev. sites
low <- adonis(low.dist ~ elevation, data = distance[distance$elevation < 2300,])

#=========================================================================================
# PERMANOVA taking into account distance
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Outlier trees and non-EM OTUs for removal 
#-----------------------------------------------------------------------------------------
#--Outliers to remove for EM fungi
em.out <- c("LB025", "LB005", "LB034", "LB038", "LB040", "LB050")
#--Remove non-EM otus
root.out <- c("otu58", "otu125", "otu153", "otu44", "otu51", "otu68", "otu122", "otu135")

#--Remove non-em fungi only from dataframe
nmds.ab <- abund.data[!(abund.data$tree_number %in% root.out),]
#--Comment to include outliers and non-em fungi
#--remove rows from otu.table and comm.matrix; outliers are LB005, LB034, LB038, LB040
nmds.ab <- nmds.ab [!(nmds.ab$tree_number %in% em.out),]

#--make soil data numeric
env <- c("elevation","ph.su","EC.ds.m","ca.ppm","mg.ppm","na.ppm", "k.ppm","zn.ppm",
         "fe.ppm","mn.ppm","cu.ppm","ni.ppm", "po4.p.ppm","so4.s.ppm","b.ppm","esp",
         "cec.meq.100g","northing","easting")
for (i in env) {
  nmds.ab[[i]] <- as.numeric(nmds.ab [[i]])
}

#-----------------------------------------------------------------------------------------
# Create distance matrix of community data
#-----------------------------------------------------------------------------------------
#--isolate otu data
comm.matrix <- nmds.ab[2:157]
#--Remove otu with no occurrences
#comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]
#--comment to add singletons
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--remove OTUs with less than 5 occurrences
#comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 5)]
#--distance matrix using morisita index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)

#--PERMANOVA: Terms added sequentially
adonis(formula = comm.dist.horn ~ average.temp.warm.quarter
       + po4.p.ppm
       + plant.comm
       + northing + easting,
       data = nmds.ab, permutations = 1000, strata = nmds.ab$group)

#--PERMANOVA: Combination of variables
adonis(formula = comm.dist.horn ~ average.temp.warm.quarter
       * po4.p.ppm
       * plant.comm
       * northing * easting,
       data = nmds.ab, permutations = 1000, strata = nmds.ab$group)

#--PERMANOVA: just distance
#--PERMANOVA: Combination of variables
adonis(formula = comm.dist.horn ~ northing + easting,
       data = nmds.ab, permutations = 1000, strata = nmds.ab$group)
