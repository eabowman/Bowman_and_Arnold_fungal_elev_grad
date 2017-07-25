## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate whether the two fungal communities were sampled deeply enough.

#=========================================================================================
# Species accumulation curves
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------

#--read in EM data
otu.data <- read.csv(paste0(dat.dir, 'SCM_EM_otu_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)
#--Read in FE data
fe.data <- read.csv(paste0(dat.dir, 'SCM_FE_site_x_species_matrix.csv'), as.is = T,
                    header = T)

#-----------------------------------------------------------------------------------------
# EM species accumulation curve based on OTU count
#-----------------------------------------------------------------------------------------
#<< With singeltons included, Overall EM >>-----------------------------------------------
#--Isolate OTU data
ov.em <- otu.data[5:length(otu.data)]
#--remove columns equal to 0, no species occurrence 
ov.em <- ov.em[which(colSums(ov.em) != 0)]
#--abundance for each otu, double check for discrepancies
summarise_each(ov.em,funs(sum))

#<< With singeltons included, Low elevation group, EM >>----------------------------------
#--Isolate OTU data, low elevation
le.em <- otu.data[which(otu.data$group == "Low elev."), ]
le.em <- le.em[5:length(le.em)]
#--remove columns equal to 0, no species occurence
le.em <- le.em[which(colSums(le.em) != 0)]

#<< With singeltons included, High elevation group, EM >>---------------------------------
#--Isolate OTU data, high elevation
he.em <- otu.data[which (otu.data$group == "High elev."), ]
he.em <- he.em[5:length(le.em)]
#--remove columns equal to 0, no species occurence
he.em <- he.em[which(colSums(he.em) != 0)]

#<< With singeltons removed, Overall EM >>------------------------------------------------
#--Using previous dataframe, remove singletons
ovns.em <- ov.em[which(colSums(ov.em) != 1)]
#--abundance for each otu, double check for discrepancies
colSums(ovns.em)

#<< With singeltons removed, Low elevation group, EM >>-----------------------------------
#--Using previous dataframe, remove singletons
lens.em <- le.em[which(colSums(le.em) != 1)]
#--abundance for each otu, double check for discrepancies
colSums(lens.em)

#<< With singeltons removed, High elevation group, EM >>----------------------------------
#--Using previous dataframe, remove singletons
hens.em <- he.em[which(colSums(he.em) != 1)]
#--abundance for each otu, double check for discrepancies
colSums (hens.em)

#-----------------------------------------------------------------------------------------
# FE species accumulation curve based on OTU count
#-----------------------------------------------------------------------------------------
#<< With singeltons included, Overall FE >>-----------------------------------------------
#--Isolate OTU data
ov.fe <- fe.data [2:41]
#--remove columns equal to 0, no species occurence
ov.fe <- ov.fe[which(colSums (ov.fe) != 0)]

#<< With singeltons included, Low elevation FE >>-----------------------------------------
#--Isolate OTU data
le.fe <- fe.data[which(fe.data$group == "Low elev."), ]
le.fe <- le.fe[2:41]
#--remove columns equal to 0, no species occurence
le.fe <- le.fe[which(colSums(le.fe) != 0)]

#<< With singeltons included, High elevation FE >>-----------------------------------------
#--Isolate OTU data
he.fe <- fe.data[which(fe.data$group == "High elev."), ]
he.fe <- he.fe[2:41]
#--remove columns equal to 0, no species occurence
he.fe <- he.fe[which(colSums(he.fe) != 0)]

#<< With singeltons removed, Overall FE >>------------------------------------------------
#--Using previous dataframe, remove singletons
ovns.fe <- ov.fe[which(colSums(ov.fe) != 1)]

#<< With singeltons removed, Low elevation FE >>------------------------------------------
#--Using previous dataframe, remove singletons
lens.fe <- le.fe[which(colSums(le.fe) != 1)]

#<< With singeltons removed, High elevation FE >>------------------------------------------
#--Using previous dataframe, remove singletons
hens.fe <- he.fe[which(colSums(he.fe) != 1)]

#-----------------------------------------------------------------------------------------
# EM and FE species accumulation curve based on OTU count, combined plots
#-----------------------------------------------------------------------------------------
#--Export as jpeg
jpeg(filename = paste0 (fig.dir,"Figure 3.jpeg"),
    width = 1200, height = 1200, quality = 100)
#--Combine into one figure
par(mfrow = c(3,2), "mar"=c(6, 5, 5, 3))

#<< Combine overall EM plots, with singletons and without >>------------------------------
plot(specaccum(ov.em, sample = min(rowSums(ov.em),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,150), yaxt = "n")
axis(2, at = seq(0, 150, by = 50), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(ovns.em, sample = min(rowSums(ovns.em),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,150))

#<< Combine overall FE plots, with singletons and without >>------------------------------
plot(specaccum(ov.fe, sample = min(rowSums (ov.fe),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,40), yaxt = "n")
axis(2, at = seq(0, 40, by = 10), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(ovns.fe, sample = min(rowSums (ovns.fe),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,40))

#<< Combine low elev. EM plots, with singletons and without >>----------------------------
plot(specaccum(le.em, sample = min(rowSums(le.em),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,100), yaxt = "n")
axis(2, at = seq(0, 100, by = 20), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(lens.em, sample = min(rowSums(lens.em),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,100))

#<< Combine low elev. FE plots, with singletons and without >>----------------------------
plot(specaccum(le.fe, sample = min(rowSums (le.fe),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,40), yaxt = "n")
axis(2, at = seq (0, 40, by = 10), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(lens.fe, sample = min(rowSums(lens.fe),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,40))

#<< Combine high elev. EM plots, with singletons and without >>----------------------------
plot(specaccum(he.em, sample = min(rowSums(he.em), permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,100), yaxt = "n")
axis(2, at = seq (0, 100, by = 20), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(hens.em, sample = min(rowSums (hens.em),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,100))

#<< Combine high elev. FE plots, with singletons and without >>----------------------------
plot(specaccum(he.fe, sample = min(rowSums(he.fe),  permutations = 999)),
     xlab = "Samples", ylab = "OTUs", cex.lab = 2, cex.axis = 2,
     ylim = c(0,40), yaxt = "n")
axis(2, at = seq(0, 40, by = 10), cex.axis = 2, las = 2)
par(new=TRUE)
plot(specaccum(hens.fe, sample = min(rowSums (hens.fe),  permutations = 999)),
     axes = FALSE, xlab = "", ylab = "", col = "darkgrey", ylim = c(0,40))

#--Close export of figure
dev.off()
