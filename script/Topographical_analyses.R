## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate whether within site topography affected ectomycorrhizal or
## endophytic fungi. Convergent indicates that trees were located in lowlying areas where
## nutrients and water pool, whereas divergent indicates trees that were located in higher
## areas where nutrients and water naturally flow away from.

#=========================================================================================
# Ectomycorrhizal fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------

#--Read in EM data frames
em.stat <- read.csv(paste0(dat.dir,'SCM_EM_per_tree_stats.csv'),
                    as.is = T, header = T)
em.div <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_landscape.csv'), as.is = T,
                        header = T)
abund.data <- read.csv(paste0(dat.dir, 'SCM_EM_root_based_site_x_species_matrix.csv'),
                       as.is = T, header = T)
em.data <- read.csv(paste0(dat.dir, 'SCM_EM_raw.csv'), as.is = T, header = T)

#--Read in FE data frames
fe.stat <- read.csv(paste0(dat.dir,'SCM_FE_per_tree_stats.csv'),
                    as.is = T, header = T)
fe.div <- read.csv(paste0(dat.dir, 'SCM_FE_div_pooled_by_landscape.csv'), as.is = T,
                   header = T)
fe.data <- read.csv(paste0(dat.dir, 'SCM_FE_site_x_species_matrix.csv'), as.is = T,
                    header = T)
fe.tax <- read.csv(paste0(dat.dir, 'SCM_FE_raw.csv'), as.is = T, header = T)
#-----------------------------------------------------------------------------------------
# Ectomycorrhizal fungi - Abundance 
#-----------------------------------------------------------------------------------------
em.abun <- em.stat
#--Uncomment to remove outliers
out <- c("LB028", "LB013", "LB061")
em.abun <- em.abun[which(!em.abun$tree_number %in% out), ]

#--T test of difference of means between the two goups
lm.em.ab <- lm(log(abundance) ~ topography, data = em.abun)
em.anova <- anova(lm.em.ab)
em.anova

#--plot abundance as a function of elevation
em.ab <- ggplot(data = em.abun, aes (x = topography, 
                                     y = abundance)) +
  geom_boxplot() + 
  ylab("Abundance") +
  xlab("Topography") +
  theme_linedraw (base_size = 15) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

rm (lm.em.ab, em.abun, el, out, em.anova)

#-----------------------------------------------------------------------------------------
# Ectomycorrhizal fungi - Diversity 
#-----------------------------------------------------------------------------------------
# << Grouped by within site topography >>-------------------------------------------------
#--change topographic categories names to be capitalized
em.div [em.div$topography == 'convergent', 'topography'] <- 'Convergent'
em.div [em.div$topography == 'divergent', 'topography'] <- 'Divergent'

#--T test of difference of means between the two goups
lm.em.div <- lm (fishers.alpha ~ topography, data = em.div)
em.anova <- anova(lm.em.div)
em.anova

#--Uncomment to remove outlier
#em.div <- em.div[-12,]
#--plot abundance as a function of elevation
em.div <- ggplot(data = em.div, aes (x = topography, 
                                     y = fishers.alpha)) +
  geom_boxplot() + 
  ylab("Fisher's alpha") +
  xlab("Topography") +
  ylim (0,30) +
  theme_linedraw(base_size = 15) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#--export graphs
ggsave ("em_topo_ab_div.jpeg", arrangeGrob (em.ab, em.div, nrow = 1, ncol = 2),
        width = 9, height = 5, device = "jpeg", path = fig.dir, dpi = 1000)

rm (lm.em.div, em.div, em.anova)

#-----------------------------------------------------------------------------------------
# Ectomycorrhizal fungi - ANOSIM and NMDS 
#-----------------------------------------------------------------------------------------
#--Outliers to remove if you wish
em.out <- c("LB025", "LB005", "LB034", "LB038", "LB040", "LB050")
#--List of non-EM otus
root.out <- c("otu58", "otu125", "otu153", "otu44", "otu51", "otu68", "otu122", "otu135")

#<< Morisita index, root tip data >>------------------------------------------------------
#--Remove non-em roots
nmds.ab <- abund.data[!colnames(abund.data) %in% root.out]
#--uncomment to remove rows from otu.table and comm.matrix;
#--outliers are LB005, LB034, LB038, LB040
nmds.ab <- nmds.ab[!(nmds.ab$tree_number %in% em.out),]

#--isolate otu data
comm.matrix <- nmds.ab[5:length(nmds.ab)]
#--make columns numeric
for (i in colnames(comm.matrix)) {
  comm.matrix[[i]] <- as.numeric(comm.matrix[[i]])
}
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]
#--uncomment to remove singletons
#comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "morisita", binary = F)
#--NMDS analysis
horn.abund <- metaMDS(comm.dist.horn, dist = "bray", try = 1000)

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(horn.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep(NA, length(rownames(comm.matrix))),
                        p.group = nmds.ab$group)
color.vec <- sapply(color.vec$p.group, function (x) 
  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
points(horn.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
ordihull(horn.abund, groups = nmds.ab$topography)

#--ANOSIM 
em.horn.anosim <- anosim(comm.matrix, grouping = nmds.ab$topography, distance = "horn")

rm (em.out, em.horn.anosim, comm.matrix, horn.abund, comm.dist.horn, nmds.ab, root.out)

#-----------------------------------------------------------------------------------------
# Ecotmycorrhizal fungi - Taxonomic composition 
#-----------------------------------------------------------------------------------------
#--Make count table of em classes
topo.em <- table(em.data$topography,
                 em.data$taxonomy_class)
#--Remove taxa with less than 1 occurence
rare.em <- "Fungi incertae sedis"
topo.em <- topo.em[,which(!colnames(topo.em) %in% rare.em)]

#--Chi square test
em.chi <- chisq.test(topo.em)

#--Remove taxa with less than 1 occurence
rare.em <- c("Fungi incertae sedis", "Eurotiomycetes", "Sordariomycetes")
em.chisq <- em.data[which(!em.data$taxonomy_class %in% rare.em), ]
#--Bar graph
ggplot(data = em.chisq, 
       aes(x = topography,
           fill = taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  xlab("Topography") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Greys") +
  guides (fill=guide_legend(title=NULL)) +
  theme_linedraw(base_size = 15)

#--export graph
ggsave ("em_topo_tax_comp.jpeg", last_plot(), width = 5, height = 5, device = "jpeg",
        path = fig.dir, dpi = 1000)

rm (em.chi, rare.em, em.chisq, topo.em)

#=========================================================================================
# Endophytic fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Endophytic fungi - Abundance 
#-----------------------------------------------------------------------------------------
#--remove those trees with 0 abundance
fe.abun <- fe.stat[which(!fe.stat$abundance == 0), ]

#--T test of difference of means between the two goups
lm.fe.ab <- lm(logit(abundance) ~ topography, data = fe.abun)
fe.anova <- anova(lm.fe.ab)

#--plot abundance as a function of elevation
fe.ab <- ggplot(data = fe.abun, aes (x = topography, 
                                     y = abundance)) +
  geom_boxplot() + 
  ylab("Abundance") +
  xlab("Topography") +
  ylim(0.0, 0.2) +
  theme_linedraw (base_size = 15) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

rm (lm.fe.ab, fe.anova)

#-----------------------------------------------------------------------------------------
# Endophytic fungi - Diversity 
#-----------------------------------------------------------------------------------------
# << Grouped by pooled by site and topographic group >>-----------------------------------
#--T test of difference of means between the two goups
lm.fe.div <- lm(log(fishers.alpha) ~ topography, data = fe.div)
fe.anova <- anova(lm.fe.div)

#--plot abundance as a function of elevation
fe.div <- ggplot(data = fe.div, aes (x = topography, 
                                            y = fishers.alpha)) +
  geom_boxplot() + 
  ylab("Fisher's alpha") +
  xlab("Topography") +
  ylim(0.0, 8.0) +
  theme_linedraw(base_size = 15) +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#--export graphs
ggsave ("fe_topo_ab_div.jpeg", arrangeGrob (fe.ab, fe.div, nrow = 1, ncol = 2),
        width = 9, height = 5, device = "jpeg", path = fig.dir, dpi = 1000)

rm (lm.fe.div, fe.div, fe.div.pooled, fe.anova)

#-----------------------------------------------------------------------------------------
# Endophytic fungi - ANOSIM and NMDS 
#-----------------------------------------------------------------------------------------

#<< Morisita index >>---------------------------------------------------------------------
nmds.fe <- fe.data
#--Uncomment to remove rows from otu.table and comm.matrix; outliers are LB043, LB044
#fe.out <- c("LB043", "LB044")
#nmds.fe <- nmds.fe[!(nmds.fe$tree_number %in% fe.out),]

#--isolate otu data
comm.matrix <- nmds.fe[2:41]
#--make OTU columns numeric
for(i in colnames(comm.matrix)) {
  comm.matrix[[i]] <- as.numeric(comm.matrix[[i]])
}
#--Remove otu with no occurrences
comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 1)]
#--Uncomment to remove singletons
#comm.matrix <- comm.matrix[which(colSums(comm.matrix) >= 2)]
#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "horn", binary = TRUE)
#--NMDS analysis
morisita.abund <- metaMDS(comm.dist.horn, dist = "bray", try = 1000)

#--Plot NMDS of EM community based on Jaccard index and OTU abundance
plot(morisita.abund, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5, yaxt = "n")
axis (2, at = seq (-0.4, 0.4, by = 0.2), cex.axis = 1.5, las = 2)
# colors for points
color.vec <- data.frame(color = rep(NA, length(rownames(comm.matrix))),
                        p.group = nmds.fe$group)
color.vec <- sapply(color.vec$p.group, function (x) 
  if (x == 'High elev.') {color.vec = 'black'} else {color.vec = 'grey'})
points(morisita.abund, display = "sites", cex = 2, pch = 20,
       col = color.vec,
       bg = color.vec)
ordihull(morisita.abund, groups = nmds.fe$topography)

#--ANOSIM 
fe.horn.anosim <- anosim(comm.matrix, grouping = nmds.fe$topography, distance = "horn")

#--remove temporary data frames, vectors, etc...
rm (comm.matrix, h1, obs, horn.abund, comm.dist.horn, nmds.fe)

#-----------------------------------------------------------------------------------------
# Endophytic fungi - Taxonomic composition 
#-----------------------------------------------------------------------------------------
fe.chisq <- fe.tax
#--Make count table of em classes
topo.fe <- table(fe.chisq$Topography,
                 fe.chisq$Taxonomy)
#--Removed those classes withonly 1 occurence
rare.fe <- c("Lecanoromycetes", "Orbiliomycetes", "Tremellomycetes")
topo.fe <- topo.fe[,which(!colnames(topo.fe) %in% rare.fe)]

#--Chi square test
fe.chi <- chisq.test(topo.fe)

#--Remove rare classes
rare.fe <- c("Lecanoromycetes", "Orbiliomycetes", "Tremellomycetes")
fe.chisq <- fe.chisq[which(!fe.chisq$Taxonomy %in% rare.fe), ]
#--remove unknowns for bar chart
fe.chisq <- fe.chisq[-which(fe.chisq$Taxonomy == "Unknown"), ]
#--Bar graph
fe <- ggplot(data = fe.chisq, 
             aes(x = Topography,
                 fill = Taxonomy)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  xlab("Topography") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Greys") +
  guides (fill=guide_legend(title=NULL)) +
  theme_linedraw(base_size = 12)

#--export graphs
ggsave ("fe_topo_tax_comp.jpeg", width = 5, height = 5,
        device = "jpeg", path = fig.dir, dpi = 1000)

#--Remove temporary data frames
rm (tax.chisq, em.chi, fe.chi, org, fe.tax, fe.table, fe.chisq, rare.fe, rare.em,
    em.chisq, em, fe)