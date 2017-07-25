## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the diversity of ectomycorrhizal and endophytic abundance as a
## function of elevation as part of a study of these two fungal communities along an
## elevation gradient in the Santa Catalina Mountains

#=========================================================================================
# Ectomycorrhizal fungi - Fisher's alpha diversity 
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------------------------------
#<< Read in pooled diversity >> ----------------------------------------------------------
#--diversity pooled by site and topography
em.div <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_landscape.csv'),
                   header = T, as.is = T)
#--diversity poooled only by site (for mean and sd)
em.div.site <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_site.csv'),
                        as.is = T, header = T)
#--results table
div.fa <- data.frame(symbiont = c('EM', 'FE'), df.1 = NA, df.2 = NA, f.stat = NA, p = NA,
                     mean = NA, sd = NA)

#-----------------------------------------------------------------------------------------
# Anova of EM diversity and box plot
#-----------------------------------------------------------------------------------------
#--Anova of difference between means of the two goups (high and low elevatio group)
lm.em.div <- lm(log(fishers.alpha) ~ group, data = em.div)
em.anova <- anova(lm.em.div)
em.anova
#--add data to results table
div.fa[div.fa$symbiont == 'EM', 'df.1'] <- em.anova$Df[1]
div.fa[div.fa$symbiont == 'EM', 'df.2'] <- em.anova$Df[2]
div.fa[div.fa$symbiont == 'EM', 'f.stat'] <- em.anova$`F value`[1]
div.fa[div.fa$symbiont == 'EM', 'p'] <- em.anova$`Pr(>F)`[1]

#--Plot
#--row 12 removed because it is an outlier
em.div.out <- em.div[-12,]

em.div <- ggplot(em.div.out, aes(x = group,
                                y = fishers.alpha)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Elevation group') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))

em.div

#--output plot
ggsave('Figure2a.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#--Mean and sd of Fisher's alpha per site
div.fa[div.fa$symbiont == 'EM', 'mean'] <- mean(em.div.site$fishers.alpha)
div.fa[div.fa$symbiont == 'EM', 'sd'] <- sd(em.div.site$fishers.alpha)

#=========================================================================================
# Endophytic fungi - Diversity
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data
#-----------------------------------------------------------------------------------------
#<< Read in pooled diversity >> ----------------------------------------------------------
#--diversity pooled by site and topography
fe.div <- read.csv(paste0(dat.dir, 'SCM_FE_div_pooled_by_landscape.csv'),
                   header = T, as.is = T)
#--diversity pooled only by site (for mean and sd)
fe.div.site <- read.csv (paste0(dat.dir, 'SCM_FE_div_pooled_by_site.csv'),
                         as.is = T, header = T)

#--exclude rows with n=s; outliers
fe.div <- fe.div[-which(fe.div$Elevation == 2421 & fe.div$topography == 'Convergent'),]
#-----------------------------------------------------------------------------------------
# ANOVA of FE diversity and boxplot
#-----------------------------------------------------------------------------------------
#--Difference between fisher's alpha means of the two goups
lm.fe.div <- lm(log (fishers.alpha) ~ group, data = fe.div)
fe.anova <- anova(lm.fe.div)
fe.anova
#--add data to results table
div.fa[div.fa$symbiont == 'FE', 'df.1'] <- fe.anova$Df[1]
div.fa[div.fa$symbiont == 'FE', 'df.2'] <- fe.anova$Df[2]
div.fa[div.fa$symbiont == 'FE', 'f.stat'] <- fe.anova$`F value`[1]
div.fa[div.fa$symbiont == 'FE', 'p'] <- fe.anova$`Pr(>F)`[1]

#--Plot
fe.div <- ggplot(fe.div, aes(x = group,
                                 y = fishers.alpha)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Elevation group') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))

fe.div

#--output plot
ggsave('Figure2b.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#--mean and sd of Fisher's alpha per site
div.fa[div.fa$symbiont == 'FE', 'mean'] <- mean(fe.div.site$fishers_alpha)
div.fa[div.fa$symbiont == 'FE', 'sd'] <- sd(fe.div.site$fishers_alpha)

#--write results table out
write.csv(div.fa, paste0(res.dir, 'Diversity_FishersAlpha.csv'), row.names = F)

#=========================================================================================
# Shannon's diversity
#=========================================================================================

#<< Read in OTU tables for calculating Shannon's diversity index >> ----------------------
#--EM site x species matrix
otu.data <- read.csv(paste0(dat.dir,'SCM_EM_otu_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)

#--FE site x species matrix
fe.data <- read.csv(paste0(dat.dir,'SCM_FE_site_x_species_matrix.csv'),
                    as.is = T, header = T)

#-----------------------------------------------------------------------------------------
# shannon's diversity for topographically similar trees within each site (n=20)
#-----------------------------------------------------------------------------------------

#--create data table
shannon <- data.frame(elevation = rep(unique(fe.data$elevation), each = 2),
                      topography = rep(c('Convergent', 'Divergent'), 10),
                      group = rep(c('High elev.','High elev.','Low elev.','Low elev.',
                                    'Low elev.','High elev.','High elev.','High elev.',
                                    'Low elev.', 'Low elev.'), each = 2),
                      em = NA, fe = NA)
#--calculate diversity
for (e in unique(shannon$elevation)) {
  for (t in c('Convergent', 'Divergent')) {
    shannon[shannon$elevation == e & shannon$topography == t, 'fe'] <-
      mean(diversity(fe.data[fe.data$elevation == e & fe.data$topography == t,
                                 2:42]))
    shannon[shannon$elevation == e & shannon$topography == t, 'em'] <-
      mean(diversity(otu.data[otu.data$elevation == e & otu.data$topography == t,
                                  5:length(otu.data)]))
  }
}

#--transform 
shannon$log.em <- log(shannon$em)
shannon$log.fe <- log(shannon$fe)

#-----------------------------------------------------------------------------------------
# shannon's diversity for each site (n=10)
#-----------------------------------------------------------------------------------------
#--create data table
shannon.site <- data.frame(elevation = unique(fe.data$elevation),
                           em = NA, fe = NA)
#--calculate diversity
for (e in shannon.site$elevation) {
  shannon.site[shannon.site$elevation == e, 'fe'] <-
    mean(diversity(fe.data[fe.data$elevation == e, 2:41]))
  shannon.site[shannon.site$elevation == e, 'em'] <-
    mean(diversity(otu.data[otu.data$elevation == e, 5:length(otu.data)]))
}

#--transform 
shannon.site$log.em <- log(shannon.site$em)
shannon.site$log.fe <- log(shannon.site$fe)

#-----------------------------------------------------------------------------------------
# ANOVA and boxplot of EM and FE diversity using Shannon's diversity
#-----------------------------------------------------------------------------------------
#<< ANOVA for EM >> ----------------------------------------------------------------------
shannon.em <- anova(lm(shannon$em ~ shannon$group))
shannon.em

#<< boxplot for EM >> --------------------------------------------------------------------
boxplot(shannon$em ~ shannon$group,
        ylab = "Shannon's diversity index", xlab = 'Elevation group',
        cex.lab = 1.5, cex.axis = 1)

#<< ANOVA for FE >> ----------------------------------------------------------------------
shannon.fe <- kruskal.test(shannon$fe ~ shannon$group)
shannon.fe

#<< boxplot for FE >> --------------------------------------------------------------------
boxplot(shannon$fe ~ shannon$group,
        ylab = "Shannon's diversity index", xlab = 'Elevation group',
        cex.lab = 1.5, cex.axis = 1)

#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#<< EM >> --------------------------------------------------------------------------------
#--mean and sd of Shannon's diversity
shannon.mean.em <- mean(shannon.site$em)
shannon.sd.em <- sd(shannon.site$em)

#<< FE >> --------------------------------------------------------------------------------
#--mean and sd of Shannon's diversity
shannon.sd.fe <- mean(shannon.site$fe)
shannon.sd.fe <- sd(shannon.site$fe)

write.csv()