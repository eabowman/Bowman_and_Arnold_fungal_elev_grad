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
em.div <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_landscape_OTU_based.csv'),
                   header = T, as.is = T)
#--diversity poooled only by site (for mean and sd)
em.div.site <- read.csv(paste0(dat.dir, 'SCM_EM_div_pooled_by_site_OTU_based.csv'),
                        as.is = T, header = T)
 #--results table
div.res <- data.frame(symbiont = c('EM.fa','EM.si','FE.fa','FE.si'),
                     df.1 = NA, df.2 = NA, f.stat = NA, p = NA, r2 = NA,
                     mean = NA, sd = NA)

#-----------------------------------------------------------------------------------------
# EM diversity analyses
#-----------------------------------------------------------------------------------------
#<< Fisher's alpha >> -----------
#--remove outlier
em.div.out <- em.div[-13,]
#--Linear regression of Fisher's alpha as a function of elevation (m)
lm.em.fa <- lm(log(fishers.alpha) ~ Elevation, data = em.div.out)
em.fa.anova <- summary(lm.em.fa)
em.fa.anova

#--add data to results table
div.res[div.res$symbiont == 'EM.fa', 'df.1'] <- em.fa.anova$fstatistic[2]
div.res[div.res$symbiont == 'EM.fa', 'df.2'] <- em.fa.anova$fstatistic[3]
div.res[div.res$symbiont == 'EM.fa', 'f.stat'] <- em.fa.anova$fstatistic[1]
div.res[div.res$symbiont == 'EM.fa', 'p'] <- em.fa.anova$coefficients[2,4]
div.res[div.res$symbiont == 'EM.fa', 'r2'] <- em.fa.anova$r.squared[1]

#<< Shannon's diversity index >> -----------
#--Linear regression of Fisher's alpha as a function of elevation (m)
lm.em.si <- lm(shannon.index ~ Elevation, data = em.div)
em.si.anova <- summary(lm.em.si)
em.si.anova

#--add data to results table
div.res[div.res$symbiont == 'EM.si', 'df.1'] <- em.si.anova$fstatistic[2]
div.res[div.res$symbiont == 'EM.si', 'df.2'] <- em.si.anova$fstatistic[3]
div.res[div.res$symbiont == 'EM.si', 'f.stat'] <- em.si.anova$fstatistic[1]
div.res[div.res$symbiont == 'EM.si', 'p'] <- em.si.anova$coefficients[2,4]
div.res[div.res$symbiont == 'EM.si', 'r2'] <- em.si.anova$r.squared[1]

#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#--Mean and sd of Fisher's alpha per site
div.res[div.res$symbiont == 'EM.fa', 'mean'] <- mean(em.div.site$fishers.alpha)
div.res[div.res$symbiont == 'EM.fa', 'sd'] <- sd(em.div.site$fishers.alpha)

#--mean and sd of Shannon's diversity
div.res[div.res$symbiont == 'EM.si', 'mean'] <- mean(em.div.site$shannon.index)
div.res[div.res$symbiont == 'EM.si', 'sd'] <- sd(em.div.site$shannon.index)

#-----------------------------------------------------------------------------------------
# Linear regression of EM species richness as a function elevation
#-----------------------------------------------------------------------------------------
lm.em.sr <- lm(em.div$species.richness ~ em.div$Elevation)
em.sr.anova <- summary(lm.em.sr)
em.sr.anova
plot(em.div$Elevation, em.div$species.richness)

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
# FE diversity analyses
#-----------------------------------------------------------------------------------------
#<< Fisher's alpha >> -----------
fit.fe.fa <- lm(log(fe.div$fishers.alpha) ~ poly(fe.div$Elevation, 2, raw=TRUE))
fe.fa.regression <- summary(fit.fe.fa)
fe.fa.regression

#--add data to results table
div.res[div.res$symbiont == 'FE.fa', 'df.1'] <- fe.fa.anova$fstatistic[2]
div.res[div.res$symbiont == 'FE.fa', 'df.2'] <- fe.fa.anova$fstatistic[3]
div.res[div.res$symbiont == 'FE.fa', 'f.stat'] <- fe.fa.anova$fstatistic[1]
div.res[div.res$symbiont == 'FE.fa', 'p'] <- fe.fa.anova$coefficients[2,4]
div.res[div.res$symbiont == 'FE.fa', 'r2'] <- fe.fa.anova$r.squared[1]

#<< Shannon's diversity index >> -----------
fit.fe.si <- lm(log(fe.div$shannon.index) ~ poly(fe.div$Elevation, 2, raw=TRUE))
fe.si.regression <- summary(fit.fe.si)
fe.si.regression

#--add data to results table
div.res[div.res$symbiont == 'FE.si', 'df.1'] <- fe.si.anova$fstatistic[2]
div.res[div.res$symbiont == 'FE.si', 'df.2'] <- fe.si.anova$fstatistic[3]
div.res[div.res$symbiont == 'FE.si', 'f.stat'] <- fe.si.anova$fstatistic[1]
div.res[div.res$symbiont == 'FE.si', 'p'] <- fe.si.anova$coefficients[2,4]
div.res[div.res$symbiont == 'FE.si', 'r2'] <- fe.si.anova$r.squared[1]
#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#--mean and sd of Fisher's alpha per site
div.res[div.res$symbiont == 'FE.fa', 'mean'] <- mean(fe.div.site$fishers_alpha)
div.res[div.res$symbiont == 'FE.fa', 'sd'] <- sd(fe.div.site$fishers_alpha)

#--mean and sd of Shannon's diversity index per site
div.res[div.res$symbiont == 'FE.si', 'mean'] <- mean(fe.div.site$shannon_index)
div.res[div.res$symbiont == 'FE.si', 'sd'] <- sd(fe.div.site$shannon_index)

#--write results table out
write.csv(div.res, paste0(res.dir, 'Diversity_FishersAlpha.csv'), row.names = F)

#-----------------------------------------------------------------------------------------
# Linear regression of FE species richness as a function elevation
#-----------------------------------------------------------------------------------------
lm.fe.sr <- lm(fe.div$species.richness ~ fe.div$Elevation)
fe.sr.anova <- anova(lm.fe.sr)
fe.sr.anova

#-----------------------------------------------------------------------------------------
# Plot with shannon and FA diversity
#-----------------------------------------------------------------------------------------
#<< Ectomycorhizal >> -------------------------------------------------------------------
#--Plot
#--row 12 removed because it is an outlier
em.div.out <- em.div[-12,]

em.div.plot <- ggplot(em.div.out) +
  geom_point(aes(x = Elevation, y = fishers.alpha, color = plant.community),
                 size = 5, shape = 0, stroke = 2) +
  geom_point(aes(x = Elevation, y = shannon.index * 40 / 2.5, color = plant.community),
                 size = 5, shape = 2, stroke = 2) +
  theme_bw() +
  xlab('Elevation (m)') +
  ylab("Fisher's alpha") +
  scale_y_continuous(sec.axis = sec_axis(~ .* 40 / 640 , name = "Shannon's index"),
                     limits = c(0, 40)) +
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Shannon's index")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        legend.position = c(0.25,0.89),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  scale_color_manual(values = c('o.p' = '#d7191c', 'p' = '#fdae61',
                                'p.f' = '#67E98A', 'p.f.md' = '#2c7bb6'),
                     labels = c('Pine-oak', 'Pine', 'Pine-Douglas fir',
                                'Pine-Douglas fir-mixed deciduous'),
                     name = 'Plant community') +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 5))) +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

em.div.plot

#--output plot
ggsave('Figure_1b.pdf', plot = em.div.plot, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#<< Endophyte >> -------------------------------------------------------------------
#--Plot
fe.div.out <- fe.div[-5,]

fe.div.plot <- ggplot(fe.div.out) +
  geom_point(aes(x = Elevation, y = fishers.alpha, color = plant.community),
             size = 5, shape = 0, stroke = 2) +
  geom_point(aes(x = Elevation, y = shannon.index * 10 / 2.1, color = plant.community),
             size = 5, shape = 2, stroke = 2) +
  theme_bw() +
  xlab('Elevation (m)') +
  ylab("Fisher's alpha") +
  scale_y_continuous(sec.axis = sec_axis(~ . * 10 / 43 , name = "Shannon's index"),
                     limits = c(0, 10)) +
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Shannon's index")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        legend.position = c(0.25,0.88),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  scale_color_manual(values = c('o.p' = '#d7191c', 'p' = '#fdae61',
                                'p.f' = '#67E98A', 'p.f.md' = '#2c7bb6'),
                     labels = c('Pine-oak', 'Pine', 'Pine-Douglas fir',
                                'Pine-Douglas fir-mixed deciduous'),
                     name = 'Plant community') +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 5))) +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

fe.div.plot

#--output plot
ggsave('Figure_2b.pdf', plot = fe.div.plot, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)
