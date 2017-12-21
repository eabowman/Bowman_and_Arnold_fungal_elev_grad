## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the abundance of ectomycorrhizal and endophytic abundance as a
## function of elevation as part of a study of these two fungal communities along an
## elevation gradient in the Santa Catalina Mountains

#=========================================================================================
# Ectomycorrhizal fungi - Abundance 
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------
#--read in EM data frame
em.stat <- read.csv(paste0(dat.dir,'SCM_EM_per_tree_stats.csv'),
                    as.is = T, header = T)
#--change group category names
em.stat[em.stat$group == 'High elev.', 'group'] <- 'High'
em.stat[em.stat$group == 'Low elev.', 'group'] <- 'Low'

#--read in FE data frame
fe.abun <- read.csv(paste0(dat.dir,'SCM_FE_per_tree_stats.csv'),
                    as.is = T, header = T)
#--add elevation groups to fe.abun dataframe
fe.abun[fe.abun$elevation.m > 2300, 'group'] <- 'High'
fe.abun[fe.abun$elevation.m < 2300, 'group'] <- 'Low'


#--make results table
abund.results <- data.frame(symbiont = c('EM', 'FE'), per.tree.mean = NA, per.tree.sd = NA)

#-----------------------------------------------------------------------------------------
# Linear regression and plot of EM abundance
#-----------------------------------------------------------------------------------------
#<< Transformed abundance data linear regression per tree >> -----------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
em.stat[em.stat$abundance == 0, 'abundance'] <- 0.000001
#--Linear regression of the two goups
lm.em.ab <- lm(log(abundance) ~ elevation.m, data = em.stat)
em.ab <- summary(lm.em.ab)

#--add results to results table
abund.results[abund.results$symbiont == 'EM', 'df.1'] <- em.ab$fstatistic[2]
abund.results[abund.results$symbiont == 'EM', 'df.2'] <- em.ab$fstatistic[3]
abund.results[abund.results$symbiont == 'EM', 'R.sq'] <- em.ab$r.squared[1]
abund.results[abund.results$symbiont == 'EM', 'f.stat'] <- em.ab$fstatistic[1]
abund.results[abund.results$symbiont == 'EM', 'p'] <- em.ab$coefficients[2,4]

#<< Plot of mean abundance and sd by elevation >> ----------------------------------------
#--make dataframe with mean and sd per elevation
em.mean <- data.frame(elevation = rep(unique(em.stat$elevation.m),2),
                      topography = c(rep('Convergent', 10),rep('Divergent',10)),
                      mean.em = NA, p.c = NA)
#--add mean
elev <- unique(em.mean$elevation)
topo <- c('Convergent', 'Divergent')
for(i in elev){
  for(t in topo) {
    em.mean[em.mean$elevation == i & em.mean$topography == t, 'mean.em'] <-
      mean(em.stat[em.stat$elevation.m == i & em.stat$topography == t,
                   'abundance'], na.rm = T)
    em.mean[em.mean$elevation == i, 'p.c'] <- 
      unique(em.stat[em.stat$elevation.m == i, 'p.c'])
  }
}

#--plot
em.ab <- ggplot(em.mean, aes(x = elevation,
                             y = mean.em,
                             color = p.c)) +
  geom_point(size = 5, shape = 0, stroke = 2) +
  theme_bw() +
  xlab("Elevation (m)") +
  ylab("Abundance (tips·g-1)") +
  #geom_errorbar(aes(ymin = mean.em - sd.em, ymax = mean.em + sd.em)) +
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
em.ab

#--output plot
ggsave('Figure_1a.pdf', plot = em.ab, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#-----------------------------------------------------------------------------------------
# T-test of difference in abundance between high and low elevations and Plot
# Previous analysis
#-----------------------------------------------------------------------------------------
#<< T-test >> ----------------------------------------------------------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
# em.stat[em.stat$abundance == 0, 'abundance'] <- 0.000001

#--T-test
# em.ttest <- t.test(log(em.stat$abundance) ~ em.stat$group)
# #--add results to table
# abund.results[abund.results$symbiont == 'EM', 't-stat'] <- em.ttest$statistic[[1]]
# abund.results[abund.results$symbiont == 'EM', 'df'] <- em.ttest$parameter[[1]]
# abund.results[abund.results$symbiont == 'EM', 'p-value'] <- em.ttest$p.value[[1]]

# #--Plot
# em.abund <- ggplot(em.stat, aes(x = elevation.m,
#                                    y = abundance,
#                                    color = p.c)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   xlab('Elevation (m)') +
#   ylab("Abundance (tips·g-1)") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(size=12), axis.title = element_text(size = 14)) +
#   theme(legend.position = "none")
# 
# em.abund
# 
# #--output plot
# ggsave('Figure1a.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Per tree mean and standard deviation
#-----------------------------------------------------------------------------------------
#--mean EM abundance per tree
abund.results[abund.results$symbiont == 'EM', 'per.tree.mean'] <- mean(em.stat$abundance)
#--sd of EM abundance per tree
abund.results[abund.results$symbiont == 'EM', 'per.tree.sd'] <- sd(em.stat$abundance)

#=========================================================================================
# Endophyte fungi - Abundance 
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frame
#-----------------------------------------------------------------------------------------
#--remove this where needles were not collected
fe.abun <- fe.abun[!is.na(fe.abun$endophyte_collection_date), ]

#--Change those sites with zero abundance to NA, so that log abundance can be calc.
fe.abun[fe.abun$abundance == 0, 'abundance'] <- 0.000001
#-----------------------------------------------------------------------------------------
# Linear regression and plot of FE abundance
#-----------------------------------------------------------------------------------------
#<< Transformed per tree >> --------------------------------------------------------------
#--Linear regression using logit transformed abundance data
lm.fe.ab <- lm(logit(abundance) ~ elevation.m, data = fe.abun)
fe.ab <- summary(lm.fe.ab)

#--add results to results table
abund.results[abund.results$symbiont == 'FE', 'df.1'] <- fe.ab$fstatistic[2]
abund.results[abund.results$symbiont == 'FE', 'df.2'] <- fe.ab$fstatistic[3]
abund.results[abund.results$symbiont == 'FE', 'R.sq'] <- fe.ab$r.squared[1]
abund.results[abund.results$symbiont == 'FE', 'f.stat'] <- fe.ab$fstatistic[1]
abund.results[abund.results$symbiont == 'FE', 'p'] <- fe.ab$coefficients[2,4]

#<< Plot of mean abundance and sd by elevation >> ----------------------------------------
#--make dataframe with mean and sd per elevation
fe.mean <- data.frame (elevation = unique(fe.abun$elevation.m),
                       topography = c(rep('Convergent', 10),rep('Divergent',10)),
                       mean.fe = NA, p.c = NA)

#--add mean
elev <- unique(fe.mean$elevation)
topo <- c('Convergent','Divergent')
for(i in elev){
  for(t in topo){
    fe.mean[fe.mean$elevation == i & fe.mean$topography == t, 'mean.fe'] <- 
      mean(fe.abun[fe.abun$elevation.m == i & fe.abun$topography == t,
                   'abundance'], na.rm = T)
  fe.mean[fe.mean$elevation == i, 'p.c'] <- em.mean[em.mean$elevation == i, 'p.c']
  }
}

#--plot
fe.ab <- ggplot(fe.mean, aes(x = elevation,
                             y = mean.fe,
                             color = p.c)) +
  geom_point(size = 5, shape = 0, stroke = 2) +
  theme_bw() +
  # geom_smooth(aes(colour=NA),
  #             method=lm,   # Add linear regression line
  #             se=FALSE) +
  ylab("Isolation frequency") +
  xlab('Elevation (m)') +
  #geom_errorbar(aes(ymin = mean.fe - sd.fe, ymax = mean.fe + sd.fe)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        legend.position = c(0.78,0.88),
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
fe.ab

#--output plot
ggsave('Figure_2a.pdf', plot = fe.ab, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#-----------------------------------------------------------------------------------------
# T-test of difference in abundance between high and low elevations and Plot:
# Previous analysis
#-----------------------------------------------------------------------------------------
#<< T-test >> ----------------------------------------------------------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
# fe.abun[fe.abun$abundance == 0, 'abundance'] <- 0.000001

#--T-test
# fe.ttest <- t.test(logit(fe.abun$abundance) ~ fe.abun$group)
# #--add results to table
# abund.results[abund.results$symbiont == 'FE', 't-stat'] <- fe.ttest$statistic[[1]]
# abund.results[abund.results$symbiont == 'FE', 'df'] <- fe.ttest$parameter[[1]]
# abund.results[abund.results$symbiont == 'FE', 'p-value'] <- fe.ttest$p.value[[1]]

#--Plot
# fe.abund <- ggplot(fe.abun, aes(x = group,
#                                 y = abundance)) +
#   geom_boxplot() +
#   theme_bw() +
#   xlab('Elevation group') +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))
# 
# fe.abund

#--output plot
#ggsave('Figure1b.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Per tree mean and standard deviation
#-----------------------------------------------------------------------------------------
#--mean EM abundance per tree
abund.results[abund.results$symbiont == 'FE', 'per.tree.mean'] <-
  mean(fe.abun$abundance)*100
#--sd of EM abundance per tree
abund.results[abund.results$symbiont == 'FE', 'per.tree.sd'] <-
  sd(fe.abun$abundance)*100

#--write out results to results folder
write.csv(abund.results, paste0(res.dir, 'Bowman_and_Arnold_2016_Abundance.csv'),
          row.names = F)
