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
# T-test of difference in abundance between high and low elevations and Plot
#-----------------------------------------------------------------------------------------
#<< T-test >> ----------------------------------------------------------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
em.stat[em.stat$abundance == 0, 'abundance'] <- 0.000001

#--T-test
em.ttest <- t.test(log(em.stat$abundance) ~ em.stat$group)
#--add results to table
abund.results[abund.results$symbiont == 'EM', 't-stat'] <- em.ttest$statistic[[1]]
abund.results[abund.results$symbiont == 'EM', 'df'] <- em.ttest$parameter[[1]]
abund.results[abund.results$symbiont == 'EM', 'p-value'] <- em.ttest$p.value[[1]]

#--Plot
em.abund <- ggplot(em.stat, aes(x = group,
                                y = abundance)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Elevation group') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))

em.abund

#--output plot
ggsave('Figure1a.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Linear regression and plot of EM abundance: Previous analysis
#-----------------------------------------------------------------------------------------
#<< Transformed abundance data linear regression per tree >> -----------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
#em.stat[em.stat$abundance == 0, 'abundance'] <- 0.000001
#--Linear regression of the two goups
#lm.em.ab <- lm(log(abundance) ~ elevation.m, data = em.stat)
#em.ab <- summary(lm.em.ab)
#--add results to results table
#abund.results[abund.results$symbiont == 'EM', 'df.1'] <- em.ab$fstatistic[2]
#abund.results[abund.results$symbiont == 'EM', 'df.2'] <- em.ab$fstatistic[3]
#abund.results[abund.results$symbiont == 'EM', 'R.sq'] <- em.ab$r.squared[1]
#abund.results[abund.results$symbiont == 'EM', 'f.stat'] <- em.ab$fstatistic[1]
#abund.results[abund.results$symbiont == 'EM', 'p'] <- em.ab$coefficients[2,4]

#<< Plot of mean abundance and sd by elevation >> ----------------------------------------
#--make dataframe with mean and sd per elevation
#em.mean <- data.frame(elevation = unique(em.stat$elevation.m),
                       mean.em = NA, sd.em = NA)
#--add mean
#elev <- em.mean$elevation
#for(i in elev) {
#  em.mean[em.mean$elevation == i, 'mean.em'] <- 
#    mean(em.stat[em.stat$elevation.m == i, 'abundance'], na.rm = T)
#  em.mean[em.mean$elevation == i, 'sd.em'] <-
#    sd(em.stat[em.stat$elevation.m == i, 'abundance'], na.rm = T)
#}
#--plot
#em.ab <- qplot(em.mean$elevation, em.mean$mean.em) + 
#  geom_errorbar(aes(x= em.mean$elevation,
#      ymin = em.mean$mean.em-em.mean$sd.em, ymax = em.mean$mean.em + em.mean$sd.em),
#      width=0.25) +
#  theme_classic(base_size = 16, base_family = '') +
#  #geom_smooth(method = 'lm', se = F, col = 'black', weight = 0.2) +
#  ylab('Abundance') +
#  xlab('Elevation (m)') +
#  theme(panel.background = element_rect(color = "black"), text = element_text(size=20))

#--output plot
#ggsave('Figure1a.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

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

#-----------------------------------------------------------------------------------------
# T-test of difference in abundance between high and low elevations and Plot
#-----------------------------------------------------------------------------------------
#<< T-test >> ----------------------------------------------------------------------------
#--Change those sites with zero abundance to NA, so that log abundance can be calc.
fe.abun[fe.abun$abundance == 0, 'abundance'] <- 0.000001

#--T-test
fe.ttest <- t.test(logit(fe.abun$abundance) ~ fe.abun$group)
#--add results to table
abund.results[abund.results$symbiont == 'FE', 't-stat'] <- fe.ttest$statistic[[1]]
abund.results[abund.results$symbiont == 'FE', 'df'] <- fe.ttest$parameter[[1]]
abund.results[abund.results$symbiont == 'FE', 'p-value'] <- fe.ttest$p.value[[1]]

#--Plot
fe.abund <- ggplot(fe.abun, aes(x = group,
                                y = abundance)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Elevation group') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))

fe.abund

#--output plot
ggsave('Figure1b.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Linear regression and plot of FE abundance
#-----------------------------------------------------------------------------------------
#<< Transformed per tree >> --------------------------------------------------------------
#--Linear regression using logit transformed abundance data
#lm.fe.ab <- lm(logit(abundance) ~ elevation.m, data = fe.abun)
#fe.ab <- summary(lm.fe.ab)
#--add results to results table
#abund.results[abund.results$symbiont == 'FE', 'df.1'] <- fe.ab$fstatistic[2]
#abund.results[abund.results$symbiont == 'FE', 'df.2'] <- fe.ab$fstatistic[3]
#abund.results[abund.results$symbiont == 'FE', 'R.sq'] <- fe.ab$r.squared[1]
#abund.results[abund.results$symbiont == 'FE', 'f.stat'] <- fe.ab$fstatistic[1]
#abund.results[abund.results$symbiont == 'FE', 'p'] <- fe.ab$coefficients[2,4]

#<< Plot of mean abundance and sd by elevation >> ----------------------------------------
#--make dataframe with mean and sd per elevation
#fe.mean <- data.frame (elevation = unique(fe.abun$elevation.m),
                       mean.fe = NA, sd.fe = NA)
#--add mean
#elev <- fe.mean$elevation
#for (i in elev) {
#  fe.mean[fe.mean$elevation == i, 'mean.fe'] <- 
#    mean(fe.abun[fe.abun$elevation == i, 'abundance'], na.rm = T)
#  fe.mean[fe.mean$elevation == i, 'sd.fe'] <-
#    sd(fe.abun[fe.abun$elevation == i, 'abundance'], na.rm = T)
#}

#--plot
#fe.ab <- qplot(fe.mean$elevation, fe.mean$mean.fe) + 
#  geom_errorbar(aes(x= fe.mean$elevation,
#          ymin = fe.mean$mean.fe-fe.mean$sd.fe, ymax = fe.mean$mean.fe + fe.mean$sd.fe),
#          width=0.25) +
#  theme_classic(base_size = 16, base_family = '') +
#  geom_smooth(method = 'lm', se = F, col = 'black', weight = 0.2) +
#  ylab('Abundance') +
#  xlab('Elevation (m)') +
#  theme(panel.background = element_rect(color = "black"), text = element_text(size=20))

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


#=========================================================================================
# Ectomycorrhizal fungi - Plant community
#=========================================================================================
#--differences in abundance based on plant community
anova(lm(log(abundance) ~ p.c, data = em.stat))

em.pc <- ggplot(em.stat, aes(x = p.c,
                    y = abundance)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Plant community') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))
em.pc

#=========================================================================================
# Endophyte fungi - Plant community
#=========================================================================================
#--differences in abundance based on plant community
anova(lm(logit(abundance) ~ p.c, data = fe.abun))

for(i in unique(em.stat$elevation.m)){
  fe.abun[fe.abun$elevation.m == i, 'p.c'] <-
    unique(em.stat[em.stat$elevation.m == i, 'p.c'])
}

fe.pc <- ggplot(fe.abun, aes(x = p.c,
                             y = abundance)) +
  geom_boxplot() +
  theme_bw() +
  xlab('Plant community') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12), axis.title = element_text(size = 14))
fe.pc
