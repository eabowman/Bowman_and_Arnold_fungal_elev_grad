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
#<< Read in site by species matrix >> ----------------------------------------------------
#--site
div.site.data <- read.csv(paste0(dat.dir,'SCM_EM_div_pooled_by_site_RootTip_based.csv'),
                          as.is = T, header = T) 
#--by topography
div.topo.data <- read.csv(paste0(dat.dir,'SCM_EM_div_pooled_by_landscape_RootTip_based.csv'),
                     as.is = T, header = T)

#--results table
div.out <- data.frame(test = c('FA','Shannon'), df.1 = NA, df.2 = NA, f.stat = NA, p = NA,
                     mean = NA, sd = NA)

#-----------------------------------------------------------------------------------------
# Linear regression of diversity as a function of elevation
#-----------------------------------------------------------------------------------------
#--Fisher's alpha
lm.em.div <- lm(log(fishers.alpha) ~ Elevation, data = div.topo.data)
em.anova <- anova(lm.em.div)
em.anova
#--add data to results table
div.out[div.out$test == 'FA', 'df.1'] <- em.anova$Df[1]
div.out[div.out$test == 'FA', 'df.2'] <- em.anova$Df[2]
div.out[div.out$test == 'FA', 'f.stat'] <- em.anova$`F value`[1]
div.out[div.out$test == 'FA', 'p'] <- em.anova$`Pr(>F)`[1]

#--Shannon's diversity index
lm.si.div <- lm(shannon.diversity ~ Elevation, data = div.topo.data)
em.anova <- anova(lm.si.div)
em.anova
#--add data to results table
div.out[div.out$test == 'Shannon', 'df.1'] <- em.anova$Df[1]
div.out[div.out$test == 'Shannon', 'df.2'] <- em.anova$Df[2]
div.out[div.out$test == 'Shannon', 'f.stat'] <- em.anova$`F value`[1]
div.out[div.out$test == 'Shannon', 'p'] <- em.anova$`Pr(>F)`[1]


#<< Plot >> ------------------------------------
em.div.plot <- ggplot(div.topo.data) +
  geom_point(aes(x = Elevation, y = fishers.alpha, color = plant.community),
             size = 5, shape = 0) +
  geom_point(aes(x = Elevation, y = shannon.diversity * 4.1 / 2.5, color = plant.community),
             size = 5, shape = 2) +
  theme_bw() +
  xlab('Elevation (m)') +
  ylab("Fisher's alpha") +
  scale_y_continuous(sec.axis = sec_axis(~ .* 4.1 / 6 , name = "Shannon's index"),
                     limits = c(0, 4.1)) +
  #scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Shannon's index")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28),
        legend.position = c(0.25,0.83),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  scale_color_manual(values = c('o.p' = '#d7191c', 'p' = '#fdae61',
                                'p.f' = '#67E98A', 'p.f.md' = '#2c7bb6'),
                     labels = c('Pine-oak', 'Pine', 'Pine-Doug fir',
                                'Pine-Doug fir-mixed decid.'),
                     name = 'Plant comm.') +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 5))) +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

em.div.plot

#--output plot
ggsave('Appendix_S5_fig1a.pdf', plot = em.div.plot, device = 'pdf', path = fig.dir)

#-----------------------------------------------------------------------------------------
# Site mean and standard deviation
#-----------------------------------------------------------------------------------------
#--Mean and sd of Fisher's alpha per site
div.out[div.out$test == 'FA', 'mean'] <- mean(div.topo.data$fishers.alpha)
div.out[div.out$test == 'FA', 'sd'] <- sd(div.topo.data$fishers.alpha)

#--mean and sd of Shannon's diversity
div.out[div.out$test == 'Shannon', 'mean'] <- mean(div.topo.data$shannon.diversity)
div.out[div.out$test == 'Shannon', 'sd']<- sd(div.topo.data$shannon.diversity)

write.csv(div.out,paste0(res.dir, 'diversity_EM_root_tip_based.csv'),row.names = F)
