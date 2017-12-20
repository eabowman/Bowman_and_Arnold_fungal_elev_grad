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
otu.data <- read.csv(paste0(dat.dir, 'SCM_EM_root_based_site_x_species_matrix.csv'),
                     as.is = T, header = T)

#-----------------------------------------------------------------------------------------
# EM species accumulation curve based on root tip abundance
#-----------------------------------------------------------------------------------------
#<< With singeltons included, Overall EM >>-----------------------------------------------
#--Isolate OTU data
ov.em <- otu.data[6:length(otu.data)]
#--remove columns equal to 0, no species occurrence 
ov.em <- ov.em[which(colSums(ov.em) != 0)]
#--abundance for each otu, double check for discrepancies
summarise_each(ov.em,funs(sum))

#<< With singeltons removed, Overall EM >>------------------------------------------------
#--Using previous dataframe, remove singletons
ovns.em <- ov.em[which(colSums(ov.em) > 1)]
#--abundance for each otu, double check for discrepancies
colSums(ovns.em)

#-----------------------------------------------------------------------------------------
# EM species accumulation curve based on root tip count, combined plots
#-----------------------------------------------------------------------------------------

#<< Combine overall EM plots, with singletons and without >>------------------------------
em.all <- specaccum(ov.em, sample = min(rowSums(ov.em),  permutations = 999))
em.all.df <- data.frame(Sites=em.all$sites,
                         Richness=em.all$richness,
                         SD=em.all$sd)
em.all.ns <- specaccum(ovns.em, sample = min(rowSums(ovns.em),  permutations = 999))
em.all.ns.df <- data.frame(Sites=em.all.ns$sites,
                            Richness=em.all.ns$richness,
                            SD=em.all.ns$sd)
ggplot() +
  geom_point(data=em.all.df, aes(x=Sites, y=Richness)) +
  geom_line(data=em.all.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=em.all.df ,aes(x=Sites,
                                   ymin=(Richness-2*SD),
                                   ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=em.all.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_line(data=em.all.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=em.all.ns.df ,aes(x=Sites,
                                      ymin=(Richness-2*SD),
                                      ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,150)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=28, color = 'black'),
        axis.title = element_text(size = 28)) +
  theme(legend.position = "none")

ggsave('Appendix_S5_Fig1b.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)

#<< Combine low and highelev. EM plots, with singletons and without >>--------------------
# em.low <- specaccum(le.em, sample = min(rowSums(le.em),  permutations = 999))
# em.low.df <- data.frame(Sites=em.low$sites,
#                         Richness=em.low$richness,
#                         SD=em.low$sd)
# em.low.ns <- specaccum(lens.em, sample = min(rowSums(lens.em),  permutations = 999))
# em.low.ns.df <- data.frame(Sites=em.low.ns$sites,
#                            Richness=em.low.ns$richness,
#                            SD=em.low.ns$sd)
# em.high <- specaccum(he.em, sample = min(rowSums(he.em), permutations = 999))
# em.high.df <- data.frame(Sites=em.high$sites,
#                          Richness=em.high$richness,
#                          SD=em.high$sd)
# em.high.ns <- specaccum(hens.em, sample = min(rowSums (hens.em),  permutations = 999))
# em.high.ns.df <- data.frame(Sites=em.high.ns$sites,
#                             Richness=em.high.ns$richness,
#                             SD=em.high.ns$sd)
# ggplot() +
#   geom_point(data=em.low.df, aes(x=Sites, y=Richness), size = 3) +
#   geom_line(data=em.low.df, aes(x=Sites, y=Richness)) +
#   geom_ribbon(data=em.low.df ,aes(x=Sites,
#                                   ymin=(Richness-2*SD),
#                                   ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=em.low.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
#   geom_line(data=em.low.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
#   geom_ribbon(data=em.low.ns.df ,aes(x=Sites,
#                                      ymin=(Richness-2*SD),
#                                      ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=em.high.df, aes(x=Sites, y=Richness), shape = 5, size = 3) +
#   geom_line(data=em.high.df, aes(x=Sites, y=Richness)) +
#   geom_ribbon(data=em.high.df ,aes(x=Sites,
#                                    ymin=(Richness-2*SD),
#                                    ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=em.high.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey',
#              shape = 5, size = 3) +
#   geom_line(data=em.high.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
#   geom_ribbon(data=em.high.ns.df ,aes(x=Sites,
#                                       ymin=(Richness-2*SD),
#                                       ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   theme_bw() +
#   expand_limits(y=c(0,150)) +
#   ylab('OTUs') +
#   xlab('Samples') +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(size=28, color = 'black'),
#         axis.title = element_text(size = 28)) +
#   theme(legend.position = "none")
# 
# ggsave('Appendix_S5_fig1c.pdf', plot = last_plot(), device = 'pdf', path = fig.dir)
