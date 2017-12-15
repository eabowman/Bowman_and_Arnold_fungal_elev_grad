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
ov.em <- otu.data[6:length(otu.data)]
#--remove columns equal to 0, no species occurrence 
ov.em <- ov.em[which(colSums(ov.em) != 0)]
#--abundance for each otu, double check for discrepancies
summarise_each(ov.em,funs(sum))

#<< With singeltons included, Low elevation group, EM >>----------------------------------
#--Isolate OTU data, low elevation
le.em <- otu.data[which(otu.data$group == "Low elev."), ]
le.em <- le.em[6:length(le.em)]
#--remove columns equal to 0, no species occurence
le.em <- le.em[which(colSums(le.em) != 0)]

#<< With singeltons included, High elevation group, EM >>---------------------------------
#--Isolate OTU data, high elevation
he.em <- otu.data[which (otu.data$group == "High elev."), ]
he.em <- he.em[6:length(le.em)]
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

#<< Combine overall EM plots, with singletons and without >>------------------------------
em.all <- specaccum(ov.em, sample = min(rowSums(ov.em),  permutations = 999))
em.all.df <- data.frame(Sites=em.all$sites,
                        Richness=em.all$richness,
                        SD=em.all$sd)
em.all.ns <- specaccum(ovns.em, sample = min(rowSums(ovns.em),  permutations = 999))
em.all.ns.df <- data.frame(Sites=em.all.ns$sites,
                           Richness=em.all.ns$richness,
                           SD=em.all.ns$sd)
EM.1c <- ggplot() +
  geom_point(data=em.all.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=em.all.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=em.all.df ,aes(x=Sites,
                                  ymin=(Richness-2*SD),
                                  ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=em.all.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
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
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('Figure_1c.pdf', plot = EM.1c, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#<< Combine low elev. EM plots, with singletons and without >>----------------------------
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
#                         Richness=em.high$richness,
#                         SD=em.high$sd)
# em.high.ns <- specaccum(hens.em, sample = min(rowSums (hens.em),  permutations = 999))
# em.high.ns.df <- data.frame(Sites=em.high.ns$sites,
#                            Richness=em.high.ns$richness,
#                            SD=em.high.ns$sd)
# EM.1d <- ggplot() +
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
#   geom_point(data=em.high.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', shape = 5, size = 3) +
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
#         axis.title = element_text(size = 28), legend.position = "none") +
#   theme(axis.title.x = element_text(margin = margin(t = 30)),
#         axis.title.y = element_text(margin = margin(r = 30)))
# 
# ggsave('Figure_1d.pdf', plot = EM.1d, device = 'pdf', path = fig.dir)

#<< Combine overall FE plots, with singletons and without >>------------------------------
fe.all <- specaccum(ov.fe, sample = min(rowSums (ov.fe),  permutations = 999))
fe.all.df <- data.frame(Sites=fe.all$sites,
                         Richness=fe.all$richness,
                         SD=fe.all$sd)
fe.all.ns <- specaccum(ovns.fe, sample = min(rowSums (ovns.fe),  permutations = 999))
fe.all.ns.df <- data.frame(Sites=fe.all.ns$sites,
                            Richness=fe.all.ns$richness,
                            SD=fe.all.ns$sd)
fe.3c <- ggplot() +
  geom_point(data=fe.all.df, aes(x=Sites, y=Richness), size = 3) +
  geom_line(data=fe.all.df, aes(x=Sites, y=Richness)) +
  geom_ribbon(data=fe.all.df ,aes(x=Sites,
                                   ymin=(Richness-2*SD),
                                   ymax=(Richness+2*SD)),
              alpha=0.2) +
  geom_point(data=fe.all.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
  geom_line(data=fe.all.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
  geom_ribbon(data=fe.all.ns.df ,aes(x=Sites,
                                      ymin=(Richness-2*SD),
                                      ymax=(Richness+2*SD)),
              alpha=0.2) +
  theme_bw() +
  expand_limits(y=c(0,40)) +
  ylab('OTUs') +
  xlab('Trees sampled') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=22, color = 'black'),
        axis.title = element_text(size = 28), legend.position = "none") +
  theme(axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)))

ggsave('Figure_2c.pdf', plot = fe.3c, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#<< Combine low elev. FE plots, with singletons and without >>----------------------------
# fe.low <- specaccum(le.fe, sample = min(rowSums (le.fe),  permutations = 999))
# fe.low.df <- data.frame(Sites=fe.low$sites,
#                         Richness=fe.low$richness,
#                         SD=fe.low$sd)
# fe.low.ns <- specaccum(lens.fe, sample = min(rowSums(lens.fe),  permutations = 999))
# fe.low.ns.df <- data.frame(Sites=fe.low.ns$sites,
#                            Richness=fe.low.ns$richness,
#                            SD=fe.low.ns$sd)
# fe.high <- specaccum(he.fe, sample = min(rowSums(he.fe),  permutations = 999))
# fe.high.df <- data.frame(Sites=fe.high$sites,
#                         Richness=fe.high$richness,
#                         SD=fe.high$sd)
# fe.high.ns <- specaccum(hens.fe, sample = min(rowSums (hens.fe),  permutations = 999))
# fe.high.ns.df <- data.frame(Sites=fe.high.ns$sites,
#                            Richness=fe.high.ns$richness,
#                            SD=fe.high.ns$sd)
# 
# fe.3d <- ggplot() +
#   geom_point(data=fe.low.df, aes(x=Sites, y=Richness), size = 3) +
#   geom_line(data=fe.low.df, aes(x=Sites, y=Richness)) +
#   geom_ribbon(data=fe.low.df ,aes(x=Sites,
#                                   ymin=(Richness-2*SD),
#                                   ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=fe.low.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', size = 3) +
#   geom_line(data=fe.low.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
#   geom_ribbon(data=fe.low.ns.df ,aes(x=Sites,
#                                      ymin=(Richness-2*SD),
#                                      ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=fe.high.df, aes(x=Sites, y=Richness), shape = 5, size = 3) +
#   geom_line(data=fe.high.df, aes(x=Sites, y=Richness)) +
#   geom_ribbon(data=fe.high.df ,aes(x=Sites,
#                                    ymin=(Richness-2*SD),
#                                    ymax=(Richness+2*SD)),
#               alpha=0.2) +
#   geom_point(data=fe.high.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey', shape = 5, size = 3) +
#   geom_line(data=fe.high.ns.df, aes(x=Sites, y=Richness), colour = 'darkgrey') +
#   geom_ribbon(data=fe.high.ns.df ,aes(x=Sites,
#                                       ymin=(Richness-2*SD),
#                                       ymax=(Richness+2*SD)),
#               alpha = 0.2) +
#   theme_bw() +
#   expand_limits(y=c(0,40)) +
#   ylab('OTUs') +
#   xlab('Samples') +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(size=28, color = 'black'),
#         axis.title = element_text(size = 28), legend.position = "none") +
#   theme(axis.title.x = element_text(margin = margin(t = 30)),
#         axis.title.y = element_text(margin = margin(r = 30)))
# 
# ggsave('Figure_3d.pdf', plot = fe.3d, device = 'pdf', path = fig.dir)

