## R code written by Elizabeth A. Bowman Oct. 27, 2016
## University of Arizona, School of Plant Sciences, eabowman@email.arizona.edu
## Analyses to evaluate the taxonomic differences in the two fungal communities and 
## elevation groups (i.e., high and low).

#=========================================================================================
# Taxonomic composition: Chi-square test
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Read in data frames
#-----------------------------------------------------------------------------------------

#--Read in EM data frame
em.data <- read.csv(paste0(dat.dir, 'SCM_EM_raw.csv'), as.is = T, header = T)

#--Read in Environmental data for plant community information
env.data <- read.csv(paste0(dat.dir,'SCM_environmental_data.csv'))

#--Make results table
org <- c("em")
tax.chisq <- data.frame(org)

env.data$plant.comm <- as.character(env.data$plant.comm)

for(e in unique(em.data$elevation)){
  em.data[em.data$elevation == e, 'plant.comm'] <-
    unique(env.data[env.data$elevation == e, 'plant.comm'])
}

#=========================================================================================
# Ectomycorrhizal fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean-up data
#-----------------------------------------------------------------------------------------

#<< Make count table of em classes by morphotype data >> ---------------------------------
group = NA
taxonomy_class = NA
for(r in 1:nrow(em.data)){
  tax.n <- rep(em.data[r, 'taxonomy_class'], em.data[r, 'abundance'])
  group.n <- rep(em.data[r, 'plant.comm'], em.data[r, 'abundance'])
  group <- c(group, group.n)
  taxonomy_class <- c(taxonomy_class, tax.n)
}
group <- group[-1]
taxonomy_class <- taxonomy_class[-1]
em.morph <- data.frame(plant.comm = group,
                       taxonomy_class = taxonomy_class)

elvgr.em.morph <- table(em.morph$plant.comm,
                        em.morph$taxonomy_class)

#-----------------------------------------------------------------------------------------
# Chi-square test
#-----------------------------------------------------------------------------------------
#--Chi square test: 
em.chi <- chisq.test(elvgr.em.morph)
em.chi

#--Add results to results table; chi-square, degrees of freedom, and the p-value
tax.chisq[which(tax.chisq$org == "em"), "chi"] <- em.chi$statistic[[1]]
tax.chisq[which(tax.chisq$org == "em"), "df"] <- em.chi$parameter[[1]]
tax.chisq[which(tax.chisq$org == "em"), "p.value"] <- em.chi$p.value[[1]]

#-----------------------------------------------------------------------------------------
# 100% bar chart of taxa
#-----------------------------------------------------------------------------------------
#--Remove taxa with rare taxa
rare.em <- c("Fungi incertae sedis", "Eurotiomycetes", "Sordariomycetes")
em.chisq <- em.data[which(!em.data$taxonomy_class %in% rare.em), ]
em.chisq$elevation <- as.factor(em.chisq$elevation)
#--Bar graph
tax.em <- ggplot(data = em.chisq, 
          aes(x = elevation,
          fill = taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  xlab("Elevation (m)") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Blues") +
  guides (fill=guide_legend(title=NULL)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=28, colour = 'black'),
        axis.title = element_text(size = 28), legend.position = "bottom",
        legend.text=element_text(size=14),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x = element_text(angle = 45, hjust = 1))
tax.em

#--save plot to figures folder
ggsave('Appendix_S5_fig1d.pdf', plot = tax.em, device = 'pdf', path = fig.dir)

