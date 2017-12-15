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
#--Read in FE data frame
fe.data <- read.csv(paste0(dat.dir, 'SCM_fe_raw.csv'), as.is = T, header = T)

#--Read in environmental data
env.data <- read.csv(paste0(dat.dir,'SCM_environmental_data.csv'))

env.data$plant.comm <- as.character(env.data$plant.comm)

#--Add plant community data to em.data and fe.data
for(e in unique(em.data$elevation)) {
  for(x in em.data$elevation){
    if(x == e){
      em.data[em.data$elevation == e, 'plant.comm'] <-
        unique(env.data[env.data$elevation == e, 'plant.comm'])
    }
  }
}

tree_no <- c("LB004","LB005","LB006","LB007","LB008","LB009","LB010","LB011","LB012",
             "LB012 adjacent","LB013","LB014","LB015","LB016","LB017","LB018","LB019",
             "LB020","LB021","LB023","LB024","LB025","LB026","LB027","LB029","LB030",
             "LB031","LB032","LB033","LB034","LB035","LB037","LB038","LB039","LB040",
             "LB041","LB042","LB043","LB044","LB045","LB046","LB047","LB048","LB049",
             "LB050","LB051","LB052","LB053","LB054","LB055","LB056","LB057","LB058",
             "LB059","LB060","LB062","LB063","LB064","LB065","LB066")
for(t in tree_no) {
  for(x in env.data$tree_number){
    if(x == t){
      fe.data[fe.data$Tree_number == t, 'plant.comm'] <-
        unique(env.data[env.data$tree_number == t, 'plant.comm'])
    }
  }
}

fe.op <- c('LB022','LB061')
fe.pf <- 'LB036'
fe.data[fe.data$Tree_number %in% fe.op, 'plant.comm'] <- 'o.p'
fe.data[fe.data$Tree_number %in% fe.pf, 'plant.comm'] <- 'p.f'


#--Make results table
org <- c("em", "fe")
tax.chisq <- data.frame(org)
#=========================================================================================
# Ectomycorrhizal fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean-up data
#-----------------------------------------------------------------------------------------
#--Make count table of em classes
pc.em <- table(em.data$plant.comm,
                  em.data$taxonomy_class)
#--Remove taxa with equal to or less than 1 occurence
rare.em <- c("Fungi incertae sedis", 'Eurotiomycetes')
pc.em <- pc.em[ ,which(!colnames(pc.em) %in% rare.em)]

#-----------------------------------------------------------------------------------------
# Chi-square test
#-----------------------------------------------------------------------------------------
#--Chi square test
em.chi <- chisq.test(pc.em)

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
em.chisq$elevation <-  as.factor(em.chisq$elevation)
#--Bar graph
tax.em <- ggplot(data = em.chisq, 
          aes(x = elevation,
          fill = taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  xlab("Elevation (m)") +
  #ggtitle("Proportion of Classes by Topography") +
  #scale_fill_brewer(palette = "Blues") +
  scale_fill_manual(values = c('#c6dbef', '#9ecae1', '#6baed6', '#4292c6')) +
  guides(fill=guide_legend(title=NULL)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=24, colour = 'black'),
        axis.title = element_text(size = 24), legend.position = 'bottom',
        legend.text=element_text(size=14),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x = element_text(angle = 45, hjust = 1))

tax.em

#--save plot to figures folder
ggsave('Figure_1f.pdf', plot = tax.em, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)

#=========================================================================================
# Endophytic fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean-up data
#-----------------------------------------------------------------------------------------
#--Make count table of em classes
pc.fe <- table(fe.data$plant.comm,
                  fe.data$Taxonomy)
#--For removing those classes with only 1 occurence 
rare.fe <- c("Lecanoromycetes", "Orbiliomycetes", "Tremellomycetes","Pezizomycetes",
             "Unknown")
pc.fe <- pc.fe[,!colnames(pc.fe) %in% rare.fe]

#-----------------------------------------------------------------------------------------
# Chi-square test
#-----------------------------------------------------------------------------------------
#--Chi square test
fe.chi <- chisq.test(pc.fe)
fe.chi

#--Add results to results table; chi-square, degrees of freedom, and the p-value
tax.chisq[which(tax.chisq$org == "fe"), "chi"] <- fe.chi$statistic[[1]]
tax.chisq[which(tax.chisq$org == "fe"), "df"] <- fe.chi$parameter[[1]]
tax.chisq[which(tax.chisq$org == "fe"), "p.value"] <- fe.chi$p.value[[1]]

#--write out results table
write.csv(tax.chisq, paste0(res.dir, 'Taxonomic_analyses.csv'),
          row.names = F)

#-----------------------------------------------------------------------------------------
# 100% bar chart of taxa
#-----------------------------------------------------------------------------------------
#--remove unknowns for bar chart
fe.chisq <- fe.data[-which(fe.data$Taxonomy == "Unknown"), ]
#--Remove rare classes
fe.chisq <- fe.chisq[which(!fe.chisq$Taxonomy %in% rare.fe), ]
fe.chisq$elevation <- as.factor(fe.chisq$elevation)
fe.chisq$Taxonomy <- as.factor(fe.chisq$Taxonomy)
fe.chisq$Taxonomy <- factor(fe.chisq$Taxonomy,
                            levels = c('Dothideomycetes', 'Leotiomycetes',
                                       'Eurotiomycetes', 'Sordariomycetes'))

#--Bar graph
tax.fe <- ggplot(data = fe.chisq, 
             aes(x = elevation,
                 fill = Taxonomy)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of \n sequences per class") +
  xlab("Elevation (m)") +
  theme_bw() +
  #ggtitle("Proportion of Classes by Topography") +
  #scale_fill_brewer(palette = "Blues") +
  scale_fill_manual(values = c('#9ecae1', '#6baed6', '#2171b5', '#084594')) +
  guides(fill=guide_legend(title=NULL)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        ,axis.text = element_text(size=24, colour = 'black'),
        axis.title = element_text(size = 24), legend.position = "bottom",
        legend.text=element_text(size=14), 
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x = element_text(angle = 45, hjust = 1))

tax.fe

#--save plot to figures folder
ggsave('Figure_2e.pdf', plot = tax.fe, device = 'pdf', path = fig.dir,
       height = 6.85, width = 9.75)
