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
elvgr.em <- table(em.data$group,
                  em.data$taxonomy_class)
#--Remove taxa with equal to or less than 1 occurence
rare.em <- "Fungi incertae sedis"
elvgr.em <- elvgr.em[ ,which(!colnames(elvgr.em) %in% rare.em)]

#-----------------------------------------------------------------------------------------
# Chi-square test
#-----------------------------------------------------------------------------------------
#--Chi square test
em.chi <- chisq.test(elvgr.em)

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
#--Bar graph
tax.em <- ggplot(data = em.chisq, 
          aes(x = group,
          fill = taxonomy_class)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  xlab("Elevation group") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Greys") +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)

#--save plot to figures folder
ggsave('Figure4a.jpeg', plot = tax.em, device = 'jpeg', path = fig.dir)

#=========================================================================================
# Endophytic fungi
#=========================================================================================
#-----------------------------------------------------------------------------------------
# Clean-up data
#-----------------------------------------------------------------------------------------
#--Make count table of em classes
elvgr.fe <- table(fe.data$Group,
                  fe.data$Taxonomy)
#--For removing those classes with only 1 occurence 
rare.fe <- c("Lecanoromycetes", "Orbiliomycetes", "Tremellomycetes")
elvgr.fe <- elvgr.fe[,!colnames(elvgr.fe) %in% rare.fe]

#-----------------------------------------------------------------------------------------
# Chi-square test
#-----------------------------------------------------------------------------------------
#--Chi square test
fe.chi <- chisq.test(elvgr.fe)
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
#--Bar graph
tax.fe <- ggplot(data = fe.chisq, 
             aes(x = Group,
                 fill = Taxonomy)) + 
  geom_bar(position = "fill") + 
  ylab("Proportion of sequences per class") +
  xlab("Elevation group") +
  #ggtitle("Proportion of Classes by Topography") +
  scale_fill_brewer(palette = "Greys") +
  guides (fill=guide_legend(title=NULL)) +
  theme_classic(base_size = 12)

#--save plot to figures folder
ggsave('Figure4b.jpeg', plot = tax.fe, device = 'jpeg', path = fig.dir)
