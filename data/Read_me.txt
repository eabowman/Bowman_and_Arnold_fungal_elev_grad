Explanations of column names and contents for each data file.

SCM_EM_div_pooled_by_landscape.csv
##################################
elevation: meters above sea level
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
species.richness: number of unique OTUs per site and topographic position
n: total number of isolates per site and topographic position
fisher.alpha: diversity calculated as fisher's alpha using species.richness and n

SCM_EM_div_pooled_by_site.csv
#############################
elevation: meters above sea level
species.richness: species richness, the number of unique OTUs per site
n:  total number of isolates per site
fisher.alpha: diversity calculated as fisher's alpha using species.richness and n

SCM_EM_otu_based_site_x_species_matrix.csv
##########################################
tree_number: sample unit identifier
otu1-99: occurrence of OTUs at each sample unit
elevation: meters above sea level
species: tree species communities sampled
ph.su: pH measure in standard units
EC.ds.m: electrical conductivity in deciSiemens per meter
ca.ppm: calcium in parts per million
mg: magnesium in parts per million
na.ppm: sodium in parts per million
k.ppm: potassium in parts per million
zn.ppm: zinc in parts per million
fe.ppm: iron in parts per million
mn.ppm: manganese in parts per million
cu.ppm: copper in parts per million
ni.ppm: nickel in parts per million
po4.p.ppm: phosphate in parts per million
so4.s.ppm: sulfate in parts per million
b.ppm: boron in parts per million
esp: exchangeable sodium percentage
cec.meq.100g: cation exchange capacity in milliequivalents per 100 grams of soil
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
northing: UTM coordinate system
easting: UTM coordinate system
zone: zome UTM coordinates in, zone 12
p.c: plant community, abbreviated
annual.prec.cm: annual precipitation in centimeters
average.temp.warm.quarter: average temperature in the warmest quarter of the year
average.temp.cold.quarter: average temperature in the coldest quarter of the year
plant.comm: plant community, full names

SCM_EM_per_tree_stats.csv
#########################
tree_number: sample unit identifier
site: site name
elevation.m: meters above sea level
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
p.c: plant community, abbreviated
dbh_cm: diameter breast height of trees in centimeters
morphotypes: number of root tips per morphotype group (i.e., EM root tips grouped together
based on physical characteristics of the mantle)
colonized_rt_tips: total EM colonized root tipes per tree
root_weight.g: dry weight of roots from each root core
species_richness: number of unique OTUs per tree
s: number of unique OTUs per tree
n:  total number of isolates per tree
fisher_alpha: diversity calculated as fisher's alpha using species_richness/s and n
abundance:  abundance of EM per tree calculated as colonized_rt_tips/root_weight.g

SCM_EM_raw.csv
##############
otu.97: EM OTU grouped based on 97% sequence similarity
otu.95: EM OTU grouped based on 95% sequence similarity
otu.90: EM OTU grouped based on 90% sequence similarity
sequencing_number: ID for each sample for sequencing
site: site name
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
elevation: meters above sea level
species: tree species communities sampled
microsite: root core location relative to tree; 3 root cores were collected per tree
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
sequence: DNA ITS sequence
taxonomy_class: taxonomic classification at class level
taxonomy_order: taxonomic classification at order level
taxonomy_genus: taxonomic classification at genus level
taxonomy_epithet: taxonomic classification epithet
original_sample_number: original sample number
physical.description: description of morphotype
abundance: number of colonized root tips per sequence (i.e., number of colonized root tips
per morphotype that that sequence came from)
comments: comments on sequences

SCM_EM_root_based_site_x_species_matrix.csv
###########################################
tree_number: sample unit identifier
otu1-99: occurrence of morphotype at each sample unit
elevation: meters above sea level
species: tree species communities sampled
ph.su: pH measure in standard units
EC.ds.m: electrical conductivity in deciSiemens per meter
ca.ppm: calcium in parts per million
mg: magnesium in parts per million
na.ppm: sodium in parts per million
k.ppm: potassium in parts per million
zn.ppm: zinc in parts per million
fe.ppm: iron in parts per million
mn.ppm: manganese in parts per million
cu.ppm: copper in parts per million
ni.ppm: nickel in parts per million
po4.p.ppm: phosphate in parts per million
so4.s.ppm: sulfate in parts per million
b.ppm: boron in parts per million
esp: exchangeable sodium percentage
cec.meq.100g: cation exchange capacity in milliequivalents per 100 grams of soil
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
northing: UTM coordinate system
easting: UTM coordinate system
zone: zome UTM coordinates in, zone 12
p.c: plant community, abbreviated
annual.prec.cm: annual precipitation in centimeters
average.temp.warm.quarter: average temperature in the warmest quarter of the year
average.temp.cold.quarter: average temperature in the coldest quarter of the year
plant.comm: plant community, full names

SCM_environmental_data.csv
##########################
elevation: meters above sea level
tree_number: sample unit identifier
species: tree species communities sampled
ph.su: pH measure in standard units
EC.ds.m: electrical conductivity in deciSiemens per meter
ca.ppm: calcium in parts per million
mg: magnesium in parts per million
na.ppm: sodium in parts per million
k.ppm: potassium in parts per million
zn.ppm: zinc in parts per million
fe.ppm: iron in parts per million
mn.ppm: manganese in parts per million
cu.ppm: copper in parts per million
ni.ppm: nickel in parts per million
po4.p.ppm: phosphate in parts per million
so4.s.ppm: sulfate in parts per million
b.ppm: boron in parts per million
esp: exchangeable sodium percentage
cec.meq.100g: cation exchange capacity in milliequivalents per 100 grams of soil
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
northing: UTM coordinate system
easting: UTM coordinate system
zone: zome UTM coordinates in, zone 12
p.c: plant community, abbreviated

SCM_FE_div_pooled_by_landscape.csv
##################################
elevation: meters above sea level
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
species.richness: number of unique OTUs per site and topographic position
n: total number of isolates per site and topographic position
fisher.alpha: diversity calculated as fisher's alpha using species.richness and n

SCM_FE_div_pooled_by_site.csv
#############################
elevation: meters above sea level
n:  total number of isolates per site
s: species richness, the number of unique OTUs per site
fisher_alpha: diversity calculated as fisher's alpha using s and n
abundance: isolated frequency by site calculated as number of isolates with fungal growth
divided by 96*6 (i.e., 96 needle segments plated per tree and 6 trees per site)

SCM_FE_per_tree_stats.csv
#########################
tree_number: sample unit identifier
endophyte_collection_date: date needles collected; NA indicates that we were unable to
obtain needles from that particular tree
elevation.m: meters above sea level
site: site name
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
dbh.cm: diameter breast height in centimeters
abundance: isolation frequency by tree calculated as number of isolates with fungal growth
per 96 needle segments plated
species_richness: species richness, the number of unique OTUs per tree
s: species richness, the number of unique OTUs per tree
n: total number of isolates per site
fisher_alpha: diversity calculated as fisher's alpha using s/species_richness and n

SCM_FE_raw.csv
##############
Sequencing_number: identifier for sequencing isolates
OTU: OTUs clustered at 95% sequence similarity
Tree_number: sample unit identifier
site: site name
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
species: species of trees needles collected from
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
Sequence: ITS sequence
Taxonomy: taxonomic classification at the class level

SCM_FE_site_x_species_matrix.csv
################################
tree_number: sample unit identifier
OTU1-OTU9: occurrence of OTU at each sample unit
elevation: meters above sea level
species: tree species communities sampled from
ph.su: pH measure in standard units
EC.ds.m: electrical conductivity in deciSiemens per meter
ca.ppm: calcium in parts per million
mg: magnesium in parts per million
na.ppm: sodium in parts per million
k.ppm: potassium in parts per million
zn.ppm: zinc in parts per million
fe.ppm: iron in parts per million
mn.ppm: manganese in parts per million
cu.ppm: copper in parts per million
ni.ppm: nickel in parts per million
po4.p.ppm: phosphate in parts per million
so4.s.ppm: sulfate in parts per million
b.ppm: boron in parts per million
esp: exchangeable sodium percentage
cec.meq.100g: cation exchange capacity in milliequivalents per 100 grams of soil
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
northing: UTM coordinate system
easting: UTM coordinate system
zone: zome UTM coordinates in, zone 12
p.c: plant community, abbreviated
annual.prec.cm: annual precipitation in centimeters
average.temp.warm.quarter: average temperature in the warmest quarter of the year
average.temp.cold.quarter: average temperature in the coldest quarter of the year
plant.comm: plant community, full names

SCM_soil.csv
############
elevation: meters above sea level
tree_number: sample unit identifier
species: tree species communities sampled from
ph.su: pH measure in standard units
EC.ds.m: electrical conductivity in deciSiemens per meter
ca.ppm: calcium in parts per million
mg: magnesium in parts per million
na.ppm: sodium in parts per million
k.ppm: potassium in parts per million
zn.ppm: zinc in parts per million
fe.ppm: iron in parts per million
mn.ppm: manganese in parts per million
cu.ppm: copper in parts per million
ni.ppm: nickel in parts per million
po4.p.ppm: phosphate in parts per million
so4.s.ppm: sulfate in parts per million
b.ppm: boron in parts per million
esp: exchangeable sodium percentage
cec.meq.100g: cation exchange capacity in milliequivalents per 100 grams of soil
topography: topographic position of trees, convergent positions are low lying positions
where nutrients and water pool, divergent positions are higher positons where nutrients
and water flow away from
group: elevation group, Low elev. (< 2300 masl), High elev. (> 2300 masl)
northing: UTM coordinate system
easting: UTM coordinate system
zone: zome UTM coordinates in, zone 12
p.c: plant community, abbreviated
annual.prec.cm: annual precipitation in centimeters
average.temp.warm.quarter: average temperature in the warmest quarter of the year
average.temp.cold.quarter: average temperature in the coldest quarter of the year