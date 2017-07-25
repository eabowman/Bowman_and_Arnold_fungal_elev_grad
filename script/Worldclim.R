## Script created by Liz Bowman August 9, 2016
## for downloading and isolating climate data from site coordinates

#=========================================================================================
# Climate data
#=========================================================================================

#--First, you need to go to http://www.worldclim.org/bioclim and download climate data 
#--for the quadrant where you sites are. Click on Download tab at top of page and download
#--past, current, or future climate data.
#--Be sure to download by tile (there should be a link at the top of the page for it). 
#--Bioclim data downloads a ton of different variables (http://www.worldclim.org/bioclim)

#-----------------------------------------------------------------------------------------
# download climate data from Worldclim using raster package
#-----------------------------------------------------------------------------------------

#--install and set library
#--uncomment install packages function if you don't have this package installed
#install.packages ("raster")
library (raster)

# << BIOCLIM VARIABLES >> ----------------------------------------------------------------
#--load world climate data: for 19 bio variables
#--no need to change this if you downloaded tile 12
data.1 <- raster (paste0 (clim.dir, 'bio_12/bio1_12.bil'),native = F)
data.2 <- raster (paste0 (clim.dir, 'bio_12/bio2_12.bil'),native = F)
data.3 <- raster (paste0 (clim.dir, 'bio_12/bio3_12.bil'),native = F)
data.4 <- raster (paste0 (clim.dir, 'bio_12/bio4_12.bil'),native = F)
data.5 <- raster (paste0 (clim.dir, 'bio_12/bio5_12.bil'),native = F)
data.6 <- raster (paste0 (clim.dir, 'bio_12/bio6_12.bil'),native = F)
data.7 <- raster (paste0 (clim.dir, 'bio_12/bio7_12.bil'),native = F)
data.8 <- raster (paste0 (clim.dir, 'bio_12/bio8_12.bil'),native = F)
data.9 <- raster (paste0 (clim.dir, 'bio_12/bio9_12.bil'),native = F)
data.10 <- raster (paste0 (clim.dir, 'bio_12/bio10_12.bil'),native = F)
data.11 <- raster (paste0 (clim.dir, 'bio_12/bio11_12.bil'),native = F)
data.12 <- raster (paste0 (clim.dir, 'bio_12/bio12_12.bil'),native = F)
data.13 <- raster (paste0 (clim.dir, 'bio_12/bio13_12.bil'),native = F)
data.14 <- raster (paste0 (clim.dir, 'bio_12/bio14_12.bil'),native = F)
data.15 <- raster (paste0 (clim.dir, 'bio_12/bio15_12.bil'),native = F)
data.16 <- raster (paste0 (clim.dir, 'bio_12/bio16_12.bil'),native = F)
data.17 <- raster (paste0 (clim.dir, 'bio_12/bio17_12.bil'),native = F)
data.18 <- raster (paste0 (clim.dir, 'bio_12/bio18_12.bil'),native = F)
data.19 <- raster (paste0 (clim.dir, 'bio_12/bio19_12.bil'),native = F)

#<< Individual variables >> --------------------------------------------------------------
#--load world climate data: Precipitation, Tmax, Tmin, Tmean
#--there is a file for each month
#--no need to change this if you downloaded for tile 12
#--Precipitation
prec.1 <- raster (paste0 (clim.dir, 'prec_12/prec1_12.bil'),native = F)
prec.2 <- raster (paste0 (clim.dir, 'prec_12/prec2_12.bil'),native = F)
prec.3 <- raster (paste0 (clim.dir, 'prec_12/prec3_12.bil'),native = F)
prec.4 <- raster (paste0 (clim.dir, 'prec_12/prec4_12.bil'),native = F)
prec.5 <- raster (paste0 (clim.dir, 'prec_12/prec5_12.bil'),native = F)
prec.6 <- raster (paste0 (clim.dir, 'prec_12/prec6_12.bil'),native = F)
prec.7 <- raster (paste0 (clim.dir, 'prec_12/prec7_12.bil'),native = F)
prec.8 <- raster (paste0 (clim.dir, 'prec_12/prec8_12.bil'),native = F)
prec.9 <- raster (paste0 (clim.dir, 'prec_12/prec9_12.bil'),native = F)
prec.10 <- raster (paste0 (clim.dir, 'prec_12/prec10_12.bil'),native = F)
prec.11 <- raster (paste0 (clim.dir, 'prec_12/prec11_12.bil'),native = F)
prec.12 <- raster (paste0 (clim.dir, 'prec_12/prec12_12.bil'),native = F)
#--Minimum temperature
Tmin.1 <- raster (paste0 (clim.dir, 'tmin_12/tmin1_12.bil'),native = F)
Tmin.2 <- raster (paste0 (clim.dir, 'tmin_12/tmin2_12.bil'),native = F)
Tmin.3 <- raster (paste0 (clim.dir, 'tmin_12/tmin3_12.bil'),native = F)
Tmin.4 <- raster (paste0 (clim.dir, 'tmin_12/tmin4_12.bil'),native = F)
Tmin.5 <- raster (paste0 (clim.dir, 'tmin_12/tmin5_12.bil'),native = F)
Tmin.6 <- raster (paste0 (clim.dir, 'tmin_12/tmin6_12.bil'),native = F)
Tmin.7 <- raster (paste0 (clim.dir, 'tmin_12/tmin7_12.bil'),native = F)
Tmin.8 <- raster (paste0 (clim.dir, 'tmin_12/tmin8_12.bil'),native = F)
Tmin.9 <- raster (paste0 (clim.dir, 'tmin_12/tmin9_12.bil'),native = F)
Tmin.10 <- raster (paste0 (clim.dir, 'tmin_12/tmin10_12.bil'),native = F)
Tmin.11 <- raster (paste0 (clim.dir, 'tmin_12/tmin11_12.bil'),native = F)
Tmin.12 <- raster (paste0 (clim.dir, 'tmin_12/tmin12_12.bil'),native = F)
#--Maximum temperature
Tmax.1 <- raster (paste0 (clim.dir, 'tmax_12/tmax1_12.bil'),native = F)
Tmax.2 <- raster (paste0 (clim.dir, 'tmax_12/tmax2_12.bil'),native = F)
Tmax.3 <- raster (paste0 (clim.dir, 'tmax_12/tmax3_12.bil'),native = F)
Tmax.4 <- raster (paste0 (clim.dir, 'tmax_12/tmax4_12.bil'),native = F)
Tmax.5 <- raster (paste0 (clim.dir, 'tmax_12/tmax5_12.bil'),native = F)
Tmax.6 <- raster (paste0 (clim.dir, 'tmax_12/tmax6_12.bil'),native = F)
Tmax.7 <- raster (paste0 (clim.dir, 'tmax_12/tmax7_12.bil'),native = F)
Tmax.8 <- raster (paste0 (clim.dir, 'tmax_12/tmax8_12.bil'),native = F)
Tmax.9 <- raster (paste0 (clim.dir, 'tmax_12/tmax9_12.bil'),native = F)
Tmax.10 <- raster (paste0 (clim.dir, 'tmax_12/tmax10_12.bil'),native = F)
Tmax.11 <- raster (paste0 (clim.dir, 'tmax_12/tmax11_12.bil'),native = F)
Tmax.12 <- raster (paste0 (clim.dir, 'tmax_12/tmax12_12.bil'),native = F)

#--import csv data frame with site name, lat, long
#--make sure column names are site, lat, and long for next set of commands
sites <- read.csv(paste0(dat.dir,'Insert_file_name_here'), as.is = T, header = T)

#--combine lat and long
xy <- cbind (sites$long, sites$lat)
#--create object of class SpatialPoints from coordinates
sp <- SpatialPoints(xy)
#--extract information from downloaded files
#--Bio10 and Bio11 are mean temp of warmest quarter and "" of coldest quarter respectively
#--Bio12 is annual precipitation given in millimeters
sites_cli <- as.data.frame (cbind (sites,
             "BIO10" = extract (data.10, sp, method = 'bilinear'),
             "BIO11" = extract (data.11, sp, method = 'bilinear'),
             "BIO12" = extract (data.12, sp, method = 'bilinear')))
#--Divide Bio10 and Bio11 by 10 as Worldclim stores the data multiplied by 10
sites_cli$BIO10_red <- sites_cli$BIO10 / 10
sites_cli$BIO11_red <- sites_cli$BIO11 / 10

#--temperature min, temperature max, and precipitation
sites_cli2 <- as.data.frame (cbind (sites,
                                   "Prec1" = extract (prec.1, sp, method = 'bilinear'),
                                   "Prec2" = extract (prec.2, sp, method = 'bilinear'),
                                   "Prec3" = extract (prec.3, sp, method = 'bilinear'),
                                   "Prec4" = extract (prec.4, sp, method = 'bilinear'),
                                   "Prec5" = extract (prec.5, sp, method = 'bilinear'),
                                   "Prec6" = extract (prec.6, sp, method = 'bilinear'),
                                   "Prec7" = extract (prec.7, sp, method = 'bilinear'),
                                   "Prec8" = extract (prec.8, sp, method = 'bilinear'),
                                   "Prec9" = extract (prec.9, sp, method = 'bilinear'),
                                   "Prec10" = extract (prec.10, sp, method = 'bilinear'),
                                   "Prec11" = extract (prec.11, sp, method = 'bilinear'),
                                   "Prec12" = extract (prec.12, sp, method = 'bilinear'),
                                   "Tmax1" = extract (Tmax.1, sp, method = 'bilinear'),
                                   "Tmax2" = extract (Tmax.2, sp, method = 'bilinear'),
                                   "Tmax3" = extract (Tmax.3, sp, method = 'bilinear'),
                                   "Tmax4" = extract (Tmax.4, sp, method = 'bilinear'),
                                   "Tmax5" = extract (Tmax.5, sp, method = 'bilinear'),
                                   "Tmax6" = extract (Tmax.6, sp, method = 'bilinear'),
                                   "Tmax7" = extract (Tmax.7, sp, method = 'bilinear'),
                                   "Tmax8" = extract (Tmax.8, sp, method = 'bilinear'),
                                   "Tmax9" = extract (Tmax.9, sp, method = 'bilinear'),
                                   "Tmax10" = extract (Tmax.10, sp, method = 'bilinear'),
                                   "Tmax11" = extract (Tmax.11, sp, method = 'bilinear'),
                                   "Tmax12" = extract (Tmax.12, sp, method = 'bilinear'),
                                   "Tmin1" = extract (Tmin.1, sp, method = 'bilinear'),
                                   "Tmin2" = extract (Tmin.2, sp, method = 'bilinear'),
                                   "Tmin3" = extract (Tmin.3, sp, method = 'bilinear'),
                                   "Tmin4" = extract (Tmin.4, sp, method = 'bilinear'),
                                   "Tmin5" = extract (Tmin.5, sp, method = 'bilinear'),
                                   "Tmin6" = extract (Tmin.6, sp, method = 'bilinear'),
                                   "Tmin7" = extract (Tmin.7, sp, method = 'bilinear'),
                                   "Tmin8" = extract (Tmin.8, sp, method = 'bilinear'),
                                   "Tmin9" = extract (Tmin.9, sp, method = 'bilinear'),
                                   "Tmin10" = extract (Tmin.10, sp, method = 'bilinear'),
                                   "Tmin11" = extract (Tmin.11, sp, method = 'bilinear'),
                                   "Tmin12" = extract (Tmin.12, sp, method = 'bilinear')))

#--Take average of precipitation and add it to the sites.cli file
sites_cli$prec <- sites_cli2$Prec1 + sites_cli2$Prec2 + sites_cli2$Prec3 +
  sites_cli2$Prec4 + sites_cli2$Prec5 + sites_cli2$Prec6 + sites_cli2$Prec7 +
  sites_cli2$Prec8 + sites_cli2$Prec9 + sites_cli2$Prec10 + sites_cli2$Prec11 +
  sites_cli2$Prec12
#--Take average of tmin and tmax; add it to the sites.cli file
sites_cli$Tmin <- ((sites_cli2$Tmin1 + sites_cli2$Tmin2 + sites_cli2$Tmin3 +
  sites_cli2$Tmin4 + sites_cli2$Tmin5 + sites_cli2$Tmin6 + sites_cli2$Tmin7 +
  sites_cli2$Tmin8 + sites_cli2$Tmin9 + sites_cli2$Tmin10 + sites_cli2$Tmin11 +
  sites_cli2$Tmin12)/10)/12
sites_cli$Tmax <- ((sites_cli2$Tmax1 + sites_cli2$Tmax2 + sites_cli2$Tmax3 +
  sites_cli2$Tmax4 + sites_cli2$Tmax5 + sites_cli2$Tmax6 + sites_cli2$Tmax7 +
  sites_cli2$Tmax8 + sites_cli2$Tmax9 + sites_cli2$Tmax10 + sites_cli2$Tmax11 +
  sites_cli2$Tmax12)/10)/12

#--write file
#write.csv(sites_cli,paste0(res.dir,'Insert_new_file_name_here'), row.names = F)
