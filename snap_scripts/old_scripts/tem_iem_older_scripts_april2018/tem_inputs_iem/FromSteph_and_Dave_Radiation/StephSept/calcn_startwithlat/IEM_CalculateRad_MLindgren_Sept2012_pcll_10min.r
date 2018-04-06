library(raster)
library(rgdal)

# This script is for calculating the RSDS for the TEM model in the IEM project
# modified code from Steph McAfee, by Michael Lindgren September 14, 2012

# SET THE WORKING DIRECTORY
setwd('/big_storage/malindgren/AIEM/RSDS/girr')

# ADD THE FUNCTION calcNr.r
source('/workspace/UA/malindgren/projects/iem/PHASE2_DATA/code/FromSteph_and_Dave_Radiation/StephSept/calcn_startwithlat/calcNr.R')

# read in a raster layer that has the same extent and projection system as the desired output Rad files
tmp = raster("/big_storage/malindgren/AIEM/RSDS/cru_ts20/pacific_centered_grids/sunshine_cru_ts20_01_1961_1990.tif")

# project that raster to WGS1984 lat/long
xy = tmp #projectExtent(tmp,crs="+proj=longlat +datum=WGS84")

# find the resolution of latitudes from the grid cell center of the first row to the grid cell center of the last row
# this finds lat from the cell# 1 in row 1 and cell# 1 in row 2227 of xy and divides the later by the number of rows 
# between 1 and 2227, which is 2226, thus returning the cell resolution
yres = (yFromRow(xy,1)-yFromRow(xy,nrow(xy)))/(nrow(xy)-1)

lat.ll.m = as.matrix(coordinates(xy)[,2],nrow(xy),ncol(xy))

# Modified version of the lat Calc
lat = raster(extent(tmp), nrows=nrow(tmp), ncols=ncol(tmp), crs=projection(tmp))
lat = setValues(lat, lat.ll.m)

# mask the raster
lat = mask(lat,tmp)

# convert lat from degrees to radian
lat = lat*(pi/180)

# set the min and max values
# lat = setMinMax(lat)

# trim the lat raster
lat = trim(lat)

# MAKE A VECTOR OF MONTH INDICATORS
month = c('01','02','03','04','05','06','07','08','09','10','11','12')

# MAKE A VECTOR THAT DESCRIBES THE NUMBER OF DAYS IN EACH MONTH
ndm <- matrix(c(31,28,31,30,31,30,31,31,30,31,30,31))

# USE THAT MATRIX TO PULL TOGETHER A MATRIX LISTING THE FIRST AND LAST DAY OF EACH MONTH
d <- matrix(NA,12,2)
for (m in 1:12) {
d[m,2] <- sum(c(ndm[1:m]))   # Sum up number of days in months to find last day of month
} # close that loop

d[,1] <- (d[,2]-ndm)+1 # Subtract number of days in a month from the end date and add one to get Julian Day of the first day of month 
rm(ndm,m) # Delete ndm b/c not needed any more; clean up looping variable
D = apply(d,1,mean)
rm(d)
# This version of the file calculates N for the mean of the days of the year in each month;
for (m in 1:12) {  # month loop
    j=D[m]
    N = calcNr(j)
    rm(j)
       # Write the output data to a geotiff
        writeRaster(N,filename=paste("rad_10min_global_pcll_",month[m],'.tif',sep=''),format = 'GTiff',options='COMPRESS=LZW',datatype='FLT4S',overwrite=T)
        rm(N)
} # close the month loop