library(raster)
library(rgdal)

#PLEASE READ THE INTERNAL DOCUMENTATION IN calcPEThamon.r FOR MORE INFORMATION ON THE CALCULATION

#SET THE WORKING DIRECTORYhamon
setwd('/users/stephaniemcafee/documents/PET/ak771m/N')

#ADD THE FUNCTION calcNr.r
source('/users/stephaniemcafee/documents/PETscripts/functions/calcNr.r')

#read in a temperature raster and create a file with the latitude of each point in each grid box -- convert that to radians
tmp = raster('/users/stephaniemcafee/documents/data/snapcruts31/tas50_09/tas_mean_C_cru_TS31_01_1950.tif')

#reproject raster to lat/long
xy = projectExtent(tmp,crs="+proj=longlat +datum=WGS84")

#find the resolution of latitudes from the grid cell center of the first row to the grid cell center of the last row
yres = (yFromRow(xy,1)-yFromRow(xy,nrow(xy)))/(nrow(xy)-1)

#the first row is the highest latitude and it works down from there; make a matrix of those values
val = matrix(NA,nrow(xy),ncol(xy))
val[1,] = yFromRow(xy,1)
for (i in 2:nrow(xy)) {
	val[i,] = val[i-1,]-yres}

#apply those values to a projected raster
lat = raster(extent(tmp),nrows = nrow(tmp),ncols = ncol(tmp),crs = projection(tmp))
lat = setValues(lat,val)

#mask the raster
lat = mask(lat,tmp)

#clean up
rm(xy,val,i,yres,tmp)

#convert lat from degrees to radian
lat = lat*(pi/180)

#set the min and max values
lat = setMinMax(lat)

#trim the lat raster
lat = trim(lat)


#MAKE A VECTOR OF MONTH INDICATORS
month = c('01','02','03','04','05','06','07','08','09','10','11','12')

#MAKE A VECTOR THAT DESCRIBES THE NUMBER OF DAYS IN EACH MONTH
ndm <- matrix(c(31,28,31,30,31,30,31,31,30,31,30,31))


#USE THAT MATRIX TO PULL TOGETHER A MATRIX LISTING THE FIRST AND LAST DAY OF EACH MONTH
d <- matrix(NA,12,2)
for (m in 1:12) {
d[m,2] <- sum(c(ndm[1:m]))   #Sum up number of days in months to find last day of month
} #close that loop
d[,1] <- (d[,2]-ndm)+1 #Subtract number of days in a month from the end date and add one to get Julian Day of the first day of month 
rm(ndm,m) #Delete ndm b/c not needed any more; clean up looping variable
D = apply(d,1,mean)
rm(d)
#This version of the file calculates N for the mean of the days of the year in each month;
for (m in 1:12) {  #month loop
    j=D[m]
    N = calcNr(j)
    rm(j)
       #Write the output data to a geotif
        writeRaster(N,filename=paste('ak_N_771m_',month[m],'.tif',sep=''),format = 'GTiff',options='COMPRESS=LZW',datatype='FLT4S',overwrite=T)
        rm(N)
} #close the month loop

