################################################################################
#This file uses the function calcRa.r to calculate extraterrestrial solar 
#radiation (solar radiation at the top of the atmosphere) for a daily time
#step that is then averaged into monthly.
#NB: This does not need to be done for each year of the time series because Ra 
#is a function of day of year and latitude, neither of which vary between years.
#It then produces a text file of Ra values for each month.  
#NB:  The text files produced here have no header for easier processing in R but
#may need to be modified for usein ArcMap.
################################################################################

#PLEASE READ THE INTERNAL DOCUMENTATION IN calcRa.r FOR MORE INFORMATION ON THE CALCULATION
library(raster)
library(sp)
library(maptools)
library(rgeos)


#SET THE WORKING DIRECTORY
setwd("/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/girr_radiation_cmip3_process/")

outDir<-"/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/girr_radiation_cmip3_process/AKCAN_2km/"

#ADD THE FUNCTION calcRa.r
source("/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin/calcRa.r")

# use this as the template.map
template.map <- raster("/Data/malindgren/cru_november_final/ar5/GFDL-CM3/cld/downscaled/cld_metric_GFDL-CM3_historical_r1i1p1_03_1955.tif")

# here we make only the cells with actual values a new value of zero.  the rest remain out of bounds.
ind <- which(getValues(template.map) > -9999)
values(template.map)[ind] <- 0

#project raster
#template.map.projExt <- projectExtent(template.map, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#template.map <- projectRaster(template.map, template.map.projExt, method="ngb")

# create a new map where the cells values are the latitudes of the cell center
template.map.lat <- raster(matrix(coordinates(template.map)[,1],nrow=nrow(template.map),ncol=ncol(template.map)),xmn=xmin(template.map),xmx=xmax(template.map),ymn=ymin(template.map),ymx=ymax(template.map),crs=projection(template.map))


# take that template map and turn all values that are NA in the template.map back to NA
ind <- which(is.na(getValues(template.map)))
values(template.map.lat)[ind] <- NA

# calculate the radians here:
#  radians = degrees*pi/180
template.map.rad <- (template.map.lat*pi)/180

# this is a testing line that will write out the current file to an ascii grid to see if reading the damn thing in as a table will solve the issues I am having...
writeRaster(template.map.rad, filename=paste(outDir,"template_map_akalb_rad.asc",sep=""), overwrite=TRUE)

# here we read that same file back in as a table where we skip the 6 lines at the top of the arc ascii grid which is a header
template.map.rad <- t(read.table(paste(outDir,"template_map_akalb_rad.asc",sep=""), sep=" ", skip=6))

#MAKE A VECTOR THAT DESCRIBES THE NUMBER OF DAYS IN EACH MONTH
ndm <- matrix(c(31,28,31,30,31,30,31,31,30,31,30,31))

#USE THAT MATRIX TO PULL TOGETHER A MATRIX LISTING THE FIRST AND LAST DAY OF EACH MONTH
n <- matrix(NA,12,2)
for (m in 1:12) {
n[m,2] <- sum(c(ndm[1:m]))   #Sum up number of days in months to find last day of month
} #close that loop
n[,1] <- (n[,2]-ndm)+1 #Subtract number of days in a month from the end date and add one to get Julian Day of the first day of month 
rm(ndm,m) #Delete ndm b/c not needed any more; clean up looping variable

#DEFINE A VECTOR OR MONTH DESIGNATIONS FOR READING IN FILES AND NAMING OUTPUT FILES 
month<-c(1:12)

#RUN THE FUNCTION AND SAVE OUTPUT TO TEXT FILES
#This version of the file calculates Ra for each day within the month and then averages daily radiation to get a monthly value, rather than using a "representative" day.
for (m in 1:12) {  #month loop
    #Identify the Julian Days associated with that month
    md <-n[m,1]:n[m,2]
    #Pre-allocate the Ra matrix
    ra <- matrix(NA,nrow(template.map.rad),ncol(template.map.rad))
    #Also pre-allocate a temporary array to hold daily values of Ra
    temp <- array(NA,c(length(md),nrow(template.map.rad),ncol(template.map.rad)))
    for (t in 1:length(md))    {  #loop over the number of days in each month
       r <- matrix(NA,nrow(template.map.rad),ncol(template.map.rad)) #Make a matrix for  each day
       #Run the function  for each day
       r <- calcRa(md[t],template.map.rad,nrow(template.map.rad),ncol(template.map.rad))
      #Slot the daily matrix r into array temp and delete r
       temp[t,,] = r
        rm(r)
    }      #End looping over number of days per month
    #Average the daily Ra values into a single month value
    ra <- apply(temp,c(2,3),mean, na.rm = TRUE)
    ra <- round(ra,3)
    rm(temp)   #clean up the temporary matrix
    #Set up a geotiff file to store data
    out <- raster(ra,xmn=xmin(template.map),xmx=xmax(template.map),ymn=ymin(template.map), ymx=ymax(template.map), crs=projection(template.map))
    writeRaster(out,filename=paste(outDir,"rad_",month[m],"_1km_iem_akalb.tif",sep=""),options=c('COMPRESS=LZW'), overwrite=TRUE)

    # ra.akalb.projExt <- projectExtent(out, crs=projection(template.map))
    # ra.akalb <- projectRaster(out, ra.akalb.projExt, method='bilinear', filename=paste(outDir,"rad_",month[m],"_1km_iem.tif",sep=""), options=c('COMPRESS=LZW'),overwrite=TRUE)
    #Write data to the file
    #writeRaster(out,filename=paste('ak_Ra_',month[m],'.tif',sep=''),format = 'GTiff',options='COMPRESS=LZW',datatype='FLT4S',overwrite=T)
    rm(ra,out)  #delete files that will be regenerated for the next month
    print(paste(" >>> completed: ",m,sep=""))

} #close month loop
