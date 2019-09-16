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
require( raster )

#SET THE WORKING DIRECTORY
setwd('/workspace/UA/malindgren/repos/downscale/snap_scripts/old_scripts/tem_iem_older_scripts_april2018/tem_inputs_iem/FromSteph_and_Dave_Radiation')

#ADD THE FUNCTION calcRa.r
source('./calcRa.r')

fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/NCAR-CCSM4/rcp45/tas/tas_mean_C_ar5_NCAR-CCSM4_rcp45_01_2006.tif'
rst = raster( fn )
#READ IN THE FILE THAT CONTAINS INFORMATION ABOUT LATITUDE IN RADIANS FOR EACH GRID CELL
#skip=6: This skips the first 6 lines of header that provide basic metadata for the data (# rows, # cols, xllcorner, yll corner cellsize,a nd the no data value)
#na.strings = '-9999"; this defines the value -9999 as NA for easier indexing
# lat <- as.matrix(read.table('/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/radiance/radians.txt',skip=6,na.strings='-9999'))
lat <- as.matrix(read.table('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/fix_rsds/radians_ml_new_version.txt',skip=6,na.strings='-9999'))
#FIND THE SIZE OF THE LATITUDE MATRIX -- WILL NEED FOR PRE-ALLOCATING SUBSEQUENT MATRICES 
nrow <- nrow(lat)  #1320
ncol <- ncol(lat)  #2015

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
    ra <- matrix(NA,nrow,ncol)
    #Also pre-allocate a temporary array to hold daily values of Ra
    temp <- array(NA,c(length(md),nrow,ncol))
    for (t in 1:length(md))    {  #loop over the number of days in each month
       r <- matrix(NA,nrow,ncol) #Make a matrix for  each day
       #Run the function  for each day
  	   r <- calcRa(md[t],lat,nrow,ncol)
 	     #Slot the daily matrix r into array temp and delete r
       temp[t,,] = r
 	     rm(r)
    }      #End looping over number of days per month
    #Average the daily Ra values into a single month value
    ra <- apply(temp,c(2,3),mean, na.rm = TRUE)
    ra <- round(ra,3)
    rm(temp)   #clean up the temporary matrix
    #Set up a geotiff file to store data
    # out <- raster(ra,xmn=-2173223.206087799, xmx=1728212.22687, ymn=108069.78588, ymx = 2748069.78588, crs = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    out <- raster(ra, crs=crs(rst))
    extent(out) <- extent(rst) 
                                                                                                              # +proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
    #Write data to the file
    writeRaster(out,filename=paste('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/fix_rsds/smcafee_v2/','ak_Ra_',month[m],'.tif',sep=''),format='GTiff', options='COMPRESS=LZW', datatype='FLT4S', overwrite=T)
    rm(ra,out)  #delete files that will be regenerated for the next month
} #close month loop
