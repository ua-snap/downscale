require(raster)
require(sp)
require(maptools)
require(rgdal)

# -------------------------------------------------------------------------------------------------------
# surface_downwelling_shortwave_flux_in_air

# alias: surface_downwelling_shortwave_flux

# The surface called "surface" means the lower boundary of the atmosphere. "shortwave" means shortwave radiation. 
# Downwelling radiation is radiation from above. It does not mean "net downward". 
# Surface downwelling shortwave is the sum of direct and diffuse solar radiation incident on the surface, and is sometimes called "global radiation". 
# When thought of as being incident on a surface, a radiative flux is sometimes called "irradiance". 
# In addition, it is identical with the quantity measured by a cosine-collector light-meter and sometimes called "vector irradiance". 
# In accordance with common usage in geophysical disciplines, "flux" implies per unit area, called "flux density" in physics.
# -------------------------------------------------------------------------------------------------------

# read in the dat file downloaded from the CRU TS20 website
sunp <- read.table("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/original_data/TS20/grid_10min_sunp.dat")

# change the colnames to something more useable
colnames(sunp) <- c("lat","lon","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# create a shapefile
sunp.pts <- SpatialPoints(coordinates(sunp[,2:1]))

sunp.spdf <- SpatialPointsDataFrame(sunp.pts, as.data.frame(sunp))

colNums <- 3:ncol(sunp)

template <- raster("/Data/Base_Data/Climate/Canada/CRU_Climatology/Temperature/cru_tmean_1/")
values(template) <- NA

# this loop 
for(i in colNums){
	print(i)
	
	if(nchar(i) == 1){ d=paste("0",i-2, sep="")	}else{ d=i-2 }
	
	rasterize(sunp.spdf, template, field=colnames(sunp)[i] ,filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/output_rasters/sunp_cru_10min_",d,"_1961_1990.tif", sep=""), overwrite=T)

}

####################################################################################################################
# read in the NetCDF as a stack
cld.ts31 <- stack("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/original_data/TS31/cru_ts_3_10.1901.2009.cld.dat.nc")
cld.ts31.tmp <- brick("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/original_data/TS31/cru_ts_3_10.1901.2009.cld.dat.nc")

# read in the cru climatology
# this is in SUNP percentage units
sunp.ts20.clim <- stack(list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/output_rasters", pattern="tif", full.names=TRUE))

# all this does is bring the layernames in as something more descript.  The stack layernames are not as good as the brick
# so I give the stack the brick names
cld@layernames <- cld.ts31.tmp@layernames
rm(cld.ts31.tmp) # remove the unneeded brick

# AKCanada Extent for spatial subsetting
AKCanada.akalb <- raster("/workspace/Shared/Michael/ArcGIS_Shared/clipper/alaska_west_canada_clipper_forDownscaling.tif")

# project the extent 
AKCanada.wgs84 <- projectExtent(AKCanada.akalb, projection(cld.ts31))

# grab the cell ids at the intersection of the new extent and the original data
cld.ts31.desiredCells <- cellsFromExtent(cld.ts31, extent(AKCanada.wgs84), expand=FALSE)

# extract those cells and create a new stack 
cld.ts31.desiredCells.e <- cld[cld.ts31.desiredCells,drop=F]

#####  TS20 climatology data  #########################################################################
# grab the cell ids at the intersection of the new extent and the original data
sunp.ts20.clim.desiredCells <- cellsFromExtent(sunp.ts20.clim, extent(AKCanada.wgs84), expand=FALSE)
# now spatially subset those cells to a new stack
sunp.ts20.clim.desiredCells.e <- sunp.ts20.clim[sunp.ts20.clim.desiredCells,drop=F]

# would be a good idea to do the sunshine percent to cloud percent conversion
# I think this is working correctly now it is units of percent without the percent sign
cld.ts20.clim.desiredCells.e <- 100 - ts20.clim.desiredCells.e 

# temporally subset the stack object to be 1961-1990 for all 12 months
# these numbers were pre-calculated outside of the code
cld.ts31.subset <- subset(cld.ts31.desiredCells.e, 721:1080)

# create an empty brick using the extent from the other stack and no layers in it yet.
cld.ts31.clim <- brick(extent(cld.ts31.subset), nrows=nrow(cld.ts31.subset), ncols=ncol(cld.ts31.subset), crs=projection(cld.ts31.subset), nl=1)

# create climatology
for(j in 1:12){
	#monthLayers<-seq(j, nlayers(cld.subset), by= 12)
	cld.clim.tmp <- mean(subset(cld.ts31.subset, seq(j, nlayers(cld.ts31.subset), by=12)))
	cld.ts31.clim <- addLayer(cld.ts31.clim, cld.clim.tmp)
}

yearList <- 1901:2009 

# create anomalies HISTORICAL
# make an empty brick object
cld.ts31.anom <- brick(extent(cld.ts31.desiredCells.e), nrows=nrow(cld.ts31.desiredCells.e), ncols=ncol(cld.ts31.desiredCells.e), crs=projection(cld.ts31.desiredCells.e), nl=nlayers(cld.ts31.desiredCells.e))

# this is the centroid values of the cru ts20 That is used as the output points to
#  interpolate to from the 0.5 degree to the 10 min spatial resolution of the ts20 climatology 
out_xy <- coordinates(cld.ts20.clim.desiredCells.e)

# calculate and write out the anomalies
for(i in 1:12){ # go through a list of months
	print(paste("	MONTH WORKING = ",i))
	# use the sequence command to create an index of all of the same month during the entire timeseries
	monthList <- seq(i, nlayers(cld.ts31.desiredCells.e), 12)
	
	# here we select the month we want from the ts20 climatology
	cld.ts20.clim.current <- subset(cld.ts20.clim.desiredCells.e,i,drop=T)


	# this line will create a 2 character month for filenaming procedures (if needed)
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}
	
	# here we subset the cld
	cld.ts31.clim.current <- subset(cld.ts31.clim, i, drop=T)
	
	count=0 # a counter to iterate through the years
	
	for(j in monthList){
		print(paste("anomalies iter ",j))
		count=count+1
		cld.current <- subset(cld.desiredCells.e, j, drop=T)

		# this line generates a proportional anomaly
		anom <- cld.current/cld.clim.current

		# write out the 0.5 resolution anomalies (temporary file)
		writeRaster(anom, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/anomalies/anomalies_original_resolution/","cld_pct_proportionalAnom_cru_ts31_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
		

		# # # # # # INTERPOLATION SECTION # # # # # # 
		in_xy <- coordinates(anom)
		z_in <- getValues(anom)

		xyz <- cbind(in_xy,z_in)
		xyz.na <- na.omit(xyz)

		in_xy <- xyz.na[,1:2]
		z_in <- xyz.na[,3]

		anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(cld.ts20.clim.desiredCells.e)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(cld.ts20.clim.desiredCells.e),linear=F))

		# transpose the data here
		nc.interp <- t(anom.spline$z)[,nrow(anom.spline$z):1]

		# turn that into a new raster object
		nc.interp.r <- raster(nc.interp, xmn=xmin(cld.ts20.clim.desiredCells.e), xmx=xmax(ts20.clim.desiredCells.e), ymn=ymin(ts20.clim.desiredCells.e), ymx=ymax(ts20.clim.desiredCells.e), crs=projection(ts20.clim.desiredCells.e))
		
		# write out the interpolated to 10min resolution anomalies (temporary file)
		writeRaster(anom, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/anomalies/interpolated_10min/","cld_pct_proportionalAnom_cru_ts31_10min_", month,"_",yearList[count],".tif",sep=""), overwrite=T)

		cld.ts31.downscaled <- cld.ts20.current*anom

		# write out the downscaled cloud raster
		writeRaster(anom, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/anomalies/interpolated_10min/","cld_pct_cru_ts31_10min_", month,"_",yearList[count],".tif",sep=""), overwrite=T)

	}
}



#########  RUN ABOVE AS A TEST!!!  ##################  FIX BELOW!


# here we need to create a filelist that follows the chronology of the timeseries
fileList <- character()
for(i in yearList){
	print(i)
	for(j in 1:12){
		if(nchar(j)<2){ month=paste("0",j,sep="")}else{month=paste(j,sep="")}
		fileList <- append(fileList, paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/anomalies/","cld_anom_cru_ts31_", month,"_",i,".tif",sep=""), after=length(fileList))
	}
}

cld.anom <- stack(fileList)

# lets create a new empty brick using the information from the 10min data
cld.anom.interp <- brick(ts20.clim.desiredCells.e, values=F, nl=1)

# here we need to interpolate the 0.5 resolution anomalies to something that matches the cru ts20 10 min data
# get the centroid xy
# here we create the interpolated anomalies using a spline interpolation down to 10min resolution
for (f in 1:nlayers(cld.anom)){
	print(basename(fileList[f]))
	inFile <- unlist(strsplit(basename(fileList[f]), "_"))
	outFile <- paste(inFile[1],inFile[2], inFile[3], inFile[4], "10min_interp", inFile[5], inFile[6], sep="_")
	# this is input xy and z data to use in interpolation
	in_xy <- coordinates(subset(cld.subset.climatology, f, drop=T))
	z_in <- getValues(subset(cld.subset.climatology, f, drop=T))

	xyz <- cbind(in_xy,z_in)
	xyz.na <- na.omit(xyz)

	in_xy <- xyz.na[,1:2]
	z_in <- xyz.na[,3]

	# this is the centroid values of the cru ts20
	out_xy <- coordinates(ts20.clim.desiredCells.e)

	anom.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(ts20.clim.desiredCells.e)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(ts20.clim.desiredCells.e),linear=F))

	# transpose the data here
	nc.interp <- t(anom.spline$z)[,nrow(anom.spline$z):1]

	# turn that into a new raster object
	nc.interp.r <- raster(nc.interp, xmn=xmin(ts20.clim.desiredCells.e), xmx=xmax(ts20.clim.desiredCells.e), ymn=ymin(ts20.clim.desiredCells.e), ymx=ymax(ts20.clim.desiredCells.e), crs=projection(ts20.clim.desiredCells.e))
	# write out the new raster here
	writeRaster(nc.interp.r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/anomalies/interpolated_10min/", outFile, sep=""), overwrite=TRUE)

	# and add that layer to a new stack as well
	cld.anom.interp <- addLayer(cld.anom.interp, nc.interp.r)
}


# here I am going to downscale the historical data to 10 min (CLOUD COVER)

for(m in 1:12){
	monthList <- seq(i, nlayers(cld.anom.interp), 12)
	print(paste("working on: ", f,sep=""))

	cld.anom.interp.current <- subset(cld.anom.interp, monthList)

	cld.ts20.clim <- subset()



}

# next loop will involve some calculation between the new cloud cover data and the 

	# next we need to turn the girr to nirr
	# nirr = girr * (0.251 + (0.509*(1.0 - clds/100.0)));
	nirr = cld.subset.v * (0.251 + (0.509*(1.0 - clds/100.0)))

3. CALCULATE THE HISTORICAL ANOMALIES
4. CALCULATE THE FUTURE ANOMALIES
5. INTERPOLATE THE HISTORICAL / FUTURE ANOMALIES TO MATCH THE TS20
6. MULTIPLY/ADD THE ANOMALIES TO THE CRU TS20 CLIMATOLOGY 


* WHAT ABOUT THE UNITS?



# get the values of the stack
cld.subset.v <- getValues(cld.subset)

# apply the function from Dave to create a nirr value using apply()
nirr = cld.subset.v * (0.251 + (0.509*(1.0 - cld.subset.v/100.0)))

# now we return the modified values to the stack.
values(cld.subset2) <- cld.subset.v

# what we have produced in the above code is the NIRR (essentially the same as RSDS accd to Fengming/Dave)
if (clds > -0.1) {
	nirr = cld.subset.v * (0.251 + (0.509*(1.0 - clds/100.0))) # basically percent cloudcover
	}else { 
		nirr = -999.9; }


