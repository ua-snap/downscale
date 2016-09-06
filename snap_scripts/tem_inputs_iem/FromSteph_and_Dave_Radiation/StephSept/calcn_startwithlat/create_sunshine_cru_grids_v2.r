require(raster)
require(sp)
require(maptools)
require(rgdal)
require(akima)


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

# PROJ4 Reference projections
akalb <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

# template map used to create the ts20 cru raster data from the input file format
template <- raster("/Data/Base_Data/Climate/Canada/CRU_Climatology/Temperature/cru_tmean_1/")
values(template) <- NA

# the 1km output template for the final spline interpolation
out_template <- raster("/Data/Base_Data/ALFRESCO_formatted/ALFRESCO_Master_Dataset/ALFRESCO_Model_Input_Datasets/AK_CAN_Inputs/Climate/cru_TS31/historical/pr/pr_total_mm_alf_cru_TS31_01_1901.tif")


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

# read in the list of 12 monthly RADIATION (girr) data calculated from Dave McGuires C++ code translated into R by Michael Lindgren
l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/girr_radiation/version2", pattern=".tif", full.names=T)
# now stack that list of files
rad.monthly <- stack(l)

# read in the cru climatology
# this is in SUNP percentage units
sunp.ts20.clim <- stack(list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/output_rasters", pattern="tif", full.names=TRUE))

# all this does is bring the layernames in as something more descript.  The stack layernames are not as good as the brick
# so I give the stack the brick names
cld.ts31@layernames <- cld.ts31.tmp@layernames
rm(cld.ts31.tmp) # remove the unneeded brick

# AKCanada Extent for spatial subsetting
AKCanada.akalb <- raster("/workspace/Shared/Michael/ArcGIS_Shared/clipper/alaska_west_canada_clipper_forDownscaling.tif")

# project the extent 
AKCanada.wgs84 <- projectExtent(AKCanada.akalb, projection(cld.ts31))

# grab the cell ids at the intersection of the new extent and the original data
cld.ts31.desiredCells <- cellsFromExtent(cld.ts31, extent(AKCanada.wgs84), expand=FALSE)

# extract those cells and create a new stack 
cld.ts31.desiredCells.e <- cld.ts31[cld.ts31.desiredCells,drop=F]

#####  TS20 climatology data  #########################################################################
# grab the cell ids at the intersection of the new extent and the original data
sunp.ts20.clim.desiredCells <- cellsFromExtent(sunp.ts20.clim, extent(AKCanada.wgs84), expand=FALSE)
# now spatially subset those cells to a new stack
sunp.ts20.clim.desiredCells.e <- sunp.ts20.clim[sunp.ts20.clim.desiredCells,drop=F]

# would be a good idea to do the sunshine percent to cloud percent conversion
# I think this is working correctly now it is units of percent without the percent sign
cld.ts20.clim.desiredCells.e <- 100 - sunp.ts20.clim.desiredCells.e 

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

# lets reproject the ts20 Climatology into the final projection system of the output
cld.ts20.clim.extent <- projectExtent(cld.ts20.clim.desiredCells.e, crs=akalb)
cld.ts20.clim.desiredCells.e <- projectRaster(cld.ts20.clim.desiredCells.e, cld.ts20.clim.extent, res=res(rad.monthly), method='ngb')

# create anomalies HISTORICAL
# make an empty brick object
cld.ts31.anom <- brick(extent(cld.ts31.desiredCells.e), nrows=nrow(cld.ts31.desiredCells.e), ncols=ncol(cld.ts31.desiredCells.e), crs=projection(cld.ts31.desiredCells.e), nl=nlayers(cld.ts31.desiredCells.e))
cld.ts31.downscaled.brick <- brick(extent(cld.ts20.clim.desiredCells.e), nrows=nrow(cld.ts20.clim.desiredCells.e), ncols=ncol(cld.ts20.clim.desiredCells.e), crs=projection(cld.ts20.clim.desiredCells.e))

# this is the centroid values of the cru ts20 That is used as the output points to
#  interpolate to from the 0.5 degree to the 10 min spatial resolution of the ts20 climatology 
out_xy <- coordinates(cld.ts20.clim.desiredCells.e)
# the output coords for the 1km output
out_1km <- coordinates(out_template)

# calculate and write out the anomalies
for(i in 1:12){ # go through a list of months
	print(paste("	MONTH WORKING = ",i))
	# use the sequence command to create an index of all of the same month during the entire timeseries
	monthList <- seq(i, nlayers(cld.ts31.desiredCells.e), 12)
	
	# here we select the month we want from the ts20 climatology
	cld.ts20.clim.current <- subset(cld.ts20.clim.desiredCells.e,i,drop=T)
	
	# this is where we subset the RADIATION montlies to the current month for the final calculation
	rad.monthly.current <- subset(rad.monthly,i,drop=T)
	# rad.monthly.current.v <- getValues(rad.monthly.current)

	# this line will create a 2 character month for filenaming procedures (if needed)
	if(nchar(i)<2){ month=paste("0",i,sep="")}else{month=paste(i,sep="")}
	
	# here we subset the cld
	cld.ts31.clim.current <- subset(cld.ts31.clim, i, drop=T)

	# here we reproject the data to be in the final output reference system
	cld.ts31.clim.current.extent <- projectExtent(cld.ts31.clim.current, crs=akalb)
	cld.ts31.clim.current <- projectRaster(cld.ts31.clim.current, cld.ts31.clim.current.extent, method='ngb')

	count=0 # a counter to iterate through the years
	
	for(j in monthList){
		print(paste("anomalies iter ",j))
		count=count+1
		cld.ts31.current <- subset(cld.ts31.desiredCells.e, j, drop=T)

		# here we are going to make the final reprojection into the output reference system
		cld.ts31.current.extent <- projectExtent(cld.ts31.current, crs=akalb)
		cld.ts31.current <- projectRaster(cld.ts31.current, cld.ts31.current.extent, method='ngb')

		# this line generates a proportional anomaly
		anom <- cld.ts31.current/cld.ts31.clim.current

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
		anom.interp <- t(anom.spline$z)[,nrow(anom.spline$z):1]

		# turn that into a new raster object
		anom.interp.r <- raster(anom.interp, xmn=xmin(cld.ts20.clim.desiredCells.e), xmx=xmax(cld.ts20.clim.desiredCells.e), ymn=ymin(cld.ts20.clim.desiredCells.e), ymx=ymax(cld.ts20.clim.desiredCells.e), crs=projection(cld.ts20.clim.desiredCells.e))
		
		# write out the interpolated to 10min resolution anomalies (temporary file)
		writeRaster(anom.interp.r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/anomalies/interpolated_10min/","cld_pct_proportionalAnom_cru_ts31_10min_", month,"_",yearList[count],".tif",sep=""), overwrite=T)

		# this is where we put downscale the interpolated proportional anomalies
		cld.ts31.downscaled <- cld.ts20.clim.current*anom.interp.r

		# write out the downscaled cloud raster
		writeRaster(cld.ts31.downscaled, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/downscaled_10min/cld/","cld_pct_cru_ts31_10min_", month,"_",yearList[count],".tif",sep=""), overwrite=T)

		######  INTERPOLATION TO 1km Section #######

		# here we want to interpolate the new downscaled raster to 1km resolution used by TEM
		xyz <- cbind(coordinates(cld.ts31.downscaled), getValues(cld.ts31.downscaled))
		# get rid of the NA's
		xyz.na <- na.omit(xyz)
		xyz=xyz.na

		# here we interpolate to a 1km data set
		downscaled.spline	<- interp(x=xyz[,1],y=xyz[,2],z=xyz[,3],xo=seq(min(out_1km[,1]),max(out_1km[,1]),l=ncol(out_template)), yo=seq(min(out_1km[,2]),max(out_1km[,2]),l=nrow(out_template),linear=F))

		# transpose the data here
		downscaled.interp <- t(downscaled.spline$z)[,nrow(downscaled.spline$z):1]

		# turn that into a new raster object
		downscaled.interp.r <- raster(downscaled.interp, xmn=xmin(out_template), xmx=xmax(out_template), ymn=ymin(out_template), ymx=ymax(out_template), crs=projection(out_template))
		
		# these lines take the values that exceed 100% and truncate them back to 100
		ind <- which(values(downscaled.interp.r) > 100)
		values(downscaled.interp.r)[ind] <- 100

		# write out the downscaled cloud raster
		writeRaster(downscaled.interp.r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/downscaled_1km/cld/","cld_pct_cru_ts31_1km_", month,"_",yearList[count],".tif",sep=""), overwrite=T)

		# now we calculate the nirr from the girr and the cloudcover percent using the code from Dave McGuire
		##############################################################################################################
		#		** this is the c++ version given by Dave.  I am using R to calculate it here
		#		if (clds > -0.1) {
		#			nirr = cld.subset.v * (0.251 + (0.509*(1.0 - clds/100.0))) # basically percent cloudcover
		#			}else { 
		#				nirr = -999.9; }
		##############################################################################################################

		output.v <- getValues(downscaled.interp.r)
		
		#ind <- which(values(downscaled.interp.r) > -0.1)
		oob <- which(values(downscaled.interp.r) <= -0.1)
		
		if(length(oob)>0){
			nirr = rad.monthly.current * (0.251 + (0.509 * (1.0 - downscaled.interp.r/100)))
			values(nirr)[oob] <- -9999
		}else{
			nirr = rad.monthly.current * (0.251 + (0.509 * (1.0 - downscaled.interp.r/100)))
		}

		writeRaster(nirr, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/downscaled_1km/nirr/","nirr_pct_cru_ts31_1km_", month,"_",yearList[count],".tif",sep=""), overwrite=T)
	}
}



