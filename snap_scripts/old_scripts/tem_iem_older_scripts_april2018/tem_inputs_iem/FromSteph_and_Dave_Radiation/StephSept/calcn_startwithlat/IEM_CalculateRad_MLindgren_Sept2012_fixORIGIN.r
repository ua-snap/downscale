require(maptools)
require(sp)
require(raster)
require(rgeos)


# this little script fixes a small problem with the origin of the RAD data being different than the origin of the downscaled outputs.  this is creating an issue where I 
#  cannot perform any calculation of the 2 maps until the origins match.   

l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/girr_radiation/", pattern=".tif", full.names=T)
# now stack that list of files
rad.monthly <- stack(l)

# this is the template map used to rasterize the data to a new origin
out_template <- raster("/Data/Base_Data/ALFRESCO_formatted/ALFRESCO_Master_Dataset/ALFRESCO_Model_Input_Datasets/AK_CAN_Inputs/Climate/cru_TS31/historical/pr/pr_total_mm_alf_cru_TS31_01_1901.tif")

for (i in 1:nlayers(rad.monthly)){
	print(i)
	rad.monthly.current <- subset(rad.monthly, i , drop=T)
	xyz <- cbind(coordinates(rad.monthly.current), getValues(rad.monthly.current))
	pts <- SpatialPoints(xyz[,1:2])
	pts.spdf <- SpatialPointsDataFrame(pts, as.data.frame(xyz), match.ID=TRUE,proj4string=projection(out_template))
	rasterize(pts.spdf, out_template, field=pts.spdf$V3, fun=mean,filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/outputs/girr_radiation/version2/",rad.monthly@layernames[i],".tif",sep=""), overwrite=TRUE)
}
