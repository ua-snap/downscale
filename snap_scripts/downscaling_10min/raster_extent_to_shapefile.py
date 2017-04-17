# make extent shapefile from raster
# use the auto cli-generator fire.
import fire 

class ExtentToShapefile( object ):
	@staticmethod
	def bounds_to_extent( bounds ):
		'''
		take input list of bounds and return an extent list of 
		coordinate pairs. Can be used to generate a shapely Polygon.
		'''
		l,b,r,t = bounds
		return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]
	@staticmethod
	def extent_to_polygon( extent ):
		from shapely.geometry import Polygon
		return Polygon( extent )
	def extent_to_shapefile( self, fn, output_filename ):
		''' 
		convert rasterio-readable raster extent to a Polygon
		shapefile.

		ARGUMENTS:
		----------
		fn = raster path
		output_filename = shapefile path to create

		RETURNS:
		--------
		path to new shapefile
		'''	
		import rasterio, os
		import geopandas as gpd

		# get poly
		rst = rasterio.open( fn )
		ext = self.bounds_to_extent( rst.bounds )
		pol = self.extent_to_polygon( ext )

		# make shapefile
		df = gpd.GeoDataFrame.from_dict( { 1:{ 'id':1, 'geometry':pol } }, orient='index' )
		gdf = gpd.GeoDataFrame( df, crs=rst.crs, geometry='geometry' )

		# cleanup
		if os.path.exists( output_filename ):
			os.unlink( output_filename )

		gdf.to_file( output_filename )
		return output_filename

if __name__ == '__main__':
	fire.Fire( ExtentToShapefile )


# # # # CLI EXAMPLE -- FIRE is COOOL! # # # # #
# python raster_extent_to_shapefile.py extent_to_shapefile 
# 		--fn /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/sunp_cru_cl20_akcan_01_1961-1990_GCLL_trim.tif 
# 		--output_filename /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/sunp_cru_cl20_akcan_01_1961-1990_GCLL_trim.shp
#
# or:
# python raster_extent_to_shapefile.py extent_to_shapefile /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/sunp_cru_cl20_akcan_01_1961-1990_GCLL_trim.tif /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/akcan_10min_template/sunp_cru_cl20_akcan_01_1961-1990_GCLL_trim.shp