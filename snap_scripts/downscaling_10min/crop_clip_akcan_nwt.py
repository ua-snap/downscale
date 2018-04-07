
def crop_clip( shp_fn, rst_fn, out_fn ):
	'''
	crop/clip to the shapefile we want using gdalwarp.
	'''
	import subprocess, os
	try:
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	# proj4string = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	proj4string = 'EPSG:3338'
	subprocess.call(['gdalwarp', '-q', '-overwrite', '-tap','-t_srs', proj4string, '-co', 'COMPRESS=LZW', '-tr', '15000', '15000', 
						'-srcnodata', '-9999', '-dstnodata', '-9999', '-wo', 'CUTLINE_ALL_TOUCHED=TRUE', '-crop_to_cutline', '-cutline',
						shp_fn, rst_fn, out_fn ])
		subprocess.call(['gdalwarp', '-q', '-overwrite', '-co', 'COMPRESS=LZW', '-tr', '15000', '15000', 
						'-srcnodata', '-9999', '-dstnodata', '-9999', '-wo', 'CUTLINE_ALL_TOUCHED=TRUE', '-crop_to_cutline', '-cutline',
						shp_fn, rst_fn, out_fn ])
	return out_fn

def wrap( x ):
	''' wrapper for clean multiprocessing call to pool.map '''
	return crop_clip( *x )

if __name__ == '__main__':
	import os
	from pathos.mp_map import mp_map

	# setup args
	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_10min_nwt'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_10min_nwt_clip'
	ncpus = 16
	subdomain_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/masks/mask_10min_nwt_clipper.shp'

	# BUILD ARGS
	args_list = [ ( subdomain_fn, os.path.join( root, fn ), os.path.join( root, fn ).replace( base_path, output_path ) ) for root, subs, files in os.walk( base_path ) for fn in files if fn.endswith( '.tif' ) ]

	# RUN IN PARALLEL
	out = mp_map( wrap, args_list, nproc=ncpus )

# # # SOME NOTES ON THE CLIPPING OF THE 10' outputs on Feb 8, 2018 ...
# running this on Phobos with the new CUTLINE_ALL_TOUCHED=TRUE warp option as a test to see if we can not clip the coastline unnecessarily.
# I am not sure if this is even best practice, but I fear that by not doing it we might introduce some stuff we dont want.  The way I am 
# thinking about it now is that it is easier to re-clip than it is to bring back sometihng that is lost in a previous clipping strategy.

