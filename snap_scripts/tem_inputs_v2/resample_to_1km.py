# resample to 1km
def run( fn, out_fn ):
	import os
	dirname = os.path.dirname( out_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	os.system( 'gdalwarp -q -r bilinear -overwrite -t_srs EPSG:3338 -tr 1000 1000 -srcnodata -9999 -dstnodata -9999 -co COMPRESS=LZW {} {}'.format( fn, out_fn ) )
	return out_fn

def wrap( x ):
	return run( *x )

if __name__ == '__main__':
	import os
	from pathos.mp_map import mp_map

	# directories
	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cru_40/downscaled'
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_1km'

	# list the data
	files = [os.path.join(root, fn) for root, subs, files in os.walk( base_path ) for fn in files if fn.endswith( '.tif' ) and 'anom' not in fn ]
	out_files = [ fn.replace( base_path, output_path ) for fn in files ]

	# run it in parallel
	done = mp_map( wrap, zip( files, out_files ), nproc=32 )
