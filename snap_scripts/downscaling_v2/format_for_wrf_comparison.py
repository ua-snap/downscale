# convert SNAP tas / pr to WRF domain/res
def reproject( template_fn, in_fn, out_fn ):
	try:
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			_ = os.makedirs( dirname )
	except:
		pass
	
	_ = shutil.copy( template_fn, out_fn )
	with rasterio.open( out_fn, mode='r+' ) as tmp:
		arr = tmp.read(1)
		arr[:] = -9999
		tmp.write( arr, 1 )

	return subprocess.call(['gdalwarp','-q','-multi', '-r', 'cubicspline', '-srcnodata', '-9999', in_fn, out_fn ])

def run( x ):
	return reproject( *x )

if __name__ == '__main__':
	import rasterio,
	import os, glob, subprocess, shutil
	import multiprocessing as mp

	wrf_fn = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_testing_extent/wrf_extent_raw_TEST.tif'
	snap_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	out_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/wrf_20km'
	variables = [ 'tas', 'pr' ]

	# list the files from the akcan dir
	files = [ os.path.join(r, fn) for r,s,files in os.walk( snap_path ) for fn in files if fn.endswith('.tif') and fn.split('_')[0] in variables ]

	# make output_filenames and paths
	output_filenames = [ fn.replace(snap_path, out_path) for fn in files ]

	args = list( zip([wrf_fn for fn in files], files, output_filenames ) )

	pool = mp.Pool( 16 )
	out = pool.map( run, args )
	pool.close()
	pool.join()
