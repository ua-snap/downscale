# compare new / old versions of SNAP downscaled outputs 
def open_raster( fn, band=1 ):
	with rasterio.open( fn ) as rst:
		arr = rst.read(1)
	return arr

def compare_files( old, new ):
	old = open_raster(old)
	new = open_raster(new)
	diff = (old-new)
	return (np.min(diff),np.max(diff))

def wrap( x ):
	old, new = x
	return compare_files(old, new)

if __name__ == '__main__':
	import os, glob, rasterio
	import multiprocessing as mp
	import numpy as np

	variable = 'tas'
	model = '5ModelAvg'
	scenario = 'rcp85'

	old_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data_OLD_March2018/downscaled/{}/{}/{}'.format(model, scenario, variable)
	new_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/{}/{}/{}'.format(model, scenario, variable)

	old_files = glob.glob(os.path.join(old_dir,'*.tif'))
	new_files = [ fn.replace(old_dir, new_dir) for fn in old_files ]
	args = zip( old_files, new_files ) 

	pool = mp.Pool( 64 )
	out = pool.map( wrap, zip(old_files, new_files) )
	pool.close()
	pool.join()
