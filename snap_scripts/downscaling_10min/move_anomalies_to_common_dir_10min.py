# move the anomalies (deltas) that have been resampled to 2km
# to a new directory
def move( fn, new_fn, *args ):
	import shutil, os

	dirname, basename = os.path.split( new_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass
		
	return os.rename( fn, new_fn )

def make_args( fn ):
	dirname, basename = os.path.split( fn )
	new_dirname = dirname.replace( '/anom', '' ).replace( 'downscaled_10min', 'anomalies_10min' )
	new_fn = os.path.join( new_dirname, basename )
	return fn, new_fn

if __name__ == '__main__':
	import os
	import numpy as np
	from pathos.mp_map import mp_map

	# base path 
	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_10min'

	# list the dirs we want to move
	args = [ make_args( os.path.join( root, fn ) ) for root,s,files in os.walk( base_path ) if 'anom' in root for fn in files ]

	# move 'em
	done = mp_map( lambda x: move( *x ), args, nproc=64 )

	# remove the old 'anom' dirs
	dirs = np.unique( np.array([ root for root,s,files in os.walk( base_path ) if 'anom' in root ]) ).tolist()
	_ = [ os.system( 'rm -r {}'.format( d ) ) for d in dirs ]