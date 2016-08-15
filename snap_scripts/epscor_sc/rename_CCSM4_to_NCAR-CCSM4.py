# modify the names to match what they are supposed to match
# per TK --> we should not be CCSM4 but NCAR-CCSM4
def rename_file( in_fn, out_fn, *args, **kwargs ):
	import os, shutil
	dirname = os.path.dirname( out_fn )
	
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass	
	return shutil.copy( in_fn, out_fn )

def wrap( x ):
	return rename_file( *x )

if __name__ == '__main__':
	import os, glob, subprocess
	from pathos.mp_map import mp_map

	base_paths = ['/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models', '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models']

	out = []
	for base_path in base_paths:
		for root, subs, files in os.walk( base_path ):
			out = out + [ os.path.join( root, fn ) for fn in files if 'CCSM4' in fn ]

	out_files = [ fn.replace( 'CCSM4', 'NCAR-CCSM4' ) for fn in out ]
	args = zip( out, out_files )

	_ = mp_map( wrap, args, nproc=32 )

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# # # # Here is how I renamed the 5modelAvg to 5ModelAvg, which is better
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def rename_file( in_fn, out_fn, *args, **kwargs ):
	import os, shutil
	dirname = os.path.dirname( out_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass	
	return shutil.copy( in_fn, out_fn )

def wrap( x ):
	return rename_file( *x )

if __name__ == '__main__':
	import os, glob, subprocess
	from pathos.mp_map import mp_map

	base_paths = ['/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models', '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models']

	out = []
	for base_path in base_paths:
		for root, subs, files in os.walk( base_path ):
			out = out + [ os.path.join( root, fn ) for fn in files if '5modelAvg' in fn ]

	out_files = [ fn.replace( '5modelAvg', '5ModelAvg' ) for fn in out ]
	args = zip( out, out_files )

	_ = mp_map( wrap, args, nproc=32 )

