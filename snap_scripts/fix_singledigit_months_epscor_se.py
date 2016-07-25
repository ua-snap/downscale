# walk directories and repair single digit months to 2 digit

def get_month( fn ):
	return fn.split( '.' )[ 0 ].split( '_' )[ -2 ]
def get_year( fn ):
	return fn.split( '.' )[ 0 ].split( '_' )[ -1 ]

def fix_month( fn ):
	lookup = {'1':'01','2':'02','3':'03','4':'04',\
				'5':'05','6':'06','7':'07','8':'08','9':'09',\
				'10':'10','11':'11','12':'12'}
	month = get_month( fn )
	year = get_year( fn )

	if month in lookup.keys():
		new_month = lookup[ month ]
		new_fn = fn.replace( '_'+month+'_', '_'+new_month+'_' )
		os.rename( fn, new_fn )
	return 1

if __name__ == '__main__':
	import os, glob
	from pathos import multiprocessing as mp

	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2'
	
	out_files = []
	for root, subs, files in os.walk( base_dir ):
		if len([ i for i in files if i.endswith( '.tif' ) ]) > 0:
			out_files = out_files + [ os.path.join( root, j ) for j in files ]


	pool = mp.Pool( 32 )
	pool.map( fix_month, out_files )
	pool.close()
	pool.join()
	pool.terminate()

