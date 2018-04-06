# # # # EPSCoR SE # # # 
# this is how I changed the filenames after downscaling to the standard names
# used in SNAP downscaling from the ones stored in the CRU_TS323 data series
# # # 
def rename_file( in_fn, out_fn, *args, **kwargs ):
	import os
	return os.rename( in_fn, out_fn )

def wrap( x ):
	return rename_file( *x )

	import os, glob, shutil, subprocess
	from pathos.mp_map import mp_map

	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2/CRU_TS323/ts323/historical'

	old_vars = ['tmp', 'tmn', 'tmx']
	new_vars = ['tas', 'tasmin', 'tasmax']

	for old_var, new_var in zip(old_vars, new_vars)[:1]:
		out = []
		for root, subs, files in os.walk( base_path ):
			out = out + [ os.path.join( root, fn ) for fn in files if old_var+'_' in fn ]

		out_files = [ fn.replace( old_var+'_', new_var+'_' ) for fn in out ]
		args = zip( out, out_files )

		_ = mp_map( wrap, args, nproc=32 )
