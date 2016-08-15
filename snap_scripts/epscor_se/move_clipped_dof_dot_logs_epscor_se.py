import os, shutil
from pathos.mp_map import mp_map

in_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_dof_dot_logs'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_AUG2016/derived/grids'
out_dir_lookup = {'dot':'decadal_dot','dof':'decadal_dof','logs':'decadal_logs' }

# get all files 
file_list = []
for root, subs, files in os.walk( in_path ):
	if len(files) > 0:
		file_list = file_list + [ os.path.join( root, fn ) for fn in files ]

# make some arguments for passing to shutil.copy
args = []
for fn in file_list:
	folder_name = out_dir_lookup[ os.path.basename( fn ).split( '_' )[0] ]
	out_path = os.path.join( output_path, folder_name )	
	new_fn = fn.replace( in_path, out_path ) 
	args = args + [(fn, new_fn)]

def copy_it( x ):
	fn, new_fn = x
	out_path = os.path.dirname( new_fn )
	try:
		if not os.path.exists( out_path ):
			os.makedirs( out_path )
	except:
		pass
	# out_path.replace( '/dof', '' )
	return shutil.copy( fn, new_fn )

_ = mp_map( copy_it, args, nproc=16 )


