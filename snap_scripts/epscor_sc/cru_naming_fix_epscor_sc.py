# CRU DATA NAMING HACK UPDATER

import os, glob
import pandas as pd
from pathos.mp_map import mp_map

# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/ts323'
base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_FINAL_OCT'

# change the folder names
out = []
for root, subs, files in os.walk( base_path ):
	if root.endswith( '323' ):
	# if '5M' in os.path.split( root )[-1]:
		out.append( root )

done = [ os.rename( i, os.path.join( os.path.split(i)[0],'CRU_TS323')) for i in out ]
# done = [ os.rename( i, os.path.join( os.path.split(i)[0],'5ModelAvg')) for i in out ]

# change the filenames
out = []
for root, subs, files in os.walk( base_path ):
	fn_list = [ os.path.join(root, i) for i in files if i.endswith( '.tif' ) and '323' in i ]
	# fn_list = [ os.path.join(root, i) for i in files if i.endswith( '.tif' ) and '5M' in i ]
	if fn_list > 0:
		out = out + fn_list

def rename( fn ):
	import shutil
	dirname, basename = os.path.split( fn )
	out_fn = basename.replace( 'cru_ts323', 'CRU_TS323' ).replace( 'cru_TS323', 'CRU_TS323' )
	# out_fn = basename.replace( '5MODELAVG', '5ModelAvg' )
	_ = shutil.move( fn, os.path.join( dirname, out_fn ) )
	return out_fn

done = mp_map( rename, out, nproc=32 )


