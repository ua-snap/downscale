# CRU DATA NAMING HACK UPDATER

import os, glob
import pandas as pd
from pathos.mp_map import mp_map

base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/ts323'

out = []
for root, subs, files in os.walk( base_path ):
	fn_list = [ os.path.join(root, i) for i in files if i.endswith( '.tif' ) ]
	if fn_list > 0:
		out = out + fn_list

def rename( fn ):
	import shutil
	dirname, basename = os.path.split( fn )
	out_fn = basename.replace( 'cru_ts323', 'CRU_TS323' ).replace( 'cru_TS323', 'CRU_TS323' )
	_ = shutil.move( fn, os.path.join( dirname, out_fn ) )
	return out_fn

done = mp_map( rename, out, nproc=32 )


