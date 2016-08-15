# CRU DATA NAMING HACK UPDATER

import os, glob
import pandas as pd

base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru'

out = []
for root, subs, files in os.walk( base_path ):
	fn_list = [ os.path.join(root, i) for i in files if i.endswith( '.tif' ) ]
	if fn_list > 0:
		out = out + fn_list

out = pd.Series( out )
grouper = [ os.path.basename(i).split('_')[0] for i in out ]
grouped = out.groupby( grouper )

grouped = [ i for i in grouped ]

for group, files in grouped:
	if group in ['tmx','tmn','tmp','pre']:
		print( group )
		name_lookup = {'tmx':'tasmax', 'tmn':'tasmin', 'tmp':'tas', 'pre':'pr' }
		new_var = name_lookup[ group ]
		done = [ os.rename(i, i.replace( group, new_var )) for i in files ]
	

