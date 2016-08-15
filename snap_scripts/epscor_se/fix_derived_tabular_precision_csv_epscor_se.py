import os, glob
from copy import deepcopy
import pandas as pd
import numpy as np

path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_tabular'
output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_tabular_rounded'

if not os.path.exists( output_path ):
	os.makedirs( output_path )

for v in ['pr', 'tas']:
	for fn in glob.glob( os.path.join( path, v + '*.csv' ) ):
		output_filename = os.path.join( output_path, os.path.basename( fn ) )
		df = pd.read_csv( fn, index_col=0 )
		if v == 'pr':
			# round to an int then convert to integer
			df_rounded = deepcopy( df ).apply( lambda x: np.rint( x ) ).astype( np.int ).astype( str ).to_csv( output_filename, sep=',' )
		else:
			# round it and write it to disk
			df_rounded = deepcopy( df ).round( decimals=1 )
			df_rounded.apply( lambda x: x.apply( lambda x: '%2.1f' % x) ).to_csv( output_filename, sep=',', float_format='%11.6f')