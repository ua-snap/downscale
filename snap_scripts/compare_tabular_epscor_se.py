# COMPARE TAS AND TASMIN / TASMAX

if __name__ == '__main__':
	import os, glob
	import pandas as pd

	# inputs
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_DATA/derived_outputs/csv'

	l = pd.Series( glob.glob( os.path.join( base_path, '*.csv' ) ) )
	l2 = pd.Series( [ i for i in l if 'tasmin_' in i or 'tas_' in i ] )

	args = zip( *[ sorted(j.tolist()) for i, j in l2.groupby( l2.apply(lambda x: os.path.basename(x).split('_')[0]) ) ])

def get_variable( x ):
	return os.path.basename(x).split('_')[0]

def compare_by_decade( x ):
	a,b = x
	var_a = get_variable(a)
	var_b = get_variable(b)
	a = pd.read_csv( a, index_col=0 )
	b = pd.read_csv( b, index_col=0 )

	a.index = pd.Series(a.index).apply( lambda x: x.replace(var_a+'_', '') )
	b.index = pd.Series( b.index).apply( lambda x: x.replace(var_b+'_', '') )

	return b - a

compare_by_decade( args[0] )


