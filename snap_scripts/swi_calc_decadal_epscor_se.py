# # # # PYTHON BELOW

def list_files( base_dir, wildcard='tas_*.tif' ):
	import glob
	out = []
	for root, subs, files in os.walk( base_dir ):
		if len( [i for i in files if 'tas_' in i ] ) > 0:
			out = out + glob.glob(os.path.join( root, wildcard ) )
	return out

def make_grouper( x ):
	group = os.path.basename( x ).split('.')[0].split('_')
	group[7] = group[7][:3]+'0s'
	return '_'.join([ group[i] for i in [4,5,7] ])

def swi( files, output_path ):
	import rasterio
	import numpy as np

	months = files.apply( lambda x: x.split('_')[-2] )
	files = [ fn for fn, month in zip( files, months ) if month in ['05', '06', '07', '08', '09'] ]

	arr = np.array( [ rasterio.open(i).read(1) for i in files ] )
	arr[ arr < 0 ] = 0
	arr = np.sum( arr, axis=0 )
	
	variable, metric, units, project, model, \
		scenario, month, year = \
			os.path.basename( files[0] ).split( '.' )[0].split( '_' )

	# make sure the year is a decadal year for grouping if annual...
	# this works with both decadal and annual monthlies.
	year = year[:3]+'0s'
	rst = rasterio.open( files[0] )
	mask = rst.read_masks( 1 )
	meta = rst.meta

	_ = meta.pop( 'transform' )
	meta.update( compress='lzw' )
	output_filename = os.path.join( output_path, model, scenario, 'swi', '_'.join(['swi', 'cumulative', units, project, model, scenario, year])+'.tif' )

	try:
		dirname = os.path.dirname( output_filename )
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass

	with rasterio.open( output_filename, 'w', **meta ) as out:
		arr[ mask == 0 ] = meta[ 'nodata' ]
		out.write( arr, 1 )
	return output_filename

if __name__ == '__main__':
	import os
	import pandas as pd

	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_v2/monthly_decadals'
	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_v2/swi_decadals'
	l = pd.Series( list_files( base_dir ) )

	groups = l.apply( make_grouper )
	grouped = l.groupby( groups )
	done = grouped.apply( swi, output_path=output_path )

