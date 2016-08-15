# # # # #
# this is how we will fix the bad data if it ends up coming back to us.
# # # # #
def group_val( group, arr, rst ):
	rows = np.array( [i for i,j in group] )
	cols = np.array( [j for i,j in group] )

	# remove the data outside the domain that is acceptable
	df = pd.DataFrame( {'rows':rows , 'cols':cols } )
	df = df[ (rows < rst.height) & (rows >= 0) & (cols < rst.width) & (cols >= 0)]

	# unpack it again
	rows = df.rows.tolist()
	cols = df.cols.tolist()

	vals = arr[ ( rows, cols ) ]
	vals = vals[ (vals >= 0) ] # [ hardwired ]

	if len(vals) > 0:
		new_val = vals.mean()
	else:
		new_val = 0.0 # this is hairy (49 unsolvable pixels)
	return new_val
def run_replace( fn, mask_arr ):
	from functools import partial
	rst = rasterio.open( fn )
	meta = rst.meta
	meta.update( compress='lzw' )

	arr = rst.read( 1 )
	ind = np.where( (mask_arr == 1) & (arr < 0) )
	ind_zip = zip( *ind )

	# a little neighborhood math for the queens case
	index_groups = [[(i-1,j-1),
					(i-1,j+0),
					(i-1,j+1),
					(i+0,j-1),
					(i+0,j+1),
					(i+1,j-1),
					(i+1,j+0),
					(i+1,j+1)] for i,j in ind_zip ]

	new_vals = [ group_val( group=group, arr=arr, rst=rst ) for group in index_groups ]

	rows = [ row for row, col in ind_zip ]
	cols = [ col for row, col in ind_zip ]

	arr[ ( rows, cols ) ] = new_vals
	with rasterio.open( fn.replace( '.tif', '_fix.tif' ), 'w', **meta ) as out:
		out.write_band( 1, arr )
	return fn.replace( '.tif', '_fix.tif' )

if __name__ == '__main__':
	import rasterio, glob, os
	import numpy as np
	import pandas as pd
	from pathos.mp_map import mp_map
	from functools import partial

	mask = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/extents/IEM_Mask_1km.tif'
	mask = rasterio.open( mask ).read( 1 )

	l = glob.glob( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/girr_radiation_cmip3_process/IEM/*.tif' )
	args = [ (i, mask) for i in l ]
	done = mp_map( lambda x: run_replace( *x ), args, nproc=32)




	# for fn in l:
	# 	rst = rasterio.open( fn )

	# 	meta = rst.meta
	# 	meta.update( compress='lzw' )

	# 	arr = rst.read( 1 )
	# 	ind = np.where( (mask == 1) & (arr > -3) )
	# 	ind_zip = zip( *ind )

	# 	# a little neighborhood math for the queens case
	# 	index_groups = [[(i-1,j-1),
	# 					(i-1,j+0),
	# 					(i-1,j+1),
	# 					(i+0,j-1),
	# 					(i+0,j+1),
	# 					(i+1,j-1),
	# 					(i+1,j+0),
	# 					(i+1,j+1)] for i,j in ind_zip ]

	# 	# run_partial = partial(run, arr=arr )
	# 	new_vals = map( lambda x: partial(group_val, arr=arr )( group=x ), index_groups )

	# 	arr[ ind_zip ] = new_vals
	# 	with rasterio.open( fn.replace( '.tif', '_fix.tif' ), 'w', **meta ) as out:
	# 		out.write_band( 1, arr )




# def run( group ):
# 	rows = np.array( [i for i,j in group] )
# 	cols = np.array( [j for i,j in group] )

# 	# remove the data outside the domain that is acceptable
# 	df = pd.DataFrame( {'rows':rows , 'cols':cols } )
# 	df = df[ (rows < rst.height) & (rows >= 0) & (cols < rst.width) & (cols >= 0)]

# 	# unpack it again
# 	rows = df.rows.tolist()
# 	cols = df.cols.tolist()

# 	vals = arr[ ( rows, cols ) ]
# 	vals = vals[ (vals > 0) ] # [ hardwired ]

# 	if len(vals) != 0:
# 		new_val = vals.mean()
# 	else:
# 		new_val = rst.nodata

# 	return new_val



		# new_vals = []
		# for count, group in enumerate( index_groups ):
		# 	rows = np.array( [i for i,j in group] )
		# 	cols = np.array( [j for i,j in group] )

		# 	# remove the data outside the domain that is acceptable
		# 	df = pd.DataFrame( {'rows':rows , 'cols':cols } )
		# 	df = df[ (rows < rst.height) & (rows >= 0) & (cols < rst.width) & (cols >= 0)]

		# 	# unpack it again
		# 	rows = df.rows.tolist()
		# 	cols = df.cols.tolist()

		# 	vals = arr[ ( rows, cols ) ]
		# 	vals = vals[ (vals > 0) ] # [ hardwired ]

		# 	if len(vals) != 0:
		# 		new_val = vals.mean()
		# 	else:
		# 		new_val = rst.nodata
			
		# 	new_vals.append( new_val )


# a little neighborhood math for the queens case
# ul = (-1,-1)
# uc = (-1, 0)
# ur = (-1,+1)
# l = (0,-1)
# r = (0,+1)
# ll = (+1,-1)
# lc = (+1,0)
# lr = (+1,+1)
