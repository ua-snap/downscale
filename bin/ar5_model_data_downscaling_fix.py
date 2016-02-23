# # # # #
# Tool to downscale the CMIP5 data from the PCMDI group. 
# # # # #
def shiftgrid(lon0,datain,lonsin,start=True,cyclic=360.0):
	import numpy as np
	"""
	Shift global lat/lon grid east or west.
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Arguments        Description
	==============   ====================================================
	lon0             starting longitude for shifted grid
					 (ending longitude if start=False). lon0 must be on
					 input grid (within the range of lonsin).
	datain           original data with longitude the right-most
					 dimension.
	lonsin           original longitudes.
	==============   ====================================================
	.. tabularcolumns:: |l|L|
	==============   ====================================================
	Keywords         Description
	==============   ====================================================
	start            if True, lon0 represents the starting longitude
					 of the new grid. if False, lon0 is the ending
					 longitude. Default True.
	cyclic           width of periodic domain (default 360)
	==============   ====================================================
	returns ``dataout,lonsout`` (data and longitudes on shifted grid).
	"""
	if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
		# Use all data instead of raise ValueError, 'cyclic point not included'
		start_idx = 0
	else:
		# If cyclic, remove the duplicate point
		start_idx = 1
	if lon0 < lonsin[0] or lon0 > lonsin[-1]:
		raise ValueError('lon0 outside of range of lonsin')
	i0 = np.argmin(np.fabs(lonsin-lon0))
	i0_shift = len(lonsin)-i0
	if np.ma.isMA(datain):
		dataout  = np.ma.zeros(datain.shape,datain.dtype)
	else:
		dataout  = np.zeros(datain.shape,datain.dtype)
	if np.ma.isMA(lonsin):
		lonsout = np.ma.zeros(lonsin.shape,lonsin.dtype)
	else:
		lonsout = np.zeros(lonsin.shape,lonsin.dtype)
	if start:
		lonsout[0:i0_shift] = lonsin[i0:]
	else:
		lonsout[0:i0_shift] = lonsin[i0:]-cyclic
	dataout[...,0:i0_shift] = datain[...,i0:]
	if start:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
	else:
		lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
	dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
	return dataout,lonsout
def cru_generator( n, cru_clim_list ):
	'''
	generator that will produce the cru climatologies with a
	generator and replicate for the total number of years in n
	'''
	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	for i in range( n ):
		for count, j in enumerate( cru_clim_list ):
			yield j
def standardized_fn_to_vars( fn ):
	''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
	name_convention = [ 'variable', 'cmor_table', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
	fn = os.path.basename( fn )
	fn_list = fn.split( '.' )[0].split( '_' )
	return { i:j for i,j in zip( name_convention, fn_list )}
def downscale( src, dst, cru, src_crs, src_affine, dst_crs, dst_affine, output_filename, dst_meta, variable,\
		method='cubic_spline', operation='add', output_dtype='float32', **kwargs ):
	'''
	operation can be one of two keywords for the operation to perform the delta downscaling
	- keyword strings are one of: 'add'= addition, 'mult'=multiplication, or 'div'=division (not implemented)
	- method can be one of 'cubic_spline', 'nearest', 'bilinear' and must be input as a string.
	- output_dtype can be one of 'int32', 'float32'
	'''
	from rasterio.warp import reproject, RESAMPLING
	def add( cru, anom ):
		return cru + anom
	def mult( cru, anom ):
		return cru * anom
	def div( cru, anom ):
		# return cru / anom
		# this one may not be useful, but the placeholder is here 
		return NotImplementedError

	# switch to deal with numeric output dtypes
	dtypes_switch = {'int32':np.int32, 'float32':np.float32}

	# switch to deal with different resampling types
	method_switch = { 'nearest':RESAMPLING.nearest, 'bilinear':RESAMPLING.bilinear, 'cubic_spline':RESAMPLING.cubic_spline }
	method = method_switch[ method ]

	# reproject src to dst
	out = np.zeros( dst.shape ) 
	reproject( src,
				out,
				src_transform=src_affine,
				src_crs=src_crs,
				dst_transform=dst_affine,
				dst_crs=dst_crs,
				resampling=method )
	# switch to deal with different downscaling operators
	operation_switch = { 'add':add, 'mult':mult, 'div':div }
	downscaled = operation_switch[ operation ]( cru, out )

	# reset any > 100 values to 95 if the variable is cld or hur
	if variable == 'clt' or variable == 'hur' or variable == 'cld':
		downscaled[ downscaled > 100.0 ] = 95.0

	# give the proper fill values to the oob regions
	downscaled.fill_value = dst_meta['nodata']
	downscaled = downscaled.filled()

	# this is a geotiff creator so lets pass in the lzw compression
	dst_meta.update( compress='lzw' )
	with rasterio.open( output_filename, 'w', **dst_meta ) as out:
		out.write( downscaled.astype( dtypes_switch[ output_dtype ] ), 1 )
	return output_filename
def run( args ):
	''' 
	simple function wrapper for unpacking an argument dict 
	to the downscale function for getting around the single 
	argument pass to multiprocessing.map implementation issue.
	'''
	return( downscale( **args ) )

if __name__ == '__main__':
	import pandas as pd
	import numpy as np
	import os, sys, re, xray, rasterio, glob, argparse
	from rasterio import Affine as A
	from rasterio.warp import reproject, RESAMPLING
	from pathos import multiprocessing as mp

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='preprocess cmip5 input netcdf files to a common type and single files' )
	parser.add_argument( "-mi", "--modeled_fn", nargs='?', const=None, action='store', dest='modeled_fn', type=str, help="path to modeled input filename (NetCDF); default:None" )
	parser.add_argument( "-hi", "--historical_fn", nargs='?', const=None, action='store', dest='historical_fn', type=str, help="path to historical input filename (NetCDF); default:None" )
	parser.add_argument( "-o", "--output_path", action='store', dest='output_path', type=str, help="string path to the output folder containing the new downscaled outputs" )
	parser.add_argument( "-cbt", "--climatology_begin_time", nargs='?', const='196101', action='store', dest='climatology_begin', type=str, help="string in format YYYYMM or YYYY of the beginning month and potentially (year) of the climatology period" )
	parser.add_argument( "-cet", "--climatology_end_time", nargs='?', const='199012', action='store', dest='climatology_end', type=str, help="string in format YYYYMM or YYYY of the ending month and potentially (year) of the climatology period" )
	parser.add_argument( "-plev", "--plev", nargs='?', const=None, action='store', dest='plev', type=int, help="integer value (in millibars) of the desired pressure level to extract, if there is one." )
	parser.add_argument( "-cru", "--cru_path", action='store', dest='cru_path', type=str, help="path to the directory storing the cru climatology data derived from CL2.0" )
	parser.add_argument( "-at", "--anomalies_calc_type", nargs='?', const='absolute', action='store', dest='anomalies_calc_type', type=str, help="string of 'proportional' or 'absolute' to inform of anomalies calculation type to perform." )
	parser.add_argument( "-m", "--metric", nargs='?', const='metric', action='store', dest='metric', type=str, help="string of whatever the metric type is of the outputs to put in the filename." )
	parser.add_argument( "-dso", "--downscale_operation", action='store', dest='downscale_operation', type=str, help="string of 'add', 'mult', 'div', which refers to the type or downscaling operation to use." )
	parser.add_argument( "-nc", "--ncores", nargs='?', const=2, action='store', dest='ncores', type=int, help="integer valueof number of cores to use. default:2" )

	# parse args
	args = parser.parse_args()

	# unpack args
	modeled_fn = args.modeled_fn
	historical_fn = args.historical_fn

	# temporary
	# modeled_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/MRI-CGCM3/clt/clt_Amon_MRI-CGCM3_rcp26_r1i1p1_200601_210012.nc'
	# historical_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/MRI-CGCM3/clt/clt_Amon_MRI-CGCM3_rcp85_r1i1p1_200601_210012.nc'

	output_path = args.output_path
	climatology_begin = args.climatology_begin
	climatology_end = args.climatology_end
	plev = args.plev
	cru_path = args.cru_path
	anomalies_calc_type = args.anomalies_calc_type
	metric = args.metric
	downscale_operation = args.downscale_operation
	ncores = args.ncores

	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# [NOTE]: hardwired raster metadata meeting the ALFRESCO Model's needs for 
	# perfectly aligned inputs this is used as template metadata that 
	# is used in output generation. template raster filename below:
	# '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/
	#	TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	meta_3338 = {'affine': A(2000.0, 0.0, -2173223.206087799, 
					0.0, -2000.0, 2548412.932644147),
				'count': 1,
				'crs': {'init':'epsg:3338'},
				'driver': u'GTiff',
				'dtype': 'float32',
				'height': 1186,
				'nodata': -3.4e+38,
				'width': 3218,
				'compress':'lzw'}

	# output template numpy array same dimensions as the template
	dst = np.empty( (1186, 3218) )
	
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# condition to deal with reading in historical data if needed.
	if modeled_fn is not None and historical_fn is not None:
		print 'here'
		# parse the input name for some file metadata
		output_naming_dict = standardized_fn_to_vars( modeled_fn )

		# this is to maintain cleanliness
		variable = output_naming_dict[ 'variable' ]

		# read in both modeled and historical
		ds = xray.open_dataset( modeled_fn )
		ds = ds[ variable ].load()
		clim_ds = xray.open_dataset( historical_fn )
		clim_ds = clim_ds[ variable ].load()
		# generate climatology / anomalies
		clim_ds = clim_ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
		climatology = clim_ds.groupby( 'time.month' ).mean( 'time' )

		# find the begin/end years of the prepped files
		dates = ds.time.to_pandas()
		years = dates.apply( lambda x: x.year )
		begin_time = years.min()
		end_time = years.max()
		print 'here'
		del clim_ds
	elif historical_fn is not None and modeled_fn is None:
		# parse the input name for some file metadata
		output_naming_dict = standardized_fn_to_vars( historical_fn )
		
		# this is to maintain cleanliness
		variable = output_naming_dict[ 'variable' ]

		# read in historical
		ds = xray.open_dataset( historical_fn )
		ds = ds[ variable ].load()
		# generate climatology / anomalies
		climatology = ds.loc[ {'time':slice(climatology_begin,climatology_end)} ]
		climatology = climatology.groupby( 'time.month' ).mean( 'time' )

		# find the begin/end years of the prepped files
		dates = ds.time.to_pandas()
		years = dates.apply( lambda x: x.year )
		begin_time = years.min()
		end_time = years.max()

	else:
		NameError( 'ERROR: must have both modeled_fn and historical_fn, or just historical_fn' )

	# standardize the output pathing
	if output_naming_dict[ 'variable' ] == 'clt':
		variable_out = 'cld'
		print 'here'
	else:
		variable_out = output_naming_dict[ 'variable' ]

	output_path = os.path.join( output_path, 'ar5', output_naming_dict['model'], variable_out, 'downscaled' )
	if not os.path.exists( output_path ):
		os.makedirs( output_path )


	# if there is a pressure level to extract, extract it
	if plev is not None:
		plevel, = np.where( ds.plev == plev )
		ds = ds[ :, plevel[0], ... ]
		climatology = climatology[ :, plevel[0], ... ]

	# deal with different anomaly calculation types
	if anomalies_calc_type == 'absolute':
		anomalies = ds.groupby( 'time.month' ) - climatology
	elif anomalies_calc_type == 'proportional':
		print 'here'
		anomalies = ds.groupby( 'time.month' ) / climatology
	else:
		NameError( 'anomalies_calc_type can only be one of "absolute" or "proportional"' )
	# some setup of the output raster metadata
	time_len, rows, cols = anomalies.shape
	crs = 'epsg:4326'
	affine = A( *[np.diff( ds.lon )[ 0 ], 0.0, -180.0, 0.0, -np.diff( ds.lat )[ 0 ], 90.0] )
	count = time_len
	resolution = ( np.diff( ds.lat )[ 0 ], np.diff( ds.lon )[ 0 ] )

	# close the dataset and clean it up
	ds = None

	# shift the grid to Greenwich Centering
	dat, lons = shiftgrid( 180., anomalies[:], anomalies.lon.data, start=False )

	# metadata for input?
	meta_4326 = {'affine':affine,
				'height':rows,
				'width':cols,
				'crs':crs,
				'driver':'GTiff',
				'dtype':np.float32,
				'count':time_len,
				'compress':'lzw' }
	# build some filenames for the outputs to be generated
	months = [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12' ]
	years = [ str(year) for year in range( begin_time, end_time + 1, 1 ) ]
	# combine the months and the years
	combinations = [ (month, year) for year in years for month in months ]

	output_filenames = [ os.path.join( output_path, '_'.join([variable_out, 'metric', output_naming_dict['model'], output_naming_dict['scenario'], output_naming_dict['experiment'], month, year]) + '.tif' ) for month, year in combinations ]
	print 'here'
	# load the baseline CRU CL2.0 data 
	# [NOTE]: THIS ASSUMES THEY ARE THE ONLY FILES IN THE DIRECTORY -- COULD BE A GOTCHA
	cru_files = glob.glob( os.path.join( cru_path, '*.tif' ) )
	cru_files.sort()
	cru_stack = [ rasterio.open( fn ).read( 1 ) for fn in cru_files ]
	# this is a hack to make a masked array with the cru data
	cru_stack = [ np.ma.masked_where( cru == cru.min(), cru ) for cru in cru_stack ]
	cru_gen = cru_generator( len(output_filenames), cru_stack )
	print 'here'
	
	# cleanup some uneeded vars that are hogging RAM
	del climatology, anomalies

	# run in parallel using PATHOS
	pool = mp.Pool( processes=ncores )
	args_list = zip( np.vsplit( dat, time_len ), output_filenames, cru_gen )
	del dat, cru_gen, cru_stack

	out = pool.map( run, [{'src':src, 'output_filename':fn, 'dst':dst, 'cru':cru, 'src_crs':meta_4326[ 'crs' ], 'src_affine':meta_4326[ 'affine' ], \
							'dst_crs':meta_3338[ 'crs' ], 'dst_affine':meta_3338[ 'affine' ], 'dst_meta':meta_3338, 'operation':downscale_operation, 'variable':variable } \
							for src,fn,cru in args_list ] )
	pool.close()


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
# # # # # # # # 	 SOME EXAMPLES OF USE 		# # # # # # # # # # # # # #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# # TO RUN THE CLOUDS DOWNSCALING USE THIS EXAMPLE:
import os
import pandas as pd
import numpy as np

# change to the script repo
os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )


variable = 'clt' # AR5 naming convention cloud fraction

# to run the futures:
prepped_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped'
file_groups = [ [os.path.join(root,f) for f in files] for root, sub, files in os.walk( prepped_dir ) if len(files) > 0 and files[0].endswith('.nc') and variable in files[0] ]
variable = 'cld' # swap it back for the cru naming convention

def make_rcp_file_pairs( file_group ):
	# there is only one historical per group since these have been pre-processed to a single file and date range
	historical = [ file_group.pop( count ) for count, i in enumerate( file_group ) if 'historical' in i ]
	return zip( np.repeat( historical, len(file_group) ).tolist(), file_group )

grouped_pairs = [ make_rcp_file_pairs( file_group ) for file_group in file_groups ]
for file_group in grouped_pairs:
	for historical_fn, modeled_fn in file_group:
		if 'MRI-CGCM3_rcp26_' in modeled_fn:
			print( 'running: %s' % os.path.basename( modeled_fn ) )
			output_path = '/Data/malindgren/cru_november_final'
			climatology_begin = '1961-01'
			climatology_end = '1990-12'
			cru_path = os.path.join( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final/cru_cl20', variable, 'akcan' )
			anomalies_calc_type = 'proportional'
			metric = 'pct'
			downscale_operation = 'mult'
			ncores = '30'
			# future modeled data
			# build the args
			modeled_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/MRI-CGCM3/clt/clt_Amon_MRI-CGCM3_rcp26_r1i1p1_200601_210012.nc'
			historical_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped/MRI-CGCM3/clt/clt_Amon_MRI-CGCM3_historical_r1i1p1_185001_200512.nc'
			args_tuples = [ ( 'mi', modeled_fn ),
							( 'hi', historical_fn ),
							( 'o', output_path ),
							( 'cbt', climatology_begin ),
							( 'cet', climatology_end ),
							( 'cru', cru_path ),
							( 'at', anomalies_calc_type ),
							( 'm', metric ),
							( 'dso', downscale_operation ),
							( 'nc', ncores ) ]

			args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
			os.system( 'python clt_ar5_model_data_downscaling.py ' + args )

			del modeled_fn

	output_path = output_path
	climatology_begin = climatology_begin
	climatology_end = climatology_end
	plev = plev
	cru_path = cru_path
	anomalies_calc_type = anomalies_calc_type
	metric = metric
	downscale_operation = downscale_operation
	ncores = ncores

# 	# now historical modeled data
# 	# build the args
# 	print( 'running: %s' % os.path.basename( historical_fn ) )

# 	args_tuples = [	( 'hi', historical_fn ),
# 					( 'o', output_path ),
# 					( 'cbt', climatology_begin ),
# 					( 'cet', climatology_end ),
# 					( 'cru', cru_path ),
# 					( 'at', anomalies_calc_type ),
# 					( 'm', metric ),
# 					( 'dso', downscale_operation ),
# 					( 'nc', ncores ) ]

# 	args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])

# 	os.system( 'ipython -- ar5_model_data_downscaling.py ' + args )

# # # # # # # # # # # # # # # # # # # # 
# # TO RUN THE TEMPERATURE DOWNSCALING USE THIS EXAMPLE:
# import os
# import pandas as pd
# import numpy as np

# # change to the script repo
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )

# variable = 'tas'

# # to run the futures:
# prepped_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped'
# file_groups = [ [os.path.join(root,f) for f in files] for root, sub, files in os.walk( prepped_dir ) if len(files) > 0 and files[0].endswith('.nc') and variable in files[0] ]

# def make_rcp_file_pairs( file_group ):
# 	# there is only one historical per group since these have been pre-processed to a single file and date range
# 	historical = [ file_group.pop( count ) for count, i in enumerate( file_group ) if 'historical' in i ]
# 	return zip( np.repeat( historical, len(file_group) ).tolist(), file_group )

# grouped_pairs = [ make_rcp_file_pairs( file_group ) for file_group in file_groups ]
# for file_group in grouped_pairs:
# 	for historical_fn, modeled_fn in file_group:
# 		print( 'running: %s' % os.path.basename( modeled_fn ) )
# 		output_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final'
# 		climatology_begin = '1961-01'
# 		climatology_end = '1990-12'
# 		cru_path = os.path.join( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final/cru_cl20',variable,'akcan' )
# 		anomalies_calc_type = 'absolute'
# 		metric = 'C'
# 		downscale_operation = 'add'
# 		ncores = '30'
# 		# future modeled data
# 		# # build the args
# 		args_tuples = [ ( 'mi', modeled_fn ),
# 						( 'hi', historical_fn ),
# 						( 'o', output_path ),
# 						( 'cbt', climatology_begin ),
# 						( 'cet', climatology_end ),
# 						( 'cru', cru_path ),
# 						( 'at', anomalies_calc_type ),
# 						( 'm', metric ),
# 						( 'dso', downscale_operation ),
# 						( 'nc', ncores ) ]

# 		args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])

# 		os.system( 'python clt_ar5_model_data_downscaling.py ' + args )

#		del modeled_fn

# 	# now historical modeled data
# 	# # build the args by pop(-ping) out the first entry which is modeled_fn
# 	args_tuples.pop(0)
# 	args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])

# 	os.system( 'ipython -- ar5_model_data_downscaling.py ' + args )

# # # # # # # # # # # # # # # # # # # # 
# # TO RUN THE RELATIVE HUMIDITY DOWNSCALING USE THIS EXAMPLE:
# import os
# import pandas as pd
# import numpy as np

# # change to the script repo
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )

# # variable we are running
# variable = 'hur'

# # to run the futures:
# prepped_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped'
# file_groups = [ [os.path.join(root,f) for f in files] for root, sub, files in os.walk( prepped_dir ) if len(files) > 0 and files[0].endswith('.nc') and variable in files[0] ]

# def make_rcp_file_pairs( file_group ):
# 	# there is only one historical per group since these have been pre-processed to a single file and date range
# 	historical = [ file_group.pop( count ) for count, i in enumerate( file_group ) if 'historical' in i ]
# 	return zip( np.repeat( historical, len(file_group) ).tolist(), file_group )

# grouped_pairs = [ make_rcp_file_pairs( file_group ) for file_group in file_groups ]
# for file_group in grouped_pairs:
# 	for historical_fn, modeled_fn in file_group:
# 		print( 'running: %s' % os.path.basename( modeled_fn ) )
# 		output_path = '/Data/malindgren/cru_november_final'
# 		plev = '1000'
# 		climatology_begin = '1961-01'
# 		climatology_end = '1990-12'
# 		cru_path = os.path.join( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_november_final/cru_cl20',variable,'akcan' )
# 		anomalies_calc_type = 'proportional'
# 		metric = 'pct'
# 		downscale_operation = 'mult'
# 		ncores = '32'
# 		# future modeled data

# 		args_tuples = [ ( 'mi', modeled_fn ),
# 						( 'hi', historical_fn ),
# 						( 'o', output_path ),
# 						( 'plev', plev ),
# 						( 'cbt', climatology_begin ),
# 						( 'cet', climatology_end ),
# 						( 'cru', cru_path ),
# 						( 'at', anomalies_calc_type ),
# 						( 'm', metric ),
# 						( 'dso', downscale_operation ),
# 						( 'nc', ncores ) ]

# 		args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])
# 		os.system( 'python ar5_model_data_downscaling.py ' + args )

# 		del modeled_fn

# 	# now historical modeled data
# 	# # build the args by pop(-ping) out the first entry which is modeled_fn
# 	args_tuples.pop(0)
# 	args = ''.join([ ' -'+flag+' '+value for flag, value in args_tuples ])

# 	os.system( 'ipython -- ar5_model_data_downscaling.py ' + args )
