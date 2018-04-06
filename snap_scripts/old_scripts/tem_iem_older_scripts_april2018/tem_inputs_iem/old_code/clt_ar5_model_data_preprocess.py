#!/usr/bin/python2

# #
# pre-processing of raw downloaded CMIP5 data from the PCMDI portal to something that is standardized for 
# use in later downscaling to ALFRESCO AK/Canada extent and resolution needs.
# # # # # 

def group_input_filenames( prefix, root_dir ):
	import fnmatch, functools, itertools, os, glob
	import pandas as pd
	''' function that wraps some ugliness regarding returning the files we want to process '''
	def find_files( dir_path, patterns ):
		"""
		Returns a generator yielding files matching the given patterns
		:type dir_path: [str]
		:type patterns: [str]
		:rtype : [str]
		:param dir_path: Directory to search for files/directories under. Defaults to current dir.
		:param patterns: Patterns of files to search for. Defaults to ["*"]. Example: ["*.json", "*.xml"]
		"""
		import itertools, functools
		path = dir_path
		if not patterns:
			path_patterns = [ "*" ]
		else:
			path_patterns = patterns

		for root_dir, dir_names, file_names in os.walk( path ):
			filter_partial = functools.partial(fnmatch.filter, file_names)

			for file_name in itertools.chain( *map( filter_partial, path_patterns ) ):
				yield os.path.join( root_dir, file_name )
	def version_grouper( x ):
		''' groupby function for grouping by filenames '''
		dir_path = os.path.dirname( x )
		fn, _ = os.path.splitext( os.path.basename( x ) )
		# remove dates from filename -- they have a hyphen
		fn_base = '_'.join([ i for i in fn.split( '_' ) if '-' not in i ])
		# return the path element that startswith 'v' this is the version attribute
		version = [ x for x in dir_path.split( os.path.sep ) if x.startswith( 'v' ) ]
		return '_'.join([ fn_base, version[0] ])
	def drop_old_versions( df ):
		rows,cols = df.shape
		if rows > 1 & rows < 3:
			version_nums = df[ df.columns[-1] ].apply( lambda x : int( x.replace( 'v', '' ) ) )
			# max( version_nums )
			return df.drop( df[ df[ df.columns[-1] ] != 'v' + str( max( version_nums ) )].index )
		elif rows > 3:
			# potentially unnecessary
			None
		else:
			return df
	# [ !ML CHANGED! ]
	# get all matches with prefix
	matches = pd.Series([ match for match in find_files( root_dir, [ prefix ] ) ])
	input_paths = matches.apply( os.path.dirname )
	# # group by version
	# grouped = dict([ group for group in matches.groupby( matches.apply( version_grouper ) )])
	# group keys to DataFrame
	fn_list = matches.apply( lambda x: os.path.splitext( os.path.basename( x ) )[0] ).apply( lambda x: x.split( '_' )[3] )
	grouped = matches.groupby( fn_list )
	final_out = { group:files for group,files in grouped }
	# keys_df = pd.DataFrame({ key:key.split( '_' ) for key in grouped.keys() }).T
	# parse the keys / values and keep only latest versions
	# keys_df_grouped = pd.concat([ drop_old_versions(i[1]) for i in keys_df.groupby( keys_df.columns[-3] ) ])
	# # make a new dictionary holding the filenames grouped the way we want
	# final_out = { k:v for k,v in grouped.iteritems() if k in keys_df_grouped.index.tolist() }
	return final_out

def get_file_years( filename ):
	path, fn = os.path.split( filename )
	fn, ext = os.path.splitext( fn )
	split = fn.split( '_' )
	dates = split[ len( split ) - 1 ] # grab last element
	begin, end = dates.split( '-' )
	return [begin, end]

def get_modelname( filename ):
	path, fn = os.path.split( filename )
	return [ i for i in path.split( '/' ) if i in models ][0]

def concat_to_nc( filelist, output_filename, dim='time', begin_time=None, end_time=None, nc_format='NETCDF4', **kwargs ):
	'''
	take list of consecutive netcdf files (made for CMIP5 data) and stack them into a 
	single larger netcdf file.  This was necessary to overcome some bugginess in how 
	MFDataset is dealing with different calendar units on different files.  This is 
	technically valid CF-Compliant metadata, but is tricky to work with.  This hack allows
	us to get around some of this unpredictable behavior.

	PARAMETERS:
	-----------
	filelist = [list] list of string file paths to the sorted netcdf files to stack together
	output_filename = [str] path to and name of the output file to be generated (.nc extension)
	dim = [str] dimension to stack on -- default is 'time'
	begin_time = [str] PANDAS style datetime string syntax -- used in xray
	end_time = [str] PANDAS style datetime string syntax -- used in xray
	format = [str] output NetCDF format desired. valid strings are:
					'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'
					default is 'NETCDF4'
	**kwargs -- potential future arguments or overloaded args to pass through (none implemented)

	RETURNS:
	--------

	output_filename as string, with the important side-effect of writing data to disk

	'''
	import xray
	with xray.concat([ xray.open_dataset( i ).load() for i in filelist ], dim ) as ds:
		# time slicer condition
		if begin_time != None and end_time != None:
			ds = ds.loc[ { dim:slice( begin_time, end_time ) } ]
		if os.path.exists( output_filename ):
			os.remove( output_filename )
		ds.to_netcdf( output_filename, mode='w', format=nc_format )
	return output_filename

def year_greater_yearlimit_workaround( xray_dataset, desired_year_begin, desired_year_end, file_year_begin, file_year_end ):
	'''
	very specific function to deal with an issue in how PANDAS deals with datetime.
	its max datetime value in 64-bit nanoseconds from somewhere near year 1100, ends
	in roughly 2200.  This is not very ideal for working with some model outputs that
	put the data in files ranging from 2006-2300.  Soo this workaround solves the issue
	by subsetting the data using some desired years (it is assumed 12 month FULL years)
	to subset to these data ranges and return a new xray.dataset.

	PANDAS date_range functionality < 2300-07-06 00:00:00


	PARAMETERS:
	-----------
	ds = xray.Dataset object with year values outside the time limits of the package -- PANDAS
	desired_year_begin = [int] 4 digit year begin
	desired_year_end = [int] 4 digit year end
	file_year_begin = [int] 4 digit year begin in file
	file_year_end = [int] 4 digit year end in file

	RETURNS:
	--------
	new xray.Dataset object subset to the years of interest using integer indexing instead of year 
	slicing with strings.
	'''
	years = np.repeat(range( file_year_begin, file_year_end+1 ), 12 ) # 12 is for 12 months
	year_min = min( years )

	if desired_year_begin < year_min:
		begin_idx = ( year_min - desired_year_begin )
	else:
		begin_idx = np.min( np.where(years == desired_year_begin ))

	end_idx = np.max( np.where( years == desired_year_end ))
	return xray_dataset[ dict( time=range( begin_idx, end_idx + 1 ) ) ]

if __name__ == '__main__':
	import os, sys, re, glob, xray
	import numpy as np
	import pandas as pd
	import argparse

	# parse the commandline arguments
	parser = argparse.ArgumentParser( description='preprocess cmip5 input netcdf files to a common type and single files' )
	parser.add_argument( "-p", "--base_path", action='store', dest='base_path', type=str, help="path to parent directory with a subdirector(ies)y storing the data" )
	parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="string of the model name to use in the filename search i.e 'GISS-E2-R', 'IPSL-CM5A-LR',..." )
	parser.add_argument( "-v", "--variable", action='store', dest='variable', type=str, help="string of the variable name to use in the filename search i.e. 'tas', 'hur',..." )

	# parse and unpack args
	args = parser.parse_args()
	model = args.model
	variable = args.variable
	base_path = args.base_path

	# start an error log for any problem files
	output_base_path = os.path.join( base_path, 'prepped' )
	if not os.path.exists( output_base_path ):
		os.makedirs( output_base_path )

	problem_files_log = open( os.path.join( output_base_path, 'error_files.txt' ), mode='w' )

	# make some filters and filter the files we want into groups
	fn_prefix_filter = variable + '_*' + model + '*'
	file_groups = group_input_filenames( fn_prefix_filter, base_path ) # [ !ML CHANGED! ]

	for files in file_groups.values():
		try:
			files = sorted( files ) # [ !ML CHANGED! ]

			fn = files[ 0 ]
			# model = get_modelname( fn )
			output_path = os.path.join( output_base_path, model, variable )
			begin_year ='185001' # hardwired
			end_year = '210012' # hardwired

			# a handler for the historical (1850-2005) and the modeled (2006-2100) naming
			if 'historical' in os.path.basename( fn ):
				begin_year_fnout = '185001' # hardwired
				end_year_fnout = '200512' # hardwired
			else:
				begin_year_fnout = '200601' # hardwired
				end_year_fnout = '210012' # hardwired

			# this logic can be fine tuned to subset the data down to only the files we need
			# for this project it is 1850-2100.
			df = pd.DataFrame([ get_file_years(fn) for fn in files ])

			# this is the way to interrogate that dataframe for the values we want
			df = df.astype( int )
			begin_idx = (np.abs(df[0] - int( begin_year ) ) ).argmin()
			end_idx = (np.abs(df[1] - int( end_year ) ) ).argmin()

			# return the files between the desired date ranges
			if begin_idx == end_idx:
				files = [ files[ begin_idx ] ]
			else:
				files = files[ begin_idx:end_idx + 1 ]

			print files
			print '\n'

			begin_year_in = str(df.ix[ begin_idx ][0])
			end_year_in = str(df.ix[ end_idx ][1])

			# set up some vars for the output naming standardization
			cmor_table = os.path.splitext( os.path.basename( fn ) )[ 0 ].split( '_' )[ 1 ]
			experiment = scenario = os.path.splitext( os.path.basename( fn ) )[ 0 ].split( '_' )[ -2 ]
			scenario = os.path.splitext( os.path.basename( fn ) )[ 0 ].split( '_' )[ -3 ]

			if not os.path.exists( output_path ):
				os.makedirs( output_path )

			# run the concatenation and the output to a new netcdf file
			# --> and we are writing in a hack to get around the darn issue with GFDL-CM3
			#   we could just run them all with the reduce workaround, but I will keep both 
			#   in hopes that the library improves.
			if 'GFDL' in model:
				ds = reduce( lambda x,y: xray.concat( [x,y], 'time'), (xray.open_dataset( i ) for i in files) )
			else:
				ds = xray.concat([ xray.open_dataset( i ).load() for i in files ], 'time' )
			
			new_ds = year_greater_yearlimit_workaround( ds, int( begin_year_fnout[:4] ), int( end_year_fnout[:4] ), int(str(begin_year_in)[:4]), int(str(end_year_in)[:4]) )
			begin_year_fnout = str(int(begin_year_fnout[:4]) + (int(begin_year_in[:4]) - int(begin_year_fnout[:4]) )) + '01' # to update the output naming
			# output name generation
			new_fn_base = '_'.join([ variable, cmor_table, model, scenario, experiment, begin_year_fnout, end_year_fnout ]) + '.nc'
			output_filename = os.path.join( output_path, new_fn_base )
			new_ds.to_netcdf( output_filename, mode='w' )
			
			# cleanup
			ds.close()
			ds = None
			new_ds.close()
			new_ds = None
			
			# legacy version of file concatenation that is not ready for prime-time due to the quirky nature of the xray library
			# in its young state.
			# concat_to_nc( files, output_filename, dim='time', begin_time=begin_year[:4], end_time=end_year[:4] )
		except:
			print '\n--> ERROR !!!\n\n%s\n\n' % files
			problem_files_log.writelines( files )
			pass

	problem_files_log.close()


# EXAMPLE OF USE:
# some setup
# import os
# models = [ 'GISS-E2-R', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
# variables = [ 'tas', 'hur' ]
# base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs' )
# for model in models:
# 	for variable in variables:
# 		os.system( ' '.join(['python','hur_ar5_model_data_preprocess.py','-p', base_path, '-m', model, '-v', variable]) )


# # # # # # #
# --> special clt prep due to ESGF being down for the past months.
# import os
# models = [ 'GISS-E2-R', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
# variables = [ 'clt' ]
# base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5_clt_nonstandard'
# os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs/downscale_cmip5/bin' )
# for model in models:
# 	for variable in variables:
# 		os.system( ' '.join(['python','clt_ar5_model_data_preprocess.py','-p', base_path, '-m', model, '-v', variable]) )
