# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# preprocess the downloaded CMIP5 data before downscale-ing
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class Preprocess( object ):
	'''
	class for preprocessing downloaded CMIP5 data from the PCMDI.

	The idea here is that this will wrap some of the ugliness involved with how 
	some modeling groups slice up the data holdings for a variable.  We @SNAP 
	would rather work with these data in a single file prepped file when passing
	into the downscaling application.

	'''
	EXT = '.nc' # hardwired, but unlikely to change given current standard
	def __init__( self, path, variable, model, scenario, experiment, years, ext=EXT, *args, **kwargs ):
		'''
		path = [str] path cotaining potentially multiple files for a single series to bwe
		variable = [str] 
		model = [str] 
		scenario = [str] 
		experiment = [str] 
		years = [tuple] integers of beginning and ending year to preprocess (1901,2100)
		'''
		self.ext = ext
		self.path = path
		self.variable = variable
		self.model = model
		self.scenario = scenario
		self.experiment = experiment
		self.years = years
		self.filelist = self.list_files( )
		self._fileyears_dict = self._get_files_years( ) #
		self.ds = self._concat_nc_list( ) #

	def list_files( self ):
		import os, glob
		globber = os.path.join( self.path, '*'.join([ self.variable, self.model, \
										self.scenario, self.experiment, self.ext ]) )
		files = glob.glob( globber )
		return sorted( files )
	def _get_files_years( self ):
		'''
		we have to get the years from the naming convention 
		standard in each filename.  This is due to some issues in the way that
		pandas deals with times > ~2200
		'''
		import pandas as pd
		yearblob = [ fn.split( '_' )[-1].split('.')[0] for fn in self.filelist ]
		years = pd.DataFrame([ i.split('-') for i in yearblob ], columns=['begin', 'end'])	
		minyear = years.begin.min()[:4]
		minmonth = years.begin.min()[-2:]
		maxyear = years.end.max()[:4]
		maxmonth = years.end.max()[-2:]
		return { 'minmonth':minmonth, 'minyear':minyear, 'maxmonth':maxmonth, 'maxyear':maxyear }
	@staticmethod
	def _year_greater_yearlimit_workaround( xarray_dataset, desired_year_begin, \
		desired_year_end, file_year_begin, file_year_end ):
		'''
		very specific function to deal with an issue in how PANDAS deals with datetime.
		its max datetime value in 64-bit nanoseconds from somewhere near year 1100, ends
		in roughly 2200.  This is not very ideal for working with some model outputs that
		put the data in files ranging from 2006-2300.  Soo this workaround solves the issue
		by subsetting the data using some desired years (it is assumed 12 month FULL years)
		to subset to these data ranges and return a new xarray.dataset.

		PANDAS date_range functionality < 2300-07-06 00:00:00

		PARAMETERS:
		-----------
		ds = xarray.Dataset object with year values outside the time limits of the package -- PANDAS
		desired_year_begin = [int] 4 digit year begin
		desired_year_end = [int] 4 digit year end
		file_year_begin = [int] 4 digit year begin in file
		file_year_end = [int] 4 digit year end in file

		RETURNS:
		--------
		new xarray.Dataset object subset to the years of interest using integer indexing instead of year 
		slicing with strings.
		'''
		import numpy as np
		years = np.repeat(range( file_year_begin, file_year_end+1 ), 12 ) # 12 is for 12 months
		year_min = min( years )

		if desired_year_begin < year_min:
			begin_idx = ( year_min - desired_year_begin )
		else:
			begin_idx = np.min( np.where(years == desired_year_begin) )

		end_idx = np.max( np.where( years == desired_year_end ))
		return xarray_dataset[ dict( time=range( begin_idx, end_idx + 1 ) ) ]
	def _concat_nc_list( self ):
		import xarray
		import pandas as pd
		import warnings
		
		with warnings.catch_warnings():
			warnings.filterwarnings( 'error' )

			# NOTE: the below line is something that used to be needed for a misbehaving model. may need again. keep it.
			# THIS should be used, but it is being a PITA
			try:
				ds = xarray.concat([ xarray.open_dataset( i ).load() for i in self.filelist ], 'time' )
				ds = ds.sel( time=slice( str(pp.years[0]), str(pp.years[1]) ) )
			except RuntimeWarning:	
				ds = reduce( lambda x,y: xarray.concat( [x,y], 'time'), (xarray.open_dataset( i ) for i in self.filelist) )

				ds = self._year_greater_yearlimit_workaround( ds, self.years[0], self.years[1], int(self._fileyears_dict['minyear']), int(self._fileyears_dict['maxyear']) )
				if ds.time.dtype =='O':
					print( '\nWARNING: Preprocess: modified time -- bad type\n' )
					ds['time'] = pd.date_range(str(self.years[0]), str(self.years[1]+1), freq='M') # MONTHLY so leap doesnt matter and neither do days
		return ds
	def write_nc( self, output_path=None, overwrite=True, nc_format='NETCDF4' ):
		'''
		output_path = [str] path to output the newly prepped file. if None, input `path` will be used.
		overwrite = [bool] True (default) will overwrite existing outputs if exist.  Error if False.
		format = [str] output NetCDF format desired. valid strings are:
						'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'
						default is 'NETCDF4'
		'''
		import xarray, os
		begin_time = str(self.years[0])
		end_time = str(self.years[1])
		# handle output path
		if output_path == None:
			output_path = os.path.join( self.path, 'prepped' )
		if not os.path.exists( output_path ):
			os.makedirs( output_path )
		# handle overwrite and output
		output_filename = os.path.join( output_path, '_'.join([ self.variable, self.model, self.scenario, self.experiment, begin_time, end_time ]) + '.nc' )
		if os.path.exists( output_filename ) and overwrite == True:
			os.remove( output_filename )
		elif os.path.exists( output_filename ) and overwrite == False:
			raise AttributeError( 'overwrite set to False, but file exists on disk' )
		self.ds.to_netcdf( output_filename, mode='w', format=nc_format )
		return output_filename
