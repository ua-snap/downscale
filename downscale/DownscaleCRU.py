# # #
# Downscale CRU Historical TS3.x data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #

class DownscaleCRU( object ):
	'''
	methods to downscale the Climatic Research Unit's (CRU) Historical 
	Time Series data using a 12-month climatology pre-processed to the final
	output domain and resolution.  Typically we use a PRISM climatology or a 
	CRU CL2.0 climatology for these purposes.

	'''
	def __init__( self, cru_ts, clim_path, template_raster_fn, base_path, climatology_begin='1961', climatology_end='1990', ncores=2, \
		absolute=True, metric='metric', variable=None, post_downscale_function=None, src_crs={'init':'epsg:4326'}, write_anomalies=True, *args, **kwargs ):
		self.cru_ts = cru_ts
		self.clim_path = clim_path
		self.template_raster_fn = template_raster_fn
		self.base_path = base_path
		self.climatology_begin = climatology_begin
		self.climatology_end = climatology_end
		self.ncores = ncores
		self.absolute = absolute
		self.metric = metric
		self.variable = variable
		self.post_downscale_function = post_downscale_function
		self.src_crs = src_crs
		self.utils = DownscalingUtils()
		self.write_anomalies = write_anomalies

	@staticmethod
	def _fn_month_grouper( fn, *args, **kwargs ):
		'''
		take a filename and return the month element of the naming convention
		'''
		return os.path.splitext( os.path.basename( fn ) )[0].split( '_' )[-2]
	def _get_varname_cru( self, *args, **kwargs ):
		'''
		take as input the cru ts3* netcdf filename and return (if possible)
		the name of the variable we want to work on from that netcdf.

		Arguments:
			nc_fn = [str] filepath to the cru ts* netcdf file used in downscaling

		Returns:
			the variable name as a string if it can be deduced, and errors if
			the variable name cannot be deduced.

		'''
		ds = xray.open_dataset( self.cru_ts )
		variables = ds.variables.keys()
		variable = [ variable for variable in variables \
						if variable not in [u'lon', u'lat', u'time'] ]
		if len( variable ) == 1:
			variable = variable[ 0 ]
		else:
			AttributeError( 'cannot deduce the variable from the file. supply nc_varname and re-run' )
		return variable
	def _get_version_cru( self, *args, **kwargs ):
		version = ''.join( os.path.basename( self.cru_ts ).split( '.' )[:2] )
		version = version.replace( 'ts', 'TS' ) # to follow convention
		return version
	def _calc_anomalies( self, variable, *args, **kwargs ):
		'''
		calculate absolute or relative anomalies given a NetCDF file
		of the Climatic Research Unit (CRU) Historical Time Series.
		'''
		import xray
		ds = xray.open_dataset( self.cru_ts )
		try:
			clim_ds = ds.loc[ {'time':slice(self.climatology_begin, self.climatology_end)} ]
			climatology = clim_ds[ variable ].groupby( 'time.month' ).mean( 'time' )
		except:
			AttributeError( 'cannot slice netcdf based on climatology years given. they must overlap.' )
		# calculate anomalies
		if self.absolute == True:
			anomalies = ds[ variable ].groupby( 'time.month' ) - climatology
		elif self.absolute == False:
			anomalies = ds[ variable ].groupby( 'time.month' ) / climatology
		else:
			AttributeError( '_calc_anomalies (cru): absolute can only be True or False' )
		return anomalies
	def _interp_downscale_wrapper( self, args_dict, *args, **kwargs  ):
		'''
		interpolate anomalies and downscale to the baseline arr
		'''
		output_filename = args_dict[ 'output_filename' ]
		args_dict.update( output_filename=output_filename.replace( 'downscaled', 'anom' ) )

		anom = self.utils.interpolate_anomalies( **args_dict )

		if isinstance( anom, basestring ):
			rst = rasterio.open( anom )
			meta = rst.meta
			meta.update( compress='lzw' )
			anom_arr = rst.read( 1 )
		elif isinstance( anom, tuple ):
			anom_arr, meta = anom
		else:
			AttributeError( '_interp_downscale_wrapper: passed wrong instance type' )

		args_dict.update( output_filename=output_filename, anom_arr=anom_arr, meta=meta )
		return self.utils.downscale( **args_dict )
	def downscale_cru_ts( self, *args, **kwargs ):
		'''
		run the CRU downscaling using the monthly climatology files given
		'''
		from pathos.mp_map import mp_map
		import glob, affine, rasterio

		nc_varname = self._get_varname_cru( )
		# handle cases where the desired varname != one parsed from file.
		if self.variable == None:
			variable = nc_varname
		else:
			variable = self.variable
		
		# build output dirs
		anomalies_path = os.path.join( base_path, variable, 'anom' )
		if not os.path.exists( anomalies_path ):
			os.makedirs( anomalies_path )

		downscaled_path = os.path.join( base_path, variable, 'downscaled' )
		if not os.path.exists( downscaled_path ):
			os.makedirs( downscaled_path )

		# template setup 
		template_raster = rasterio.open( self.template_raster_fn )
		template_meta = template_raster.meta
		template_meta.update( crs={'init':'epsg:3338'} )

		# make a mask with values of 0=nodata and 1=data
		template_raster_mask = template_raster.read_masks( 1 ) # mask of band 1 is all we need
		template_raster_mask[ template_raster_mask == 255 ] = 1

		anomalies = self._calc_anomalies( self.cru_ts, variable, absolute=self.absolute )
		anomalies_pcll, lons_pcll = self.utils.shiftgrid( 0., anomalies, anomalies.lon.data ) # grabs lons from the xray ds

		# mesh the lons and lats and unravel them to 1-D
		lo, la = [ i.ravel() for i in np.meshgrid( lons_pcll, anomalies.lat ) ]
		
		# convert into pandas.DataFrame and drop all the NaNs -- land-only dataset
		anom_df_list = [ pd.DataFrame({ 'anom':i.ravel(), 'lat':la, 'lon':lo }).dropna( axis=0, how='any' ) for i in anomalies_pcll ]
		xi, yi = np.meshgrid( lons_pcll, anomalies.lat.data )

		# argument setup -- HARDWIRED
		src_transform = affine.Affine( 0.5, 0.0, -180.0, 0.0, -0.5, 90.0 )
		src_nodata = -9999.0
			
		# output_filenames setup
		dates = ds.time.to_pandas()
		years = np.unique( dates.apply( lambda x: x.year ) ).tolist()
		# years = np.unique( self._get_years_cru( self.cru_ts ) ) # CHANGED!
		cru_ts_version = self._get_version_cru( self.cru_ts ) # works if naming convention stays same
		months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
		month_year = [ (month, year) for year in years for month in months ]

		# read in the pre-processed 12-month climatology
		clim_list = sorted( glob.glob( os.path.join( self.clim_path, '*.tif' ) ) ) # this could catch you.
		clim_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, clim_list ) }
		output_filenames = [ os.path.join( downscaled_path, '_'.join([ variable, self.metric, cru_ts_version, 'downscaled', month, str(year) ])+'.tif' )
								for month, year in month_year ]

		# set downscaling_operation based on self.absolute boolean
		if self.absolute == True:
			downscaling_operation = 'add'
		elif self.absolute == False:
			downscaling_operation = 'mult'
		else:
			AttributeError( 'downscaling operation: self.absolute must be boolean' )

		args_list = [ { 'anom_df':anom_df, 
						'meshgrid_tuple':(xi, yi), 
						'template_raster_fn':template_raster_fn, 
						'lons_pcll':lons_pcll, 
						'src_transform':src_transform, 
						'src_crs':self.src_crs, \
						'src_nodata':src_nodata, 
						'output_filename':out_fn,
						'baseline_arr':clim_dict[ self._fn_month_grouper( out_fn ) ],
						'downscaling_operation':downscaling_operation, 
						'post_downscale_function':self.post_downscale_function,
						'write_anomalies':self.write_anomalies }
							for anom_df, out_fn in zip( anom_df_list, output_filenames ) ]

		# run anomalies interpolation and downscaling in a single go.
		out = mp_map( lambda args: self._interp_downscale_wrapper( args_dict=args ), args_list, nproc=self.ncores )
		return 'downscaling complete. files output at: %s' % base_path

if __name__ == '__main__':
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	# example of use of the new DownscaleCRU / DownscalingUtils classes
	# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	# example of post_downscale_function - pass in at DownscaleCRU()
	def clamp_vals( x ):
		''' clamp the values following the relative humidity downscaling '''
		x[ (x > 100) & (x < 500) ] = 95
		return x

	# minimum required arguments
	cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
	clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
	template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
	base_path = '/atlas_scratch/malindgren/CMIP5'
	
	# run example
	down = DownscaleCRU( cru_ts, clim_path, template_raster_fn, base_path, absolute=False, ncores=32 )
	output = down.downscale_cru_ts()

