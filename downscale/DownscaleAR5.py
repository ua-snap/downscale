# -*- coding: utf8 -*-
# # #
# Downscale PCMDI AR5 data to a pre-processed climatology
#  extent, resolution, reference system
#
# Author: Michael Lindgren (malindgren@alaska.edu)
# # #
import rasterio, os
import numpy as np

class DownscaleAR5( object ):
	'''
	class used in downscaling the PCMDI CMIP3/5 modeled future or historical 
	screnarios data to a new extent resolution and AOI.
	'''

	def __init__( self, ar5_modeled=None, ar5_historical=None, base_path=None, clim_path=None, climatology_begin='1961', climatology_end='1990', \
		plev=None, absolute=True, metric='metric', variable=None, ncores=2, post_downscale_function=None, src_crs={'init':'epsg:4326'}, write_anomalies=True, \
		template_raster_fn=None, *args, **kwargs ):
		'''
		NEW METHODS FOR AR5 DOWNSCALING USING THE NEW
		API-ECOSYSTEM.
		'''
		from downscale import DownscalingUtils

		self.ar5_modeled = ar5_modeled
		self.ar5_historical = ar5_historical
		self.base_path = base_path
		self.clim_path = clim_path
		self.climatology_begin = climatology_begin
		self.climatology_end = climatology_end
		self.plev = plev
		self.absolute = absolute
		self.metric = metric
		self.variable = variable
		self.ncores = ncores
		self.utils = DownscalingUtils
		self.post_downscale_function = post_downscale_function
		self.src_crs = src_crs
		self.write_anomalies = write_anomalies
		self.template_raster_fn = template_raster_fn

	@staticmethod
	def standardized_fn_to_vars( fn ):
		''' take a filename string following the convention for this downscaling and break into parts and return a dict'''
		name_convention = [ 'variable', 'cmor_table', 'model', 'scenario', 'experiment', 'begin_time', 'end_time' ]
		fn = os.path.basename( fn )
		fn_list = fn.split( '.' )[0].split( '_' )
		return { i:j for i,j in zip( name_convention, fn_list ) }
	def _calc_anomalies( self, *args, **kwargs ):
		'''
		calculate absolute or relative anomalies given a NetCDF file
		of the Climatic Research Unit (CRU) Historical Time Series.
		'''
		import xray
		# handle modeled vs. historical
		if self.ar5_modeled != None and self.ar5_historical != None:
			# parse the input name for some file metadata HARDWIRED!
			output_naming_dict = self.standardized_fn_to_vars( self.ar5_modeled )
			variable = output_naming_dict[ 'variable' ]

			# read in both modeled and historical
			ds = xray.open_dataset( self.ar5_modeled )
			ds = ds[ variable ]
			clim_ds = xray.open_dataset( self.ar5_historical )
			# climatology
			clim_ds = clim_ds.loc[ {'time':slice(self.climatology_begin,self.climatology_end)} ]
			climatology = clim_ds[ variable ].groupby( 'time.month' ).mean( 'time' )
			del clim_ds

		elif self.ar5_historical is not None and self.ar5_modeled is None:
			output_naming_dict = standardized_fn_to_vars( self.ar5_historical )
			variable = output_naming_dict[ 'variable' ]

			# read in historical
			ds = xray.open_dataset( self.ar5_historical )
			# climatology
			climatology = ds.loc[ {'time':slice(self.climatology_begin,self.climatology_end)} ]
			climatology = climatology[ variable ].groupby( 'time.month' ).mean( 'time' )
		else:
			NameError( 'ERROR: must have both ar5_modeled and ar5_historical, or just ar5_historical' )

		if self.plev is not None:
			plevel, = np.where( ds.plev == self.plev )
			ds = ds[ :, plevel[0], ... ]
			climatology = climatology[ :, plevel[0], ... ]

		# anomalies
		if self.absolute == True:
			anomalies = ds.groupby( 'time.month' ) - climatology
		elif self.absolute == False:
			anomalies = ds.groupby( 'time.month' ) / climatology
		else:
			AttributeError( '_calc_anomalies (ar5): absolute can only be True or False' )
		return anomalies
	def _get_varname_ar5( self, *args, **kwargs ):
		'''
		take as input the AR5 netcdf filename and return (if possible)
		the name of the variable we want to work on from that netcdf.

		Arguments:
			nc_fn = [str] filepath to the AR5 ts* netcdf file used in downscaling

		Returns:
			the variable name as a string if it can be deduced, and errors if
			the variable name cannot be deduced.

		'''
		return os.path.basename( self.ar5_historical ).split( '_' )[0]
	def _calc_ar5_affine( self, *args, **kwargs ):
		'''
		this assumes 0-360 longitudes and
		WGS84 LatLong.
		'''
		import affine, xray
		ds = xray.open_dataset( self.ar5_modeled )
		lat_shape, lon_shape = ds.dims[ 'lat' ], ds.dims[ 'lon' ]
		lat_res = 180.0 / lat_shape
		lon_res = 360.0 / lon_shape
		return affine.Affine( lon_res, 0.0, 0.0, 0.0, -lat_res, 360.0 )
	@staticmethod
	def _fn_month_grouper( fn, *args, **kwargs ):
		'''
		take a filename and return the month element of the naming convention
		'''
		return os.path.splitext( os.path.basename( fn ) )[0].split( '_' )[-2]
	def _interp_downscale_wrapper( self, args_dict, *args, **kwargs  ):
		'''
		interpolate anomalies and downscale to the baseline arr
		'''
		import rasterio
		from rasterio.warp import RESAMPLING, reproject

		# unpack some of the args_dict
		output_filename = args_dict[ 'output_filename' ]
		anom_filename = output_filename.replace( 'downscaled', 'anom' )

		# UGLY!
		anom_arr = args_dict[ 'anom_arr' ]
		meshgrid_tuple = args_dict[ 'meshgrid_tuple' ]
		template_raster_fn = args_dict[ 'template_raster_fn' ]
		lons_pcll = args_dict[ 'lons_pcll' ]
		src_transform = args_dict[ 'src_transform' ]
		src_crs = args_dict[ 'src_crs' ]
		src_nodata = args_dict[ 'src_nodata' ]
		write_anomalies = args_dict[ 'write_anomalies' ]

		# read in the template raster
		template_raster = rasterio.open( template_raster_fn )
		template_meta = template_raster.meta

		if np.where( lons_pcll > 200.0 ).any() == True:
			# rotate globe back to -180.0 to 180.0 longitudes if needed
			dat, lons = self.utils.shiftgrid( 180., anom_arr, lons_pcll, start=False )
			output_arr = np.empty_like( template_raster.read( 1 ) )

		# reproject it
		reproject( dat, output_arr, src_transform=src_transform, src_crs=src_crs, src_nodata=src_nodata, \
			dst_transform=template_meta['affine'], dst_crs=template_meta['crs'],\
			dst_nodata=None, resampling=RESAMPLING.cubic_spline, SOURCE_EXTRA=1000 )

		# mask it with the internal mask in the template raster, where 0 is oob. DangerTown™
		mask = template_raster.read_masks( 1 ) == 0
		output_arr[ mask ] = template_meta[ 'nodata' ]

		# write or return anomalies...
		if write_anomalies == True:
			anom = self.utils.write_gtiff( output_arr, template_meta, anom_filename, compress=True )
		elif write_anomalies == False:
			anom = ( output_arr, template_meta )
		else:
			AttributeError( '_interp_downscale_wrapper: write_anomalies can be True or False only.' )

		# downscale -- handle output GTIff or not, this could be much cleaner
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
	def downscale_ar5_ts( self, *args, **kwargs ):
		#  * * * * * * * * * *
		# template setup
		from pathos.mp_map import mp_map
		import glob, affine, rasterio

		nc_varname = self._get_varname_ar5()
		
		# handle cases where the desired varname != one parsed from file.
		if self.variable == None:
			variable = nc_varname
		else:
			variable = self.variable
		
		print( 'downscaling: %s' % variable )

		# build output dirs
		anomalies_path = os.path.join( self.base_path, variable, 'anom' )
		if not os.path.exists( anomalies_path ):
			os.makedirs( anomalies_path )

		downscaled_path = os.path.join( self.base_path, variable, 'downscaled' )
		if not os.path.exists( downscaled_path ):
			os.makedirs( downscaled_path )

		#  * * * * * * * * * *
		# calc the anomalies
		anomalies = self._calc_anomalies()
	
		# mesh the lons and lats and unravel them to 1-D
		xi, yi = np.meshgrid( anomalies.lon.data, anomalies.lat.data )
		
		# some metadata
		src_transform = self._calc_ar5_affine()

		# flip it to the greenwich centering
		a,b,c,d,e,f,g,h,i = src_transform
		src_transform = affine.Affine( a, b, -180.0, d, e, 180.0 )

		# argument setup -- HARDWIRED
		src_nodata = None # DangerTown™

		# output_filenames setup
		dates = anomalies.time.to_pandas()
		years = np.unique( dates.apply( lambda x: x.year ) ).tolist()
		months = [ i if len(i)==2 else '0'+i for i in np.arange( 1, 12+1, 1 ).astype( str ).tolist() ]
		month_year = [ (month, year) for year in years for month in months ]

		# read in the pre-processed 12-month climatology
		clim_list = sorted( glob.glob( os.path.join( self.clim_path, '*.tif' ) ) ) # this could catch you.
		clim_dict = { month:rasterio.open( fn ).read( 1 ) for month, fn in zip( months, clim_list ) }
		
		# [!] THIS BELOW NEEDS RE-WORKING FOR THE AR5 DATA MODELED DATA # DangerTown™
		output_filenames = [ os.path.join( downscaled_path, '_'.join([ variable, self.metric, 'ar5', 'downscaled', month, str(year) ])+'.tif' )
								for month, year in month_year ]

		# set downscaling_operation based on self.absolute boolean
		if self.absolute == True:
			downscaling_operation = 'add'
		elif self.absolute == False:
			downscaling_operation = 'mult'
		else:
			AttributeError( 'downscaling operation: self.absolute must be boolean' )

		args_list = [ { 'anom_arr':anom_arr, 
						'meshgrid_tuple':(xi, yi), 
						'template_raster_fn':self.template_raster_fn, 
						'lons_pcll':anomalies.lon.data, 
						'src_transform':src_transform, 
						'src_crs':self.src_crs,
						'src_nodata':src_nodata,
						'output_filename':out_fn,
						'baseline_arr':clim_dict[ self._fn_month_grouper( out_fn ) ],
						'downscaling_operation':downscaling_operation, 
						'post_downscale_function':self.post_downscale_function,
						'write_anomalies':self.write_anomalies }
							for anom_arr, out_fn in zip( anomalies.data, output_filenames ) ]

		# run anomalies interpolation and downscaling in a single go.
		out = mp_map( lambda args: self._interp_downscale_wrapper( args_dict=args ), args_list, nproc=self.ncores )
		return 'downscaling complete. files output at: %s' % self.base_path