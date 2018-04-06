from downscale import DeltaDownscale

class DeltaDownscaleMM( DeltaDownscale ):
	def _calc_anomalies( self ):
		print('calculating anomalies')	
	def downscale( self, *args, **kwargs ):
		print( 'downscaling...' )


# FOR RUN OF THE MIN / MAX TAS DATA:
# 1. COMPUTE DELTAS FIRST ANND WRITE TO NETCDF
# 2. USE `DeltaDownscaleMinMax` for the downscaling

def delta_mm( fn, mean_fn, variable, mean_variable, output_filename ):
	'''
	simple way to compute extreme - mean deltas as 
	native model resolution and write to NetCDF4 on disk
	'''
	import xarray as xr
	ds = xr.open_dataset( fn )[ variable ]
	ds_mean = xr.open_dataset( mean_fn )[ mean_variable ]
	delta = ds - ds_mean
	delta_ds = delta.to_dataset( name=variable )
	delta_ds.to_netcdf( output_filename )
	return output_filename

fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmx.dat.nc'
mean_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmp.dat.nc'
variable = 'tmx'
mean_variable = 'tmp'
output_filename = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/test/cru_ts3.23.1901.2014.tmx_delta_tmp.dat.nc'
delta_mm( fn, mean_fn, variable, mean_variable, output_filename )

# now use the new DeltaDownscaleMM class to do the work.



