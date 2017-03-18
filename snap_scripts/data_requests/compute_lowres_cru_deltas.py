# here we will generate the anomalies -- deltas -- for the researcher that wants them at low-res
def output_deltas( fn, varname, clim_begin, clim_end, method, output_dir ):
	'''
	ARGUMENTS:
	----------
	varname = [str] abbreviated name used for variable of interest in NC file.
	clim_begin = [str] 4 digit of start year
	clim_end = [str] 4 digit of end year
	method = [str] one of 'absolute' or 'relative'
	output_dir = [str] path to put the new output NC files into.
	'''
	import xarray as xr
	import os

	ds = xr.open_dataset( fn )
	ds = ds[ varname ]
	climatology = ds.sel( time=slice( clim_begin, clim_end ) )
	climatology = climatology.groupby( 'time.month' ).mean( 'time' )

	if method == 'absolute':
		anomalies = ds.groupby( 'time.month' ) - climatology
	elif method == 'relative':
		anomalies = ds.groupby( 'time.month' ) / climatology
	else:
		NameError( 'must be one of "absolute", "relative"' )
	basename, ext = os.path.splitext( os.path.basename( fn ) )
	out_fn = os.path.join( output_dir, '{}_anomalies{}'.format(basename, ext) )
	anom_ds = anomalies.to_dataset( name='{}_anomalies'.format( varname ) )
	anom_ds.to_netcdf( out_fn )
	return out_fn

if __name__ == '__main__':
	import os, glob
	import xarray as xr

	base_dir = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323'
	output_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cru_lowres_deltas'
	clim_begin = '1961'
	clim_end = '1990'

	# list the files
	l = glob.glob( os.path.join( base_dir, '*.nc' ) )
	varnames = [ i.split('.')[-3] for i in l ]

	out = []
	for fn, varname in zip( l, varnames ):
		if varname != 'wet':
			print( 'working on: {}'.format( varname ) )
			if varname in [ 'clt', 'hur', 'pre' ]:
				method = 'relative'
			else:
				method = 'absolute'
			
			out + [ output_deltas( fn, varname, clim_begin, clim_end, method, output_dir ) ]





