# # FROM STEPHANIE BACK IN THE DAY...
# 1. Calculate the saturated vapor pressure (ew) as  ew = 6.112* exp(17.62*T/(243.12+T)    with T in [°C] and ew in [hPa] from CRU temperature at 0.5 x 0.5 degree resolution*
# 2. Calculate %RH at 0.5 x 0.5 degree resolution as e/ew where e is vapor pressure in hPa from CRU
# 3. Replace %RH >100 with 95 (see my last email for why this happens)
# 4. Interpolate %RH to 1km
# 5. Calculate vapor pressure at 1km resolution from interpolated RH and 1km temperature


def convert_to_hur( tas_arr, vap_arr ):
	esa_arr = 6.112 * np.exp( 17.62 * tas_arr/ (243.12 + tas_arr) )
	# esa_arr = 6.112 * np.exp( 22.46 * tas_arr / (272.62 + tas_arr) )
	hur_arr = vap_arr/esa_arr * 100
	return hur_arr

if __name__ == '__main__':
	# convert the vap/tas to hur for the GD CRU data
	import xarray as xr
	import numpy as np

	# filenames
	vap_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.vap.dat.nc'
	tas_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.tmp.dat.nc'

	# open the cru data
	vap = xr.open_dataset( vap_fn )
	tas = xr.open_dataset( tas_fn )
	
	# slice them to the variable we want and return the 3D array	
	v = vap.vap.data
	t = tas.tmp.data

	# # FROM STEPHANIE BACK IN THE DAY...
	# 1. Calculate the saturated vapor pressure (ew) as  ew = 6.112* exp(17.62*T/(243.12+T)    with T in [°C] and ew in [hPa] from CRU temperature at 0.5 x 0.5 degree resolution*
	# 2. Calculate %RH at 0.5 x 0.5 degree resolution as e/ew where e is vapor pressure in hPa from CRU
	h = convert_to_hur( t, v )

	# 3. Replace %RH >100 with 95 (see my last email for why this happens)
	h[ (~np.isnan(h)) & (h < 0) ] = 0
	h[ (~np.isnan(h)) & (h > 100) ] = 95

	# write this to disk:
	hur = vap.vap.copy()
	# update the DataArray attributes since we updated the data to a new variable
	hur.attrs.update( long_name='relative humidity', 
					units='pct', derived_by='SNAP - 12/8/2016', 
					derived_author='Michael Lindgren (malindgren@alaska.edu)' )

	# convert to an xarray dataset
	hur = hur.to_dataset( name='hur' )
	
	# put the data from above into the object.
	hur[ 'hur' ] = (('time', 'lat', 'lon' ), h )
	
	# update the Dataset attributes
	hur.attrs.update( COMMENTS='Variable computed by Michael Lindgren at SNAP 12/8/2016 using cru tmp and cru vap', equation_step1='esa_arr = 6.112 * np.exp( 17.62 * tas_arr / (243.12 + tas_arr) )', equation_step2='hur_arr = vap_arr/esa_arr * 100', post_process_step='values < 0 were set to 0.  values < 100 were set to 95 [per Stephanie McAfee suggestion]' )
	
	# write to disk
	output_filename = vap_fn.replace( 'vap', 'hur' ).replace( '.nc', '_snap_conversion.nc' )
	hur.to_netcdf( path=output_filename, mode='w' )

	# 4. Interpolate %RH to 1km
	# ---> RUN DOWNSCALE

	# 5. Calculate vapor pressure at 1km resolution from interpolated RH and 1km temperature

