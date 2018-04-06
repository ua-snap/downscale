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
	output_filename = vap_fn.replace( 'vap', 'hur' ).replace( '.nc', '_snap_conversion.nc' )
	hur.to_netcdf( path=output_filename, mode='w' )

	# 4. Interpolate %RH to 1km
	# ---> RUN DOWNSCALE

	# 5. Calculate vapor pressure at 1km resolution from interpolated RH and 1km temperature



T, °C	T, °F	P, kPa	P, torr	P, atm
0	32	0.6113	4.5851	0.0060
5	41	0.8726	6.5450	0.0086
10	50	1.2281	9.2115	0.0121
15	59	1.7056	12.7931	0.0168
20	68	2.3388	17.5424	0.0231
25	77	3.1690	23.7695	0.0313
30	86	4.2455	31.8439	0.0419
35	95	5.6267	42.2037	0.0555
40	104	7.3814	55.3651	0.0728
45	113	9.5898	71.9294	0.0946
50	122	12.3440	92.5876	0.1218
55	131	15.7520	118.1497	0.1555
60	140	19.9320	149.5023	0.1967
65	149	25.0220	187.6804	0.2469
70	158	31.1760	233.8392	0.3077
75	167	38.5630	289.2463	0.3806
80	176	47.3730	355.3267	0.4675
85	185	57.8150	433.6482	0.5706
90	194	70.1170	525.9208	0.6920
95	203	84.5290	634.0196	0.8342
100	212	101.3200	759.9625	1.0000
