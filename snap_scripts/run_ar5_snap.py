# AR5 Run -- January 2016

import glob, os, itertools
from downscale import DownscaleAR5

input_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data'
models = [ d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d)) ]
variables = ['hur', 'tas', 'clt']

combinations = itertools.product( [input_path], models, variables )

# static args setup
clim_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_cl20/cld/akcan'
template_raster_fn = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/templates/tas_mean_C_AR5_GFDL-CM3_historical_01_1860.tif'
base_path = '/Data/malindgren/downscale_outputs/AR5'
ncores = 16

for input_path, model, variable in combinations:
	nc_list = glob.glob( os.path.join( input_path, model, variable, '*.nc' ) )
	# unpack to historical and modeled
	historical, = [ nc for nc in nc_list if '_historical' in nc ]
	modeled = [ nc for nc in nc_list if '_historical' not in nc ]

	historical_modeled = zip( itertools.repeat( historical, len( modeled ) ), modeled )

	# run the historical first
	args = {}
	if variable == 'hur':

		def clamp_vals( x ):
			''' clamp the values following the relative humidity downscaling '''
			x[ (x > 100) & (x < 500) ] = 95
			return x

		args.update( ar5_historical=historical, base_path=base_path, clim_path=clim_path, \
					template_raster_fn=template_raster_fn, ncores=ncores, post_downscale_function=clamp_vals )
	else:
		args.update( ar5_historical=historical, base_path=base_path, clim_path=clim_path, \
					template_raster_fn=template_raster_fn, ncores=ncores )

	# run it
	print( historical )
	down = DownscaleAR5.DownscaleAR5( **args )
	output = down.downscale_ar5_ts()

	# now loop through the historical_modeled groupings
	for historical, modeled in historical_modeled:
		print( modeled )
		args = {}
		if variable == 'hur':

			def clamp_vals( x ):
				''' clamp the values following the relative humidity downscaling '''
				x[ (x > 100) & (x < 500) ] = 95
				return x

			args.update( ar5_modeled=ar5_modeled, ar5_historical=historical, base_path=base_path, clim_path=clim_path, template_raster_fn=template_raster_fn, ncores=ncores, post_downscale_function=clamp_vals )
		else:
			args.update( ar5_modeled=ar5_modeled, ar5_historical=historical, base_path=base_path, clim_path=clim_path, template_raster_fn=template_raster_fn, ncores=ncores )

		# run it
		down = DownscaleAR5.DownscaleAR5( **args )
		output = down.downscale_ar5_ts()






# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

# # move the prepped clouds data into the directory with the hur and tas
# import os, glob, shutil

# input_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped/clt_prepped'
# output_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped'

# # get the model names in the directory
# models = [ d for d in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, d)) ]

# variable = 'clt'

# for model in models:
# 	files = glob.glob( os.path.join( input_path, model, variable, '*.nc' ) )
# 	if not os.path.exists( os.path.join( output_path, model, variable ) ):
# 		os.makedirs( os.path.join( output_path, model, variable ) )
# 	_ = [ shutil.copy( fn, os.path.join( output_path, model, variable ) ) for fn in files ]













