# zip all Monthly files in a single zip -- AKCAN 10min data to CKAN

def list_files( base_path, ext='.nc' ):
	return [ os.path.join(r,fn) for r,s,files in os.walk( base_path ) 
				for fn in files if fn.endswith(ext) ]

if __name__ == '__main__':
	import os, zipfile, zlib, itertools
	import pandas as pd
	import numpy as np

	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_10min_nwt'
	output_base_path = '/CKAN_Data/Base/AK_CAN_10min/projected/AR5_CMIP5_models'

	# hardwired variables to easily walk dirs... since we already know 'em.
	models = [ '5ModelAvg','GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4' ]
	scenarios = [ 'historical','rcp45','rcp60','rcp85' ]
	variables = [ 'tas','pr' ]

	for model, scenario, variable, in itertools.product(models, scenarios, variables):
		cur_path = os.path.join( base_path, model, scenario, variable )

		files = list_files( cur_path, ext='.tif' )
		if variable == 'tas':
			units = 'mean_C'
			varname = 'Temperature'
		elif variable == 'pr':
			units = 'total_mm'
			varname = 'Precipitation'

		# assumes full years... to cut out some annoying string-y code
		years = [ int(os.path.basename(fn).split('.')[0].split('_')[-1]) for fn in files ]
		begin = str(min(years))
		end = str(max(years))

		# [ hardwired ] update output_path
		output_path = os.path.join( output_base_path, f'Projected_Monthy_and_Derived_{varname}_Products_10min_CMIP5_AR5', 'monthly' )

		if not os.path.exists(output_path):
			os.makedirs(output_path)

		out_fn = os.path.join( output_path, f'{variable}_{units}_AK_CAN_AR5_{model}_{scenario}_01_{begin}-12_{end}.zip' )
		print(out_fn)

		# open a new zipfile
		zf = zipfile.ZipFile( out_fn, mode='w' )
			
		for fn in files:
			zf.write( fn, arcname=os.path.basename(fn), compress_type=zlib.DEFLATED )
		
		# cleanup
		zf.close()
		zf = None
		del zf


# /CKAN_Data/Base/AK_CAN_10min/projected/AR5_CMIP5_models/Projected_Monthy_and_Derived_Temperature_Products_10min_CMIP5_AR5/monthly
# /CKAN_Data/Base/AK_CAN_2km/projected/AR5_CMIP5_models/Projected_Monthy_and_Derived_Temperature_Products_2km_CMIP5_AR5/monthly

