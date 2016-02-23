#!/usr/bin/python2

import os, glob

os.chdir( '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/CODE/tem_ar5_inputs' )
base_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/prepped'
output_base_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/downscaled'
cru_base_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_ts20/akcan'

for root, dirs, files in os.walk( base_dir ):
	
	if files:
		path, variable = os.path.split( root )
		path, model = os.path.split( path )
	
		# this gets rid of any .xml or .txt files that may be living amongst the NetCDF's
		files = [ fn for fn in files if fn.endswith( '.nc' ) ]

		for fn in files:
			print 'running %s' % fn
			# split out the sub_dirs to have both model_name and variable folder hierarchy
			# from the prepped folder directory
			output_dir = os.path.join( output_base_dir, model, variable )

			if not os.path.exists( output_dir ):
				os.makedirs( output_dir )

			# anomalies calculation type and cru input path condition
			if 'tas_' in os.path.basename( fn ):
				anomalies_calc_type = 'absolute'
				downscale_operation = 'add'
				cru_path = os.path.join( cru_base_dir, 'tas' )
			elif 'hur_' in os.path.basename( fn ):
				anomalies_calc_type = 'proportional'
				downscale_operation = 'mult'
				cru_path = os.path.join( cru_base_dir, 'hur' )
				plev = 1000
			else:
				NotImplementedError( "only 'hur' & 'tas' have been implemented here" )

			# condition to determine if we need to read in the historical dataset with the modeled for 
			#  anomalies calculation
			if 'historical' in fn:
				# run with only the historical file
				dates = os.path.basename( fn ).strip( '.nc' ).split( '_' )
				dates = dates[ len( dates )-2 ], dates[ len( dates )-1 ]
				begin_time, end_time = [ '-'.join([ i[:4], i[4:] ]) for i in dates ]

				if 'tas_' in fn:
					os.system( 'python hur_ar5_model_data_downscaling.py' + ' -hi ' + os.path.join( root, fn ) + ' -o ' + output_dir + ' -bt ' + begin_time + \
									' -et ' + end_time + ' -cbt ' + '1961-01' + ' -cet ' + '1990-12' + \
									' -cru ' + cru_path + ' -at ' + anomalies_calc_type + ' -m ' + 'mean' + ' -dso ' + downscale_operation )
				elif 'hur_' in fn:
					# run with only the historical file
					os.system( 'python hur_ar5_model_data_downscaling.py' + ' -hi ' + os.path.join( root, fn ) + ' -o ' + output_dir + ' -bt ' + begin_time + \
									' -et ' + end_time + ' -cbt ' + '1961-01' + ' -cet ' + '1990-12' + \
									' -plev ' + str(plev) + ' -cru ' + cru_path + ' -at ' + anomalies_calc_type + ' -m ' + 'mean' + ' -dso ' + downscale_operation )
				else:
					NotImplementedError( "only 'hur' & 'tas' have been implemented here" )
			else:
				# grab the historical file from that particular folder
				historical_fn = glob.glob( os.path.join( root, '*historical*.nc' ) )[0]
				# run with both historical and modeled files for anomalies calc.
				if 'tas_' in fn:
					os.system( 'python hur_ar5_model_data_downscaling.py' + ' -mi ' + os.path.join( root, fn ) + ' -hi ' + historical_fn + ' -o ' + output_dir + \
									' -bt ' + '2006-01' + ' -et ' + '2100-12' + ' -cbt ' + '1961-01' + ' -cet ' + '1990-12' + \
									' -cru ' + cru_path + ' -at ' + anomalies_calc_type + ' -m ' + 'mean' + ' -dso ' + downscale_operation )
				elif 'hur_' in fn:
					os.system( 'python hur_ar5_model_data_downscaling.py' + ' -mi ' + os.path.join( root, fn ) + ' -hi ' + historical_fn + ' -o ' + output_dir + \
									' -bt ' + '2006-01' + ' -et ' + '2100-12' + ' -cbt ' + '1961-01' + ' -cet ' + '1990-12' + ' -plev ' + str(plev) + \
									' -cru ' + cru_path + ' -at ' + anomalies_calc_type + ' -m ' + 'mean' + ' -dso ' + downscale_operation )
				else:
					NotImplementedError( "only 'hur' & 'tas' have been implemented here" )

