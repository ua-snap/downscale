def copy_file( in_fn, out_fn ):
	try:
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			_ = os.makedirs( dirname )
	except:
		pass

	return shutil.copy( in_fn, out_fn )

def run( x ):
	return copy_file( *x )

if __name__ == '__main__':
	import os, glob, shutil
	import multiprocessing as mp
	
	base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/iem_1km'
	out_path = '/atlas_scratch/ALFRESCO/ALFRESCO_Master_Dataset_v2_1/ALFRESCO_Model_Input_Datasets/IEM_for_TEM_inputs'

	# models
	models = [ 'MRI-CGCM3','NCAR-CCSM4' ]
	scenarios = [ 'rcp85' ] 

	# # cru
	# models = [ 'CRU-TS40' ]
	# scenarios = [ 'historical' ] 

	variables = [ 'hurs', 'rsds', 'vap' ]

	for variable in variables:
		for model in models:
			# for output naming to follow an odd convention.
			if model == 'CRU-TS40':
				group = 'iem'
			else:
				group = 'ar5'

			for scenario in scenarios:
				if model == 'CRU-TS40':
					years_suffix = '1901_2015'
				else:
					if scenario == 'historical':
						years_suffix = '1900_2005'
					else:
						years_suffix = '2006_2100'

				var_units_lookup = {'hurs':'pct', 'rsds':'MJ-m2-d1', 'vap':'hPa'}
				foldername = '{}_mean_{}_{}_{}_{}_{}_fix'.format(variable, var_units_lookup[variable], group, model, scenario, years_suffix)
				new_path = os.path.join( out_path, foldername, variable )
				
				if not os.path.exists( new_path ):
					_ = os.makedirs( new_path )

				old_path = os.path.join( base_path, model, scenario, variable )
				files = glob.glob( os.path.join( base_path, model, scenario, variable, '*.tif' ) )
				out_files = [ fn.replace(old_path,new_path) for fn in files ]

				args = list( zip( files, out_files ) )

				pool = mp.Pool( 64 )
				done = pool.map( run, args )
				pool.close()
				pool.join()



