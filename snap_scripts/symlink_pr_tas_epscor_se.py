# symlink matts pr/tas data to the folder structure in the EPSCoR directory...
# cp -rs /Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models/rcp26/CCSM4/pr/ pr/

def run( x ):
	return os.system( x )

if __name__ == '__main__':
	import os, glob, subprocess, itertools
	from pathos.mp_map import mp_map

	output_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_v2'
	ncpus = 15
	project = 'cmip5' # 'cru'
	variables = ['tas', 'pr' ]
	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	
	commands = []
	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		if scenario == 'historical':
			base_path = '/Data/Base_Data/Climate/AK_CAN_2km/historical/AR5_CMIP5_models'
		else:
			base_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models'
		
		# print( '{} {} {}'.format( variable, model, scenario ) )
		base = os.path.join( base_path, scenario, model, variable )
		out = os.path.join( output_path, model, scenario )

		if not os.path.exists( out ):
			os.makedirs( out )

		# symlink them to their new directory structure
		commands = commands + [ ' '.join([ 'cp', '-rs', base, out ]) ]
		# _ = subprocess.call([ 'cp', '-rs', base, out ])
	
	final = mp_map( run, commands, nproc=ncpus )
