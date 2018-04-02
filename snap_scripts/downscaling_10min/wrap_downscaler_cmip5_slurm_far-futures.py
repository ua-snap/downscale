# # # # # 
# wrap downscaler for running on slurm
# # # # # 
def run_model( fn, base_dir, variable, model, scenario, units, metric, begin, end ):
	import os, subprocess
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_10min/downscale_cmip5_far-future.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path,\
							 '--', '-b', base_dir, '-m', model, '-v', variable, '-s', scenario, '-u', units, '-met', metric, '-by', str(begin), '-ey', str(end)])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, glob, itertools, subprocess

	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	models = [ 'GFDL-CM3', 'IPSL-CM5A-LR', 'GISS-E2-R', 'NCAR-CCSM4', 'MRI-CGCM3' ]
	variables = [ 'tas','pr' ]
	scenarios = [ 'historical','rcp45','rcp60','rcp85' ]

	path = os.path.join( base_dir,'downscaled_10min','slurm_log' )
	if not os.path.exists( path ):
		os.makedirs( path )
	
	os.chdir( path )
	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		if variable == 'pr':
			units = 'mm'
			metric = 'total'
		else:
			units = 'C'
			metric = 'mean'

		if 'CCSM4' in model and variable == 'pr' and scenario != 'historical': # special case for the pr variable
			year_groups = [(2006, 2100),(2101,2300)]
			for year_group in year_groups:
				begin, end = year_group
				fn = os.path.join( path, 'slurm_run_downscaler_'+'_'.join([variable, model, scenario, '-'.join([str(i) for i in year_group])])+'.slurm' )
				_ = run_model( fn, base_dir, variable, model, scenario, units, metric, begin, end )
			
		else:
			begin, end = 0,0 # essentially None...  See handler in downscale script.
			fn = os.path.join( path, 'slurm_run_downscaler_'+'_'.join([variable, model, scenario])+'.slurm' )
			_ = run_model( fn, base_dir, variable, model, scenario, units, metric, begin, end )