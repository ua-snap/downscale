# # # # # 
# wrap downscaler for running on slurm
# # # # # 
def run_model( fn, base_dir, variable, model, scenario, units, metric ):
	import os, subprocess
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/tem_inputs_iem/downscale_cmip5_tem_iem.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path,\
							 '--', '-b', base_dir, '-m', model, '-v', variable, '-s', scenario, '-u', units, '-met', metric ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, glob, itertools, subprocess

	base_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016'
	models = [ 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'CCSM4' ]
	variables = [ 'tas','hur','pr' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

	path = os.path.join( base_dir,'downscaled','slurm_log' )
	if not os.path.exists( path ):
		os.makedirs( path )
	
	os.chdir( path )
	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		if variable == 'pr':
			units = 'mm'
			metric = 'total'
		elif variable == 'tas':
			units = 'C'
			metric = 'mean'
		elif variable == 'hur':
			units = 'pct'
			metric = 'mean'

		fn = os.path.join( path, 'slurm_run_downscaler_'+'_'.join([variable, model, scenario])+'.slurm' )
		_ = run_model( fn, base_dir, variable, model, scenario, units, metric )
