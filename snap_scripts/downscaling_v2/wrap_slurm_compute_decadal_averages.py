# # # # # 
# wrap compute_decadal_averages.py for running on slurm
# # # # # 

def run( fn, path, output_path, ncpus ):
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main,viz\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2/compute_decadal_averages.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path, '--', '-p', path, '-o', output_path, '-n', ncpus ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess
	
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	models = []
	variables = [ 'pr', 'tasmin', 'tasmax', 'tas', 'hurs', 'vap', 'clt', 'rsds', ]
	ncpus = 32

	input_base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	output_base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_decadals'

	slurm_path = os.path.join( base_dir, 'downscaled_decadals', 'slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
	
	os.chdir( slurm_path )

	for model in models:
		scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85', ]
		if model == 'CRU-TS40':
			scenarios = [ 'historical', ]
		for variable in variables:
			for scenario in scenarios:
				fn = os.path.join( slurm_path, 'slurm_run_compute_decadals_'+variable+'_'+scenario+'.slurm' )
				path = os.path.join( input_base_path, model, scenario, variable )
				output_path = os.path.join( output_base_path, model, scenario, variable )
				_ = run( fn, path, output_path, ncpus )
