# wrap downscaler for running on slurm
def run_model( fn, variable, model, scenario ):
	import os, subprocess
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscale_epscor_se_distributedrun.py', '--', '-m', model, '-v', variable, '-s', scenario  ])
		f.writelines( "#!/bin/sh\n#SBATCH --ntasks=32\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=32\n#SBATCH --account=snap\n#SBATCH --mail-type=all\n#SBATCH --mail-user=malindgren@alaska.edu\n#SBATCH -p main\n\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, glob, itertools, subprocess

	models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]
	variables = [ 'tasmin', 'tasmax' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]

	path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_v2/slurm_log'
	if not os.path.exists( path ):
		os.makedirs( path )
	
	os.chdir( path )
	for variable, model, scenario in itertools.product( variables, models, scenarios ):
		fn = os.path.join( path, 'slurm_run_downscaler_'+'_'.join([variable, model, scenario])+'.slurm' )
		_ = run_model( fn, variable, model, scenario )