# # # EXAMPLE RUN SINGLE MODEL
def run_model( fn, command ):
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
		f.writelines( head + '\n' + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess

	# # args setup
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	ncores = '32'
	model = 'ts323'
	scenario = 'historical'
	variables = ['pre'] #['tmp', 'pre'] #[ 'tmx','tmn' ]
	out_varnames = ['pr'] #['tas', 'pr'] #[ 'tasmax','tasmin' ]
	
	slurm_path = os.path.join( base_dir, 'downscaled_FINAL_OCT','slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )

	os.chdir( slurm_path )

	for variable, out_varname in zip( variables, out_varnames ):
		if variable == 'pre':
			metric = 'total'
			units = 'mm'
		else:
			metric = 'mean'
			units = 'C'

		clim_path = os.path.join( base_dir, 'prism', out_varname )
		output_path = os.path.join( os.path.join( base_dir, 'downscaled_FINAL_OCT', model, scenario, out_varname ) )
		
		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.' + variable + '.dat.nc'
		
		# # make a command to pass to slurm
		script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/downscale_cru_epscor_sc.py'
		command = ' '.join([ 'ipython', script_path, '--',
							'-ts', cru_ts,
							'-cl', clim_path,
							'-o', output_path,
							'-m', model,
							'-v', variable,
							'-u', units,
							'-met', metric,
							'-nc', ncores,
							'-ov', out_varname ])
		
		fn = os.path.join( slurm_path, '_'.join(['downscale', model, variable]) + '.slurm' )
		_ = run_model( fn, command )
