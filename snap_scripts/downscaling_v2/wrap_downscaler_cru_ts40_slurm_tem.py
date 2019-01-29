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
			'#SBATCH -p viz\n'
	
	with open( fn, 'w' ) as f:
		f.writelines( head + '\n' + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess

	# # args setup
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	ncores = '32'
	model = 'ts40'
	scenario = 'historical'
	variables = ['hurs'] #['hur','cld'] #['tmp','hur','pre','cld']
	out_varnames = ['hurs'] #['hur','clt'] #['tas','hur','pr','clt']
	
	slurm_path = os.path.join( base_dir, 'downscaled','slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )

	os.chdir( slurm_path )

	for variable, out_varname in zip( variables, out_varnames ):
		if variable == 'pre':
			metric = 'total'
			units = 'mm'
		elif variable == 'tmp':
			metric = 'mean'
			units = 'C'
		elif variable in ['hurs','cld','clt']:
			metric = 'mean'
			units = 'pct'

		cru_cl20_varnames = {'hurs':'reh', 'clt':'clt', 'cld':'clt'} # we only support these variables for now...
		clim_path = os.path.join( base_dir, 'climatologies','cru_cl20', '2km', cru_cl20_varnames[out_varname] )
		output_path = os.path.join( os.path.join( base_dir, 'downscaled', model, scenario, out_varname ) )
		
		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS40/cru_ts4.00.1901.2015.' + variable + '.dat.nc.gz'
		if variable == 'hurs': # since we made this variable and it lives with the raw files with a slightly diff name
			cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS40/cru_ts4.00.1901.2015.' + variable + '.dat_snap_conversion.nc'
		
		# # make a command to pass to slurm
		script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2/downscale_cru_tem.py'
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
