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
	ncores = '32'
	model = 'ts323'
	units = 'C'
	variables = ['tmx','tmn','tmp','pre']
	out_varnames = ['tasmax','tasmin', 'tas', 'pr']
	
	slurm_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2/slurm_files'

	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )

	os.chdir( slurm_path )

	for variable, out_varname in zip( variables, out_varnames ):
		if variable == 'pre':
			metric = 'total'
		else:
			metric = 'mean'

		if out_varname in ['tas','pr']:
			clim_path = os.path.join( '/Data/Base_Data/Climate/AK_CAN_2km/historical/singleBand/prism/AK_CAN_2km_PRISM/AK_CAN_geotiffs', out_varname, 'ak83albers' )
		elif out_varname in ['tasmin', 'tasmax']:
			if variable == 'tmx':
				v = 'tmax'
			if variable == 'tmn':
				v = 'tmin'
			clim_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/akcan', v )
		else:
			break
		
		output_path = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2/CRU_TS323', out_varname )
		cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.' + variable + '.dat.nc'
		
		# # make a command to pass to slurm
		command = ' '.join([ 'ipython','/workspace/UA/malindgren/repos/downscale/snap_scripts/run_cru_snap_epscor_se_CLI.py', '--',
														'-ts', cru_ts ,
														'-cl', clim_path ,
														'-o', output_path ,
														'-m', model ,
														'-v', variable ,
														'-u', units ,
														'-met', metric ,
														'-nc', ncores ,
														'-ov', out_varname ])
		
		
		fn = os.path.join( slurm_path, '_'.join(['downscale', model, variable]) + '.slurm' )
		_ = run_model( fn, command )
