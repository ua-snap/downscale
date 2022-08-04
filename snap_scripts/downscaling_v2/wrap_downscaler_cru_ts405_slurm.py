# # # EXAMPLE RUN SINGLE MODEL
def run_model( fn, command ):
	import os, subprocess
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=kmredilla@alaska.edu\n' + \
			'#SBATCH -p viz\n\n' + \
			'eval "$(conda shell.bash hook)"\nconda activate\nexport OMP_NUM_THREADS=1'
	
	with open( fn, 'w' ) as f:
		f.writelines( head + '\n' + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess

	# # args setup
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	ncores = '32'
	model = 'ts405'
	scenario = 'historical'
	variables = ['tmp', 'pre']
	out_varnames = ['tas', 'pr']
	
	slurm_path = os.path.join( base_dir, 'downscaled','slurm_log' )
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

		clim_path = os.path.join( base_dir, 'climatologies', 'prism', out_varname )
		output_path = os.path.join( os.path.join( base_dir, 'downscaled', model, scenario, out_varname ) )
		
		if not os.path.exists( output_path ):
			os.makedirs( output_path )

		cru_ts = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS405/cru_ts4.05.1901.2020.' + variable + '.dat.nc'
		
		# # make a command to pass to slurm
		script_path = '/workspace/UA/kmredilla/downscale/snap_scripts/downscaling_v2/downscale_cru.py'
		command = ' '.join([ 'python', script_path,
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
        
# these commands run after this script was executed, to get these data on CKAN:
# delete this when this process is implemented properly. 
# mv /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/pr/anom/* /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/anomalies/CRU-TS405/historical/pr
# rm -r /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/pr/anom/

# mv /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/tas/anom/* /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/anomalies/CRU-TS405/historical/tas
# rm -r /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/tas/anom/

# zip -r /workspace/CKAN/CKAN_Data/Base/AK_CAN_2km/historical/CRU_TS/Historical_Monthly_and_Derived_Precipitation_Products_2km_CRU_TS/pr_AK_CAN_2km_CRU_TS405_historical.zip /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/pr/

# zip -r /workspace/CKAN/CKAN_Data/Base/AK_CAN_2km/historical/CRU_TS/Historical_Monthly_and_Derived_Temperature_Products_2km_CRU_TS/tas_AK_CAN_2km_CRU_TS405_historical.zip /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/ts405/historical/tas/
