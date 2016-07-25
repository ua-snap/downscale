# wrap downscaler for running on slurm
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
	import subprocess, os
	
	base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_v2/CRU_TS323'
	# base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_v2'
	output_base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_v2'
	ncpus = 32
	project = 'cru'
	# project = 'cmip5' # 'cru'
	variables = [ 'tasmin', 'tasmax', 'tas', 'pr' ]
	models = [ 'ts323' ]
	# models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg' ]
	scenarios = [ 'historical' ]
	# scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	for model in models:
		for scenario in scenarios:
			for variable in variables:
				# this deals with the 2 flavors of metrics for pr somewhat cleanly...
				agg_metrics_lookup = { 'pr':['mean','total'], 'tas':['mean'], 'tasmin':['mean'], 'tasmax':['mean'] }
				agg_metrics = agg_metrics_lookup[ variable ]
				for agg_metric in agg_metrics: 
					# some stuff for output slurm log files
					path = os.path.join( output_base_path, 'slurm_log' )
					
					if not os.path.exists( path ):
						os.makedirs( path )
					
					os.chdir( path )

					fn = os.path.join( path, 'slurm_run_post_downscaling_'+'_'.join([variable, model, scenario, agg_metric])+'.slurm' )
					
					# decadals
					output_path = os.path.join( output_base_path, 'monthly_decadals' )
					# decadals = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
					decadals = ' '.join([ 'python', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_decadal_grids_epscor_se.py', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
					
					# annual seasonals
					output_path = os.path.join( output_base_path, 'annual_seasonals' )
					# annual_seasonals = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_seasonal_annual_grids_epscor_se.py', '--', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
					annual_seasonals = ' '.join([ 'python', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_seasonal_annual_grids_epscor_se.py', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
					
					# decadal seasonals
					output_path = os.path.join( output_base_path, 'decadal_seasonals' )
					# decadal_seasonals = ' '.join([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_seasonal_decadal_grids_epscor_se.py', '--', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
					decadal_seasonals = ' '.join([ 'python', '/workspace/UA/malindgren/repos/downscale/snap_scripts/compute_seasonal_decadal_grids_epscor_se.py', '-b', base_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])

					command = '\n\n'.join([ decadals, annual_seasonals, decadal_seasonals ]) + '\n'
					run_model( fn, command )
