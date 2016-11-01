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
	
	base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data'
	scripts_directory = '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc'
	ncpus = 32
	variables = [ 'pr', 'tasmin', 'tasmax', 'tas' ]
		
	# # cmip5
	project = 'ar5'
	models = [ 'GFDL-CM3','IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4', '5ModelAvg' ]
	scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	agg_metric = 'mean'

	# # cru
	# project = 'cru'
	# models = [ 'ts323' ]
	# scenarios = [ 'historical' ]

	for model in models:
		for scenario in scenarios:
			for variable in variables:
				# this deals with the 2 flavors of metrics for pr somewhat cleanly...
				# agg_metrics_lookup = { 'pr':['mean','total'], 'tas':['mean'], 'tasmin':['mean'], 'tasmax':['mean'] }
				# agg_metrics = agg_metrics_lookup[ variable ]
				# for agg_metric in agg_metrics: 
				# some stuff for output slurm log files
				path = os.path.join( base_dir, 'derived_grids', 'slurm_log' )
				
				if not os.path.exists( path ):
					os.makedirs( path )
				
				os.chdir( path )

				fn = os.path.join( path, 'slurm_run_post_downscaling_'+'_'.join([variable, model, scenario, agg_metric])+'.slurm' )

				# if pr -- we want annual/seasonal totals...
				if variable == 'pr':
					agg_metric = 'total'

				# annuals
				annuals = ' '.join([ 'python', os.path.join( scripts_directory, 'derived_grids_annuals.py' ),
										'-b', base_dir, '-m', model , '-s', scenario, '-p', project,
										'-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
				# annual seasonals
				annual_seasonals = ' '.join([ 'python', os.path.join( scripts_directory, 'derived_grids_annual_seasonals.py' ),
												'-b', base_dir, '-m', model , '-s', scenario, '-p', project,
												'-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])
				
				# turn back to 'mean' for remainder of series
				agg_metric = 'mean'

				# decadal annuals
				decadal_annuals = ' '.join([ 'python', os.path.join( scripts_directory, 'derived_grids_decadal_annuals.py' ),
										'-b', base_dir, '-m', model , '-s', scenario, '-p', project,
										'-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])

				# decadal monthlies
				decadal_monthlies = ' '.join([ 'python', os.path.join( scripts_directory, 'derived_grids_decadal_monthlies.py' ),
										'-b', base_dir, '-m', model , '-s', scenario, '-p', project,
										'-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])

				# decadal seasonals
				decadal_seasonals = ' '.join([ 'python', os.path.join( scripts_directory, 'derived_grids_decadal_annual_seasonals.py' ),
										'-b', base_dir, '-m', model , '-s', scenario, '-p', project,
										'-v', variable ,'-am', agg_metric ,'-nc', str(ncpus) ])

				# combine command -- for if variable != tas
				command = '\n\n'.join([ annuals, annual_seasonals, decadal_monthlies, decadal_seasonals, decadal_annuals ]) + '\n'

				# swi - summer warmth index
				if variable == 'tas':
					# use the monthly decadals to compute swi decadal
					monthly_decadals_path = os.path.join( base_dir, 'derived_grids', 'decadal_monthlies' )
					output_path = os.path.join( base_dir, 'derived_grids', 'decadal_swi' )
					decadal_swi = ' '.join([ 'python', os.path.join( scripts_directory, 'swi_calc_decadal_epscor_sc.py' ), \
												'-b', monthly_decadals_path, '-o', output_path, '-m', model , '-s', scenario, '-p', project, '-v', variable ])
					# add swi to command
					command = '\n\n'.join([ decadals, annual_seasonals, decadal_seasonals, decadal_swi ]) + '\n'

				# run the command using slurm on ATLAS
				run_model( fn, command )