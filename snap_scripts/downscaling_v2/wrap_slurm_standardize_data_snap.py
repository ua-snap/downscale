# # # # # 
# wrap standardize_outputs_oob_epsg_snap.py for running on slurm
# --> this is a bit of a problem due to the massive file i/o so 
# 	I ended up using a serial version named:
# 	`run_standardize_serialprocess.py` * will require some hand-holding...
# # # # # 

def run_model( fn, base_dir, output_dir ):
	head = '#!/bin/sh\n' + \
			'#SBATCH --ntasks=32\n' + \
			'#SBATCH --nodes=1\n' + \
			'#SBATCH --ntasks-per-node=32\n' + \
			'#SBATCH --account=snap\n' + \
			'#SBATCH --mail-type=FAIL\n' + \
			'#SBATCH --mail-user=malindgren@alaska.edu\n' + \
			'#SBATCH -p main\n'
	
	script_path = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2/standardize_outputs_oob_epsg_snap.py'
	with open( fn, 'w' ) as f:
		command = ' '.join([ 'ipython', script_path, '--', '-b', base_dir, '-o', output_dir ])
		f.writelines( head + "\n" + command + '\n' )
	subprocess.call([ 'sbatch', fn ])
	return 1

if __name__ == '__main__':
	import os, subprocess
	
	# setup vars
	base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data'
	variables = [ 'pr', 'tasmin', 'tasmax', 'tas', 'hur', 'vap', 'clt', 'rsds' ]
	# scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
	models = ['5ModelAvg','GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4','ts323','ts40']
	
	# get list of directories to standardize
	root_list = [ root for root, subs, files in os.walk( os.path.join( base_dir, 'downscaled' ) ) if os.path.split(root)[1] in models ]
	

	slurm_path = os.path.join( base_dir, 'downscaled', 'slurm_log' )
	if not os.path.exists( slurm_path ):
		os.makedirs( slurm_path )
	
	os.chdir( slurm_path )
	
	# launch
	for root in root_list:
		fn = os.path.join( slurm_path, 'slurm_run_standardize_data_snap_'+'_'.join( root.split( '/' )[-3:] )+'.slurm' )
		_ = run_model( fn, root, root )