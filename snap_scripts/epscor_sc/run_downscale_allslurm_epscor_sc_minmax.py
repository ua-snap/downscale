if __name__ == '__main__':
	import subprocess, os, shutil
	from pathos.mp_map import mp_map

	# [ 2 ]TASMIN TASMAX DATA RUN
	# [ a ] RUN CRU first since it takes slightly longer with spatial interpolation
	done = subprocess.call([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/wrap_downscaler_cru_slurm_epscor_sc_minmax.py' ])

	# [ hack ] to work with the CCSM4 misnaming
	path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/CCSM4'
	if os.path.exists( path ):
		_ = subprocess.call([ 'mv', path, path.replace( 'CCSM4', 'NCAR-CCSM4' ) ])
		print( 'changed CCSM4 error to NCAR-CCSM4 for proper min/max handling' )
	else:
		BaseException( '{} is not present'.format( path ) )
	# end hack

	# [ b ] RUN CMIP5
	done = subprocess.call([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/wrap_downscaler_cmip5_slurm_epscor_sc_minmax.py' ])


