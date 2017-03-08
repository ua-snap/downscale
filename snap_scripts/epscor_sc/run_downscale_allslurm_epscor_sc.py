# # # run the 2 wrapper scripts that fire off jobs on ATLAS -- run this script from atlas head node.
def move( in_fn, out_fn ):
	import os, shutil
	dirname, basename = os.path.split( out_fn )
	try:
		if not os.path.exists( dirname ):
			os.makedirs( dirname )
	except:
		pass
	return shutil.move( in_fn, out_fn )

if __name__ == '__main__':
	import subprocess, os, shutil
	from pathos.mp_map import mp_map

	# [ 1 ] TAS / PR
	# [ a ]RUN CRU first since it takes slightly longer with spatial interpolation
	done = subprocess.call([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/wrap_downscaler_cru_slurm_epscor_sc.py' ])

	# [ b ]RUN CMIP5
	done = subprocess.call([ 'ipython', '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc/wrap_downscaler_cmip5_slurm_epscor_sc.py' ])

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

	# [ hack ] to move the misnamed folders files to the proper folder...
	dirname = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled'
	old = os.path.join( dirname, 'CCSM4' )
	new = os.path.join( dirname, 'NCAR-CCSM4' )
	filelist = [ (os.path.join( root, fn ), os.path.join( root.replace( old, new ) , fn ) ) for root, subs, files in os.walk( old ) for fn in files if fn.endswith( '.tif' ) ]

	out = mp_map( lambda x: move( *x ), filelist, nproc=32 )
	_ = subprocess.call([ 'rm', '-rf', old ])
	print( 'files moved' )
	# end hack