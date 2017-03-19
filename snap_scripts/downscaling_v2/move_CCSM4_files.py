# [ hack ] to move the misnamed folders files to the proper folder...
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
	import os, subprocess
	from pathos.mp_map import mp_map
	
	dirname = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	old = os.path.join( dirname, 'CCSM4' )
	new = os.path.join( dirname, 'NCAR-CCSM4' )
	filelist = [ (os.path.join( root, fn ), os.path.join( root.replace( old, new ) , fn ) ) for root, subs, files in os.walk( old ) for fn in files if fn.endswith( '.tif' ) ]

	out = mp_map( lambda x: move( *x ), filelist, nproc=32 )
	_ = subprocess.call([ 'rm', '-rf', old ])
	print( 'files moved' )
	
