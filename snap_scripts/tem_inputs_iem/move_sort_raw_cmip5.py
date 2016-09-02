def name_splitter( fn ):
	''' split filename into its parts and return in a dict (including full filename) '''
	import os
	dirname, basename = os.path.split( fn )
	basename, ext = os.path.splitext( basename )
	elems = ['variable', 'cmor_table', 'model', 'scenario', 'experiment', 'time']
	d = dict( zip( elems, basename.split( '_' ) ) )
	d[ 'fn' ] = fn
	return d

def copy_fn( fn, output_dir ):
	''' copy fn from current location to output_dir/model/scenario/variable/fn'''
	import shutil
	d = name_splitter( fn )
	dirname, basename = os.path.split( d['fn'] )
	out_fn = os.path.join( output_dir, d['model'], d['scenario'], d['variable'], basename )
	try: # mostly for parallel collisions
		if not os.path.exists( os.path.dirname( out_fn ) ):
			os.makedirs( os.path.dirname( out_fn ) )
	except:
		pass
	return shutil.copy( fn, out_fn )

if __name__ == '__main__':
	import os, glob
	from pathos.mp_map import mp_map
	from functools import partial

	base_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/raw/cmip5' # /output1/NASA-GISS/GISS-E2-R/historical/mon/atmos/Amon/r1i1p1/v20121015/hur
	output_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/raw_clean'

	filelist = []
	for root, subs, files in os.walk( base_dir ):
		if len( files ) > 0:
			filelist = filelist + [ os.path.join( root, fn ) for fn in files if fn.endswith( '.nc' ) ]
	f = partial( copy_fn, output_dir=output_dir )
	done = mp_map( f, filelist, nproc=32 )

# REACCESS WITH SYNDA ALL THE FILES WE NEED, THIS WAY WE CAN AUTOPURGE THE OLD VERSIONS WE HAVE HERE.
# project=CMIP5
# model=MRI-CGCM3 GISS-E2-R GFDL-CM3 IPSL-CM5A-LR CCSM4
# experiment=rcp26
# ensemble=r1i1p1
# variable[atmos][mon]=clt
# timeslice=1800-2150
