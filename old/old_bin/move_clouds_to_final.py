# # # #
# a tool to move the files to the needed locations in the final directory
# # # #

import os, glob, itertools, shutil
from pathos.mp_map import mp_map
from functools import partial

models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
variables = [ 'rsds' ]

# some pathing
input_path = '/Data/malindgren/cru_november_final/IEM/clouds/ar5'
output_path = '/Data/malindgren/cru_november_final/final/IEM/ar5'

# make all combinations of the output variables
combinations = itertools.product( models, scenarios, variables )

for model, scenario, variable in combinations:
	print( ' '.join([ 'runnning: ', model, scenario]) )
	l = glob.glob( os.path.join( input_path, model, variable, 'downscaled', scenario, '*.tif' ) )
	out = os.path.join( output_path, model, variable, scenario )

	# remove them if they exist with the pct in the name
	ol = glob.glob( os.path.join( out, '*rsds_*_pct_*.tif' ) )
	if len( ol ) > 0:
		_ = mp_map( lambda x: os.unlink( x ), ol, nproc=32 )

	f = partial( shutil.move, dst=out )

	_ = mp_map( f, l, nproc=32 )



# # # # CHANGE METRIC NAME # # # # # # 
# # # #
# a tool to move the files to the needed locations in the final directory
# # # #

import os, glob, itertools, shutil
from pathos.mp_map import mp_map
from functools import partial

models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
variables = [ 'rsds' ]

# some pathing
input_path = '/Data/malindgren/cru_november_final/final/IEM/ar5'

# make all combinations of the output variables
combinations = itertools.product( models, scenarios, variables )

for model, scenario, variable in combinations:
	print( ' '.join([ 'runnning: ', model, scenario]) )
	l = glob.glob( os.path.join( input_path, model, variable, scenario, '*.tif' ) )
	ol = [ i.replace( 'Wm2', 'MJ-m2-d1'	) for i in l ]
	_ = mp_map( lambda x: os.rename( x[0], x[1] ), zip( l, ol ), nproc=32 )


	
# # # # # # MOVE TO A NEW DIR STRUCTURE TO ZIP # # # # 


import os, glob, itertools, shutil
from pathos.mp_map import mp_map
from functools import partial
import zipfile

models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
variables = [ 'rsds' ]

# some pathing
input_path = '/Data/malindgren/cru_november_final/final/IEM/ar5'
output_path = '/Data/malindgren/cru_november_final/final/for_zip'

# make all combinations of the output variables
combinations = itertools.product( variables, models, scenarios )

def run( x ):
	variable, model, scenario = x
	
	if variable == 'rsds':
		metric = 'MJ-m2-d1'
	elif variable == 'vap':
		metric = 'hPa'
	else:
		AttributeError( 'variable not supported' )

	break 

	l = glob.glob( os.path.join( input_path, model, variable, scenario, '*.tif' ) )
	
	years = [ int(os.path.basename(i).split('_')[-1].split( '.' )[0]) for i in l ]
	year_begin = min( years )
	year_end = max( years )

	out_zipfile_fn = os.path.join( output_path, model, variable, '_'.join([variable, 'mean', metric, 'iem', model, scenario, str(year_begin), str(year_end) ])+'.zip' )

	try:
		if not os.path.exists( os.path.dirname( out_zipfile_fn ) ):
			os.makedirs( os.path.dirname( out_zipfile_fn ) )
	except:
		pass

	with zipfile.ZipFile( out_zipfile_fn, mode='w', allowZip64=True  ) as myzip:
		os.chdir( os.path.join( input_path, model, variable, scenario ) )
		for i in l:
			myzip.write( os.path.basename(i), compress_type=zipfile.ZIP_DEFLATED )

	return out_zipfile_fn

_ = mp_map( run, combinations, nproc=32 )


# for variable, model, scenario in combinations:
# 	print( ' '.join([ 'runnning: ', variable, model, scenario]) )

# 	if variable == 'rsds':
# 		metric = 'MJ-m2-d1'
# 	elif variable == 'vap':
# 		metric = 'hPa'
# 	else:
# 		AttributeError( 'variable not supported' )

# 	break 

# 	l = glob.glob( os.path.join( input_path, model, variable, scenario, '*.tif' ) )
	
# 	years = [ int(os.path.basename(i).split('_')[-1].split( '.' )[0]) for i in l ]
# 	year_begin = min( years )
# 	year_end = max( years )

# 	out_zipfile_fn = os.path.join( output_path, model, variable, '_'.join([variable, 'mean', metric, 'iem', model, scenario, str(year_begin), str(year_end) ])+'.zip' )

# 	try:
# 		if not os.path.exists( os.path.dirname( out_zipfile_fn ) ):
# 			os.makedirs( os.path.dirname( out_zipfile_fn ) )
# 	except:
# 		pass

# with zipfile.ZipFile( out_zipfile_fn, 'w' ) as myzip:
# 	os.chdir( os.path.join( input_path, model, variable, scenario ) )
# 	_ = [ myzip.write( os.path.basename(i), compress_type=zipfile.ZIP_DEFLATED ) for i in l ]




