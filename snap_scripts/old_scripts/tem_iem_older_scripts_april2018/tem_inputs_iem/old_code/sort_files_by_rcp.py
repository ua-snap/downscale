import os, glob, shutil
from pathos import multiprocessing as mp
import pandas as pd
import numpy as np

base_path = '/Data/malindgren/cru_november_final/IEM/ar5' 
output_base_path = '/Data/malindgren/cru_november_final/IEM/ar5'

models = [ 'IPSL-CM5A-LR', 'GISS-E2-R', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
# variables = ['rsds', 'vap' ]

for model in models:
	variables = os.listdir( os.path.join( base_path, model ) )
	_ = [ os.makedirs( os.path.join( base_path, model, variable ) ) for variable in variables if not os.path.exists( os.path.join( base_path, model, variable ) ) ]
	for variable in variables:
		print( ' '.join([model, variable]) )
		output_path = os.path.join( output_base_path, model, variable, 'downscaled' )
		cur_path = os.path.join( base_path, model, variable, 'downscaled'  )
		l = pd.Series( glob.glob( os.path.join( cur_path, '*.tif' ) ) )
				
		grouper = [ os.path.basename(i).split( '_' )[ 5 ] for i in l ]
		rcp_groups = l.groupby( grouper )
		name_group = [ group for group in rcp_groups ]

		names = [ i[0] for i in name_group ]

		_ = [ os.makedirs( os.path.join( output_path, name ) ) for name in names if not os.path.exists( os.path.join( output_path, name ) ) ]

		for count, name in enumerate( names ):
			print count
			group = name_group[ count ]
			out_group = [ os.path.join( output_path, name, os.path.basename( i ) ) for i in group[1] ]

			def run( x, y ):
				import shutil
				return shutil.move( x, y )

			pool = mp.Pool( 15 )
			out = pool.map( lambda x: run(x[0], x[1]), zip( group[1], out_group ) )
			pool.close()

