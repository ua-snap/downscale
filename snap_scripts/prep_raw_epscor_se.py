# prep cmip5 data downloaded using SYNDA for EPSCoR SE project
# author: Michael Lindgren -- June 09, 2016

import pandas as pd
import os

class Files( object ):
	def __init__( self, base_dir, *args, **kwargs ):
		'''
		list the files from the nested directory structure
		generated by SYNDA application and access of the ESGF CMIP5 data holdings
		'''
		self.base_dir = base_dir
		self.files = self.list_files( )
		self.df = self._to_dataframe( )

	def list_files( self ):
		return [ os.path.join( root, fn ) for root, subs, files in os.walk( self.base_dir ) \
					if len( files ) > 0 for fn in files if fn.endswith( '.nc' ) ]
	@staticmethod
	def _split_fn( fn ):
		return os.path.basename( fn ).split( '.' )[0].split( '_' )
	@staticmethod
	def f( x ):
		'''
		take the files dataframe and split the years
		into begin year/month and end year/month and 
		add new columns to a new dataframe
		'''
		begin, end = x[ 'years' ].split( '-' )
		x['begin_month'] = begin[4:] 
		x['begin_year'] = begin[:4]
		x['end_month'] = end[4:]
		x['end_year'] = end[:4]
		return x
	def _to_dataframe( self ):
		import pandas as pd
		out = []
		for fn in self.files:
			variable, cmor_table, model, scenario, experiment, years = self._split_fn( fn )
			out.append( {'fn':fn, 'variable':variable, 'cmor_table':cmor_table, \
						'model':model, 'scenario':scenario, 'experiment':experiment, 'years':years } )
			column_order = ['fn', 'variable', 'cmor_table', 'model', 'scenario', 'experiment', 'years']
		return pd.DataFrame( out, columns=column_order ).apply( self.f, axis=1 )


# PREP THE INPUT NETCDF FILES FROM AR5
import itertools, glob, os
from downscale import preprocess

# some setup args
base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/raw_cmip5'
prepped_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prepped_cmip5'
variables = [ 'tasmax', 'tasmin' ]
scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'GFDL-CM3', 'CCSM4' ]

if not os.path.exists( prepped_dir ):
	os.makedirs( prepped_dir )

# lets get the darn data returned that we want:
files_df = Files( base_dir )._to_dataframe( )
log = open( os.path.join( prepped_dir, 'log_file_prep.txt'), 'w' )
for variable, model, scenario in itertools.product( variables, models, scenarios ):
	# get the files we want to work with for this run
	cur_files = files_df[ (files_df.variable == variable) & (files_df.model == model) & (files_df.scenario == scenario) ]['fn'].tolist()
	if 'historical' in scenario:
		years = (1900,2005)
	else:
		years = (2006,2100)

	raw_path = os.path.dirname( cur_files[0] )
	output_path = os.path.join( prepped_dir, model, scenario, variable )
	
	if not os.path.exists( output_path ):
		os.makedirs( output_path )

	experiment = 'r1i1p1'
	try:
		pp = preprocess.Preprocess( raw_path, variable, model, scenario, experiment, years )
		pp.write_nc( output_path, True )
	except:
		log.write( 'error : %s - %s - %s - %s - %s - %s \n\n' % (raw_path, variable, model, scenario, experiment, years)  )
		pass

# close the log file
log.flush()
log.close()
