# # # # # # # # # # # # # # # # # # #
# RUN EXAMPLES -- DOWNSCALE -- SNAP
# # # # # # # # # # # # # # # # # # #

# preprocess downloaded files from SYNDA
path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/docker_move/synda/sdt/data/cmip5/output1/NASA-GISS/GISS-E2-R/rcp26/mon/atmos/Amon/r1i1p1/v20130731'
variable = 'tasmax'
model = 'GISS-E2-R'
scenario = 'rcp26'
experiment = 'r1i1p1'
years = (2006,2100)

# run it
pre = Preprocess( path, variable, model, scenario, experiment, years )
pre.write_nc( output_path='/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/prepped' ) #to disk




# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/docker_move/synda/sdt/data/cmip5/output1/NOAA-GFDL/GFDL-CM3/rcp26/mon/atmos/Amon/r1i1p1/v20120227'
# model = 'GFDL-CM3'

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# def cmip5_parse_fn( fn, elems=['variable', 'cmor_table', 'model', 'scenario', 'experiment', 'years'] ):
# 	'''
# 	naming convention anchored name splitting
# 	'''
# 	base,ext = os.path.splitext( os.path.basename( fn ) )
# 	name_split = base.split( '_' )
# 	return dict( zip( elems, name_split ) )
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *