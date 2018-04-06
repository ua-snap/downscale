# # # # # # # # # # # # # # # #
# # CHANGE THE FILENAMES OF THE CRU DATA TO MATCH TO THE VARIABLE NAMES THEY WILL BE DOWNSCALED _TO_:
# # # # # # # # # # # # # # # #
import os, glob, shutil

base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tem_data_sep2016/cru/cru_ts20'
variables = ['reh','tmp','pre']
new_names = ['hur', 'tas', 'pr']

for variable, new_name in zip(variables, new_names):
	files = glob.glob( os.path.join( base_path, variable, '*'+variable+'*.tif' ) )

	for fn in files:
		print( fn )
		dirname, basename = os.path.split( fn )
		new_dir = dirname.replace( variable, new_name )
		if not os.path.exists( new_dir ):
			os.makedirs( new_dir )
		
		# fix single digit months
		mon = basename.split('_')[-2]
		basename = basename.replace( '_'+mon+'_', '_'+months_lookup[ int(mon) ]+'_' )
		
		# change variable naming
		new_fn = os.path.join( new_dir, basename.replace( variable, new_name ) )
		shutil.copy( fn, new_fn )
