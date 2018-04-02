# # # # # # # # # # # # # # # # # # # # # # # # # 
# # MOVE ALL _RAW_ DATA TO COMMON STRUCTURED DIR
# # # # # # # # # # # # # # # # # # # # # # # # # 

def move_fn( fn, output_path ):
	import shutil, os

	dirname, basename = os.path.split( fn )
	variable, cmor_table, model, scenario, experiment, years = basename.split( '_' )
	try:
		new_path = os.path.join( output_path, model, scenario, variable )
		if not os.path.exists( new_path ):
			os.makedirs( new_path )
	except:
		pass

	return shutil.copy( fn, os.path.join( new_path, basename ) )

if __name__ == '__main__':
	import os

	base_dirs = ['/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/data/cmip5']
	output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5_nwt/cmip5_raw_restructure'

	for base_dir in base_dirs:
		files_out = []
		print( base_dir )
		files = [ os.path.join(r, fn) for r,s,files in os.walk(base_dir) for fn in files if fn.endswith( '.nc' )]
		files_out = files_out + [ fn for fn in files if 'tas_' in fn or 'pr_' in fn ] # only want tas and pr files -- NWT
		done = [ move_fn(fn, output_path) for fn in files_out ]
		

# esg-dn1.nsc.liu.se,esg.pik-potsdam.de,esgdata.gfdl.noaa.gov,esgf-data.dkrz.de,esgf-index1.ceda.ac.uk,esgf-node.ipsl.upmc.fr,esgf-node.jpl.nasa.gov,esgf-node.llnl.gov,esgf.esrl.noaa.gov,esgf.nci.org.au