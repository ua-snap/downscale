# # # # # 
# MOVE THE NEWLY DOWNLOADED TAS / PR CMIP5 data from work desktop to /Shared
# # # # # 

def move_new_dir( fn, output_dir ):
	dirname, basename = os.path.split( fn )
	elems = basename.split('.')[0].split( '_' )
	variable, cmor_table, model, scenario, experiment, years = elems
	new_dir = os.path.join( output_dir, model, scenario, variable )
	try:
		if not os.path.exists( new_dir ):
			os.makedirs( new_dir )
	except:
		pass
	return shutil.copy( fn, new_dir )

if __name__ == '__main__':
	import os, glob, shutil
	path = '/srv/synda/sdt/data'
	output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/raw_cmip5_tas_pr'

	filelist = []
	for root, subs, files in os.walk( path ):
		if len( files ) > 0:
			filelist = filelist + [ os.path.join( root, i ) for i in files if i.endswith( '.nc' ) ]

	out = [ move_new_dir( fn, output_dir ) for fn in filelist ]


# # # # # # # # 
# # CHECK FOR DUPLICATES and remove by hand. this is tedious.
# GFDL - OK
# CCSM4 - FIXED OK
# GISS-E2-R - OK
# IPSL - OK
# MRI - OK