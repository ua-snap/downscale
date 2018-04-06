# zip some of the SNAP-held raw CMIP5 data for Helene to test. -- March 2018
import os, zipfile, zlib

base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/prepped'
output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/helene_zipped_raw_prepped'
variable = 'hur_'
models = ['IPSL-CM5A-LR', 'GFDL-CM3', 'CCSM4', 'MRI-CGCM3', 'GISS-E2-R']
# scenarios = ['historical','rcp45','rcp60','rcp85']

files = [ os.path.join(r,fn) for r,s,files in os.walk(base_path) for fn in files if fn.endswith('.nc') and variable in fn ]

for model in models:
	print( 'zipping: {}'.format(model) )
	filelist = [ fn for fn in files if model in fn ]

	if model == 'CCSM4':
		model = 'NCAR-CCSM4'

	out_fn = os.path.join( output_path, 'hur_{}_cmip5_raw.zip'.format(model) )

	# open a new zipfile
	zf = zipfile.ZipFile( out_fn, mode='w' )
		
	for fn in filelist:
		zf.write( fn, arcname=os.path.basename(fn), compress_type=zlib.DEFLATED )

	# cleanup
	zf.close()
	zf = None
	del zf
