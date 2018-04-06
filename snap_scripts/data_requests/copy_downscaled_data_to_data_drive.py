# copy data from project_data to /Data/Base_Data 
def copy_file( x ):
	in_fn, out_fn = x
	dirname, basename = os.path.split( out_fn )
	try:
		if not os.path.exists( dirname ):
			_ = os.makedirs( dirname )
	except:
		pass
	return shutil.copy2( in_fn, out_fn )

if __name__ == '__main__':
	import os, shutil
	import multiprocessing as mp

	input_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	# input_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_10min'
	# input_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/alfresco_1km'
	output_dir = '/Data/Base_Data/Climate/AK_CAN_2km_v2_1'
	# output_dir = '/Data/Base_Data/Climate/AK_CAN_10min_v2_1'
	# output_dir = '/Data/Base_Data/ALFRESCO/ALFRESCO_Master_Dataset_v2_1/ALFRESCO_Model_Input_Datasets/AK_CAN_Inputs/Climate'
	# models = [ 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R']
	models = [ 'GFDL-CM3', 'NCAR-CCSM4', '5ModelAvg', 'CRU-TS40' ]

	# make args
	files = [ os.path.join(r,fn) for r,s,files in os.walk(input_dir) for fn in files if fn.endswith('.tif') and fn.split('_')[-4] in models and 'slurm' not in fn and 'rsds' in fn ]
	args = [ (fn, fn.replace(input_dir, output_dir)) for fn in files ]

	pool = mp.Pool( 64 )
	done = pool.map( copy_file, args )
	pool.close()
	pool.join()
