def list_all_data(path):
	return [ os.path.join(root, fn) for root, subdirs, files in os.walk(path) \
			if len(files) > 0 for fn in files ]

def filter_files(files, variables, models):
	out_files = []
	for variable, model in itertools.product(variables,models):
		out_files = out_files + [fn for fn in files if variable in fn and model in fn]
	return out_files

def make_args(input_base_path, output_base_path, variables, models):
	files = list_all_data(input_base_path)
	files = filter_files(files, variables, models)
	return [ (fn, fn.replace(input_base_path, output_base_path)) for fn in files ]

def copy_file(input_filename, output_filename):
	_ = subprocess.call(['cp', '-p', input_filename, output_filename])
	return output_filename

def run_copy(x):
	input_filename, output_filename = x
	dirname = os.path.dirname(output_filename)
	try:
		if not os.path.exists(dirname):
			_ = os.makedirs(dirname)
	except:
		pass
	return copy_file(input_filename, output_filename)


if __name__ == '__main__':
	# copy all delta downscaled data to /Data/Base_Data/Climate as a "truth" repository
	import os, glob, subprocess, itertools
	import pandas as pd
	import multiprocessing as mp

	# pathing
	input_base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled'
	output_base_path = '/Data/Base_Data/Climate/AK_CAN_2km_v2_2'

	# global args
	variables = ['tas', 'pr', 'clt', 'tasmin', 'tasmax', 'rsds', 'vap', 'hurs',]
	future_models = ['GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4', '5ModelAvg',]
	observed_models = ['CRU-TS40',]
	models = future_models + observed_models

	# list the data:
	args = make_args(input_base_path, output_base_path, variables, models)

	# run in parallel
	pool = mp.Pool(20)
	pool.map(run_copy, args)
	pool.close()
	pool.join()

