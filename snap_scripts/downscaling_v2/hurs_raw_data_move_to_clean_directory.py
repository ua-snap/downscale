# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# hurs raw files check and move before prepping and downscaling (next)...
# Michael Lindgren (malindgren@alaska.edu) September 2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import os, shutil
import pandas as pd

os.chdir('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data/cmip5')

files = [ os.path.join(r,fn) for r,s,files in os.walk('./') for fn in files if fn.endswith('.nc')]
split_files = [ os.path.basename(fn).split('_') for fn in files ]
dates = [fn.pop(-1).split('.')[0].split('-') for fn in split_files]
columns = ['variable','cmor_table','model', 'scenario', 'experiment']
df = pd.DataFrame(split_files, columns=columns)
df2 = pd.DataFrame(dates, columns=['begin','end'])
df['begin'] = df2['begin'].astype(int)
df['end'] = df2['end'].astype(int)
df['fn'] = files

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # copy the files we want to a new clean standardized directory structure
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def remove_files( df ):
	return df[(df.end <= 210012) & -((df.begin == 200501) & (df.end == 210012) )]['fn']

out = df.groupby(['model','scenario']).apply( remove_files )
new_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted'
variable = 'hurs'

for name, dat in out.groupby(['model','scenario']):
	model, scenario = name
	filelist = dat.tolist()
	output_dir = os.path.join(new_dir, model, scenario, variable)
	if not os.path.exists( output_dir ):
		_ = os.makedirs( output_dir )

	done = [ shutil.copy( fn, os.path.join( output_dir, os.path.basename(fn)) ) for fn in filelist ]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# fix the missing IPSL -- DUE TO ENDING IN 2300 oddly in only *some* series...
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
os.system('mkdir -p /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp26/hurs/')
os.system('cp /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp26/mon/atmos/Amon/r1i1p1/v20120114/hurs/hurs_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601-230012.nc /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp26/hurs/')
os.system('mkdir -p /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp45/hurs/')
os.system('cp /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp45/mon/atmos/Amon/r1i1p1/v20110914/hurs/hurs_Amon_IPSL-CM5A-LR_rcp45_r1i1p1_200601-230012.nc /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp45/hurs/')
os.system('mkdir -p /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp85/hurs/')
os.system('cp /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp85/mon/atmos/Amon/r1i1p1/v20111103/hurs/hurs_Amon_IPSL-CM5A-LR_rcp85_r1i1p1_200601-230012.nc /workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/cmip5/raw_download_rhs/data_raw_sorted/IPSL-CM5A-LR/rcp85/hurs/')
