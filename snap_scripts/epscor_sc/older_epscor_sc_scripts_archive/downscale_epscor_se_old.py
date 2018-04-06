# DOWNSCALE EPSCoR SE data tasmax/tasmin
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
import glob, os, rasterio
import downscale

# SETUP BASELINE
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/tmax'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# SETUP DATASET
output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/test'
future_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data/IPSL-CM5A-LR/hur/hur_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
historical_fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/snap_prepped_data/IPSL-CM5A-LR/hur/hur_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
variable = 'hur'
model = 'IPSL-CM5A-LR'
scenario = 'rcp26'
historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None )
future = downscale.Dataset( future_fn, variable, model, scenario, units=None )

# DOWNSCALE
mask = rasterio.open(baseline.filelist[0]).read_masks(1)
clim_begin = '1961'
clim_end = '1990'
ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
		metric='mean', ds_type='absolute', level=1000, level_name='plev' )# add in the mask!
ar5.downscale( output_dir=output_dir )

# CRU historical
output_dir = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/test'
historical_fn = '/Data/Base_Data/Climate/World/CRU_grids/CRU_TS323/cru_ts3.23.1901.2014.cld.dat.nc'
historical = downscale.Dataset( historical_fn, 'cld', 'cru_ts31', 'observed', units=None, interp=True )
clim_begin = '1961'
clim_end = '1990'
cru = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, metric='mean', ds_type='relative' )
cru.downscale( output_dir=output_dir )




	# filelist = [ os.path.join(root, fn) for root, subs, files in os.walk(base_dir) if len(files) > 0 for fn in files if variable in fn and scenario in fn and model in fn ]
	if len( filelist ) > 0:
		if len(filelist) > 0:
			print( filelist[0] )

			if 'historical' in scenario:
				years = (1850,2005)
			else:
				years = (2006,2100)
			raw_path = os.path.dirname( filelist[0] )
			output_path = os.path.join( '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/prepped', model, scenario, variable )
			
			if not os.path.exists( output_path ):
				os.makedirs( output_path )

			experiment = 'r1i1p1'
		
			pp = preprocess.Preprocess( raw_path, variable, model, scenario, experiment, years )
			pp.write_nc( output_path, True )
		else:
			print( 'NO FILES FOR THIS ONE!' )

		# PARSE OUT THE FILES THAT OVERLAP CHRONOLOGICALLY HERE - this is an issue fo the series being overlapping for the 	
		# same series 2005-2300 or 2006-2100...

		# parse out the filelist years:
		filelist = sorted( glob.glob( os.path.join( '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/docker_move/synda/sdt/data/cmip5/output1/NOAA-GFDL/GFDL-CM3/rcp26/mon/atmos/Amon/r1i1p1/v20120227', 'tasmax*.nc' ) ) )
		fileyears = [ i.split( '.' )[0].split( '_' )[-1].split( '-' ) for i in filelist ]
		fileyears = [ ( int(i), int(j) ) for i,j in fileyears ]

		

		try:
		except:
			print( 'BROKEN! : ' + ' '.join([variable, model, scenario, experiment, str(years[0]), str(years[1]) ]) )


# # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # #
# preprocess example -- EPSCOR Downscaling.
from downscale import preprocess

# some setup
path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/docker_move/synda/sdt/data/cmip5/output1/IPSL/IPSL-CM5A-LR/historical/mon/atmos/Amon/r1i1p1/v20110406'
output_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/prepped'
variable = 'tasmax'
scenario = 'historical'
model = 'IPSL-CM5A-LR'
experiment = 'r1i1p1'
years = (1850, 2005)

pp = preprocess.Preprocess( path, variable, model, scenario, experiment, years )
pp.write_nc( output_path, True )

# # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # ## # # # # #

import glob, os, rasterio
import downscale
from downscale import preprocess

# SETUP BASELINE
clim_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/prism/tmax'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# PREP FUTURE
path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/tasmax_test'
variable = 'tasmax'
model = 'MRI-CGCM3'
scenario = 'rcp26'
experiment = 'r1i1p1'
years = (2006,2100)

pp_test = preprocess.Preprocess( path, variable, model, scenario, experiment, years )
pp_test.write_nc( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data' )

# PREP HISTORICAL
# PREP FUTURE
path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/tasmax_test'
variable = 'tasmax'
model = 'MRI-CGCM3'
scenario = 'historical'
experiment = 'r1i1p1'
years = (1850,2005)

pp_test = preprocess.Preprocess( path, variable, model, scenario, experiment, years )
pp_test.write_nc( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data' )

# NEW PREPPED DATA FOLDER:
input_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/tasmax_test/prepped'

# SETUP DATASET
output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/tasmax_test/output_test'
future_fn = os.path.join( input_path, 'tasmax_Amon_MRI-CGCM3_rcp26_r1i1p1_200601-210012.nc' )
historical_fn = os.path.join( input_path, 'tasmax_Amon_MRI-CGCM3_historical_r1i1p1_185001-200512.nc' )
variable = 'tasmax'
model = 'MRI-CGCM3'
scenario = 'rcp26'
historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None )
future = downscale.Dataset( future_fn, variable, model, scenario, units=None )

# DOWNSCALE
mask = rasterio.open( baseline.filelist[0] ).read_masks( 1 )
clim_begin = '1961'
clim_end = '1990'

ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
		metric='mean', downscaling_operation='add', mask=mask, mask_value=0, ncpus=32, \
		src_crs={'init':'epsg:4326'}, src_nodata=-9999.0, dst_nodata=None,
		post_downscale_function=None ) # add in the mask!
ar5.downscale( output_dir=output_dir )


# /workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/tasmin_tasmax_access/docker_move/synda/sdt/data/cmip5/output1/NASA-GISS/GISS-E2-R/rcp85/mon/atmos/Amon/r1i1p1/v20121016/tasmax_Amon_GISS-E2-R_rcp85_r1i1p1_200601-202512.nc
