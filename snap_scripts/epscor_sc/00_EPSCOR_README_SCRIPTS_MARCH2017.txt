2017 NEW:
# HOW TO RUN EPSCoR Southcentral Test Case -- Data Creation
scripts_dir = '/workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_sc'
os.chdir( scripts_dir )

#### STEP 1: PREP THE RAW CMIP5 DATA -- STANDARDIZATION
	`ipython prep_raw_epscor_sc.py`

#### STEP 2: DOWNSCALE THE PREPPED CMIP5 DATA
	`ipython run_downscale_allslurm_epscor_sc.py`
	# ** wait for above to complete run then follow on to min/max **
	`ipython run_downscale_allslurm_epscor_sc_minmax.py`
	`ipython move_CCSM4_files.py`
	`ipython move_anomalies_to_common_dir.py` # move the interpolated anomalies

#### STEP3: COMPUTE 5ModelAvg
	`ipython wrap_slurm_compute_5ModelAvg_epscor_sc.py`

#### STEP4: DERIVED GRIDS
	`ipython wrap_slurm_run_post_downscaling_epscor_sc.py`
	# change the name of the CRU folder
	`cd /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled`
	`mv ts323 CRU_TS323`

#### STEP5: CROP / CLIP to EPSCoR Southcentral Test Case DOMAIN
	`ipython crop_clip_to_epcor_sc.py` # run with changing pathing from downscaled to derived_grids

#### STEP 6: RUN TABULAR CALCULATIONS OVER CLIPPED FILES
	`cd /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_MARCH2017/downscaled`
	`mv CRU_TS323 ts323`
	`ipython run_all_tabular_outputs_slurm_epscor_sc.py`
	`cd /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_MARCH2017/downscaled`
	`mv ts323 CRU_TS323`


##### [hack] STEP: CHANGE cru_ts323 to CRU_TS323
	```python
	# modify CRU names 
	import os
	from pathos.mp_map import mp_map

	paths = [ '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled', '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids', '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_MARCH2017/derived/grids', '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/EPSCOR_SC_DELIVERY_MARCH2017' ]

	# paths = [ '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/CRU_TS323' ]

	for path in paths:
		filelist = []
		root_list = []
		for root, subs, files in os.walk( path ):
			if root.endswith( 'ts323' ):
				root_list = root_list + [ root ]
			if files > 0:
				filelist = filelist + [ os.path.join( root, fn ) for fn in files if '_cru_' in fn or '_323' in fn ]

		# rename the files 
		_ = mp_map( lambda fn: os.rename( fn, fn.replace('_cru_', '_CRU_').replace( '_cru_ts323_', '_CRU_TS323_' ).replace( '_ts323_', '_TS_323_') ), filelist, ncpus=32 )

		# rename the dirs
		[ os.rename( root, root.replace('/ts323', '/CRU_TS323') ) for root in root_list ]
	```
