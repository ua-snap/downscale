# THIS IS THE WAY TO RUN EPSCOR_SE DOWNSCALING:
-------------------------------------------------------

# DOWNSCALE PROCEDURE -- this command fires off all jobs for cmip5 and cru on atlas
	ipython /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_se/run_downscale_allslurm_epscor_se.py

# CALCULATE 5ModelAvg -- tasmin/tasmax variables
	mv CCSM4/ NCAR-CCSM4 # this didnt get the right name for some reason
	ipython /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_se/compute_5ModelAvg_epscor_se.py

# CALCULATE A 5ModelAvg for the historical tas/pr since it does not exist
#  I also had to move some stuff into the proper folder structure and remove some
#  old unused models to an `other_models` folder.
#  you're welcome.
	compute_5ModelAvg_cmip5_snap_pr_tas_epscor_se.py

# symlink files from /Data for ease-of-processing
	symlink_pr_tas_epscor_se.py

# RUN DERIVED GRIDS
	# run for cmip5 and cru individually.
	# -- this runs all derived grids including swi
	ipython /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_se/wrap_slurm_run_post_downscaling_epscor_se.py 

# CROP CLIP! run for the downscaleds and derived_grids
# clip to EPSCOR_SE extent
	crop_clip_to_epscor_se.py

# RUN DERIVED TABULAR
	ipython /workspace/UA/malindgren/repos/downscale/snap_scripts/epscor_se/run_all_tabular_outputs_slurm_epscor_se.py

# I ALSO CHANGED THE NAMES OF THE ts323 folder to CRU_TS323!

# remove any slurm_files folders from the DELIVERY_JULY2016 folder

# ALSO Dont FORGET we clipped cropped the DOF, DOT, LOGS data for angie so it was consistent.
	# clip em
	crop_clip_to_epscor_se_derived_other.py
	# move em to the final directory
	move_clipped_dof_dot_log_epscor_se.py

# ALSO I fixed the precision issues in the output CSV files using this script...  Also added to 
# full code, but this is here for legacy viewing.
	fix_derived_tabular_precision_epscor_se.py

# # # # # # # # # # ADDED TO wrap_slurm_run_post_downscaling_epscor_se.py  # # # # # # # # # #
# calc seasonal decadals from seasonal annuals created above
#	decadal.py
#
# SWI calculation
#	swi_calc_decadal_epscor_se.py
# # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


