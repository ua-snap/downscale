# THIS IS THE WAY TO RUN EPSCOR_SE DOWNSCALING:
-------------------------------------------------------

# DOWNSCALE PROCEDURE
Run this script
	downscale_epscor_se_distributedrun.py
with this:
	wrap_downscaler_slurmrun_epscor_se.py

# CRU
Run this script:
	run_cru_snap_epscor_se_CLI.py
with this script:
	wrap_slurm_run_cru_snap_epscor_se_CLI.py

# change the model name of CCSM4 folder to NCAR-CCSM4
	mv CCSM4/ NCAR-CCSM4

# hackily fix the single digit months -- runs across both data sets -- this needs a proper solution.
	fix_singledigit_months_espcor_se.py

# hackily fix the issue that cru vars are called different names than cmip5
	cru_naming_fix_epscor_se.py	

# 5ModelAvg:
	compute_5ModelAvg_epscor_se.py

# 5ModelAvg was not computed for the historical data matt produced.  
#  I also had to move some stuff into the proper folder structure and remove some
#  old unused models to an `other_models` folder.
	compute_5ModelAvg_cmip5_snap_pr_tas_epscor_se.py

# symlink files from /Data for ease-of-processing
	symlink_pr_tas_epscor_se.py

# post-downscaling
	# # move cru into a folder called ts323 its 'model' and 'scenario' name
	mkdir ts323	
	mv * ts323/
	cd ts323/
	mkdir historical
	mv * historical
	wrap_slurm_run_post_downscaling_epscor_se.py

# SWI calculation
	swi_calc_decadal_epscor_se.py

# clip to EPSCOR_SE extent
	crop_clip_to_epscor_se.py

# derived tabular
	run_all_tabular_outputs_slurm_epscor_se.py

# symlink into the DELIVERY folder -- for zipping:
# I ALSO CHANGED THE NAMES OF THE ts323 folder to CRU_TS323!
	screen cp -rs /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_epscor_se/* /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_sym/derived/grids

	screen cp -rs /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_tabular/*.csv /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_sym/derived/tabular/csv

	screen cp -rs /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_tabular/*.json /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_sym/derived/tabular/json

	screen cp -rs /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cmip5_epscor_se/* /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_sym/downscaled/

	screen cp -rs /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_cru_epscor_se/* /workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/DELIVERY_JULY2016_sym/downscaled/

# remove any slurm_files folders from the DELIVERY_JULY2016_sym folder

