### STEPS TO RUNNING CMIP5 / CRU DOWNSCALING
#### March 18, 2017 -- Michael Lindgren (malindgren@alaska.edu)

1. change working directory to the scripts folder
	```python
	scripts_dir = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_10min'

	os.chdir( scripts_dir )
	```
2. run all tas/pr ...
	`ipython wrap_downscaled_cmip_slurm.py`
	`ipython wrap_downscaled_cru_slurm.py`
	`ipython wrap_downscaled_cru40_slurm.py`

# SINCE THERE IS NO MIN/MAX TEMP DATA FOR CL20 --> 10min 
# THERE IS NO MIN/MAX RUN FOR THIS GROUP.

4. [wait for previous process to be fully complete]:
	# rename CCSM4 to the proper NCAR-CCSM4
	mv CCSM4 NCAR-CCSM4

5. compute 5ModelAvg across all variables / models
	`ipython wrap_slurm_compute_5ModelAvg_cmip5.py`




# NOT RUN YET...
6. run standardization to *new* SNAP standard of crs, and oob
	`ipython wrap_slurm_standardize_data_snap.py`

