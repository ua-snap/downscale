### STEPS TO RUNNING CMIP5 / CRU DOWNSCALING
#### March 18, 2017 -- Michael Lindgren (malindgren@alaska.edu)

1. change working directory to the scripts folder
	```python
	scripts_dir = '/workspace/UA/malindgren/repos/downscale/snap_scripts/downscaling_v2'

	os.chdir( scripts_dir )
	```
2. run all non-min/max variables...
	`ipython run_downscale_wrappers.py`

3. [wait for previous process to be fully complete]:
	run tas min/max variables across all groups...
	`ipython run_downscale_wrappers_minmax.py`

4. [wait for previous process to be fully complete]:
	copy mislabeled directory tree to existing NCAR_CCSM4
	recursively, which adds each to the proper leaf of the 
	tree in the receiving set.
	`ipython move_CCSM4_files.py`
	*** MOVE THE ANOMALIES TO AN ANOMALIES FOLDER NOW....
	`ipython move_anomalies_to_common_dir.py`

5. compute 5ModelAvg across all variables / models
	`ipython wrap_slurm_compute_5ModelAvg_cmip5.py`

6. run standardization to *new* SNAP standard of crs, and oob
	`ipython wrap_slurm_standardize_data_snap.py`

