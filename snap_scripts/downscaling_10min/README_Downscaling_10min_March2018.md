# __Delta Downscaling Far-Futures 10'__ -- SNAP
#### __CMIP5 / CRU-TS40__

###### Michael Lindgren (malindgren@alaska.edu) -- March/April 2018
---
##### PROCESSING STEPS -- AKCAN 2km (run in order):
1. run `wrap_downscaler_cmip5_slurm_far-futures.py`
	- [tas,pr] cmip5 akcan 15km
2. run `wrap_downscaler_cru_ts40_slurm.py`
	- [tas,pr] cru ts40 15km
	<--[- wait for above to complete -]-->
3. run `wrap_slurm_compute_5ModelAvg_cmip5_far-futures.py`
	- [tas,pr] 5ModelAvg cmip5 15km
4. run `move_anomalies_to_common_dir_10min.py`
	- move anomalies to common dir in `project_data`
