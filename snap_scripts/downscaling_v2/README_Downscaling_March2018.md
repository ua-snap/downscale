# __Delta Downscaling__ -- SNAP
#### __CMIP5 / CRU-TS40__

###### Michael Lindgren (malindgren@alaska.edu) -- March/April 2018
---
##### PROCESSING STEPS -- AKCAN 2km (run in order):
1. run `wrap_downscaler_cmip5_slurm.py`
	- [`tas`, `pr`] cmip5 akcan 2km
2. run `wrap_downscaler_cmip5_slurm_tem.py`
	- [`hur`, `clt`] cmip5 akcan 2km
3. run `wrap_downscaler_cru_ts40_slurm.py`
	- [`tas`, `pr`] cru akcan 2km
4. run `wrap_downscaler_cru_ts40_slurm_tem.py`
	- [`hur`, `clt`] cru akcan 2km
5. run `wrap_downscaler_cmip5_slurm_minmax.py`
	- [`tasmin`,`tasmax`] cmip5 akcan 2km
6. run `wrap_downscaler_cru_ts40_slurm_minmax.py`
	- [`tasmin`,`tasmax`] cru_ts40 akcan 2km
7. run: `mv CCSM4/ NCAR-CCSM4/`
	- change folder names to the proper modelname which is not served properly from PCMDI
8. run: `mv ts40/ CRU-TS40/`
	- change the name of cru data folder to match something like the CMIP5 model naming
9. run `update_cru_ts40_naming.py`
	- update the cru_ts40 name to CRU-TS40 to better follow the CMIP5 naming scheme
10. run `move_anomalies_to_common_dir.py`
	- to move anomalies to a new common dir
11. run `generate_vap_hires_tem_iem.py`
	- [`vap`] cmip5 / cru_ts40 akcan 2km
12. run `wrap_slurm_compute_5ModelAvg_cmip5.py`
	- [`tas`,`tasmin`,`tasmax`,`pr`,`vap`,`hur`,`clt`] 5ModelAvg cmip5 akcan 2km
