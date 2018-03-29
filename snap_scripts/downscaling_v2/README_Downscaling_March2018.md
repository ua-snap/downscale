# How I am running downscaling -- March 2018
---
#### PROCESSING STEPS:
1. run wrap_downscaler_cmip5_slurm.py
	- [tas, pr] cmip5 akcan 2km
2. run wrap_downscaler_cmip5_slurm_tem.py
	- [hur, clt] cmip5 akcan 2km
3. run wrap_downscaler_cru_ts40_slurm.py
	- [tas, pr] cru akcan 2km
4. run wrap_downscaler_cru_ts40_slurm_tem.py
	- [hur, clt] cru akcan 2km
5. run run wrap_downscaler_cmip5_slurm_minmax.py
	- [tasmin,tasmax] cmip5 akcan 2km
6. run wrap_downscaler_cru_ts40_slurm_minmax.py
	- [tasmin,tasmax] cru_ts40 akcan 2km
7. run generate_vap_hires_tem_iem.py
	- [vap] cmip5 / cru_ts40 akcan 2km
8. UPDATE CCSM4 to NCAR-CCSM4 folder titles... (script needs writing here)
	- change folder names to the proper modelname which is not served properly from PCMDI
9. run wrap_slurm_compute_5ModelAvg_cmip5.py
	- [tas,tasmin,tasmax,pr,vap,hur] 5ModelAvg cmip5 akcan 2km
