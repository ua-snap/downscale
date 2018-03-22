#### How to compute Delta Downscaling with CMIP5 and CRU TS4.0 using PRISM or CRU baseline climatology. 
---

- rotate data to PCLL (if needed) -->> potentially do this ahead of time? would that even matter? 

- if projected:
	- read projected ds
	- read historical ds
	- compute climatology with historicals
- else:
	- read in [historical] ds
	- compute climatology with historicals

- if precip:
	- fix precip issues with climatology
	- fix precip issue with full dataset

- if tasmin/tasmax:
	- calculate deltas from tas
- else:	
	- calculate anomalies

- if cru: # land-only
	- interpolate convex hull to fill-in coastline np.nan values

- rotate data to GCLL for GDAL (via rasterio) to properly reproject (using most versions).
- reproject (REGRID) anomalies to baseline climatology extent/origin/res/crs
- add / mult reprojected anomalies to baseline climatology
- write to GeoTiff

