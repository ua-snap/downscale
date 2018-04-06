[![Stories in Ready](https://badge.waffle.io/ua-snap/downscale.svg?label=ready&title=Ready)](http://waffle.io/ua-snap/downscale)
[![Build Status](https://travis-ci.org/ua-snap/downscale.svg?branch=master)](https://travis-ci.org/ua-snap/downscale)

### __[ DOCUMENTATION UPDATE IN PROGRESS APRIL 2018 ]__


downscale
---------

Python Package for Simple Delta-Downscaling of Climatic Research Unit's TS 3.x OR PCMDI's CMIP5 data to a set of baseline monthly climatologies. Focus of this work is for regional downscaling of data over Alaska region, but we do aim to make it flexible enough to do this very simple downscaling over any area of interest.

#### Note:
This package developing rapidly and is expected to be somewhat problematic until full release. It is also geared specifically
for the needs of Scenarios Network for Alaska + Arctic Planning (SNAP).

##### AR5 Example:
```python
import glob, os
import downscale

# SETUP BASELINE
clim_path = './climatology'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# SETUP DATASET
future_fn = './hur_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601_210012.nc'
historical_fn = './hur_Amon_IPSL-CM5A-LR_historical_r1i1p1_185001_200512.nc'
variable = 'hur'
model = 'IPSL-CM5A-LR'
scenario = 'rcp26'
historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None )
future = downscale.Dataset( future_fn, variable, model, scenario, units=None )

# DOWNSCALE
output_dir = './outputs'
clim_begin = '1961'
clim_end = '1990'
ar5 = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, future, \
		metric='mean', ds_type='absolute', level=1000, level_name='plev' )
ar5.downscale( output_dir=output_dir )
```

##### CRU Example:
```python
import glob, os
import downscale

# SETUP BASELINE
clim_path = './climatology'
filelist = glob.glob( os.path.join( clim_path, '*.tif' ) )
baseline = downscale.Baseline( filelist )

# SETUP DATASET
historical_fn = './cru_ts3.23.1901.2014.cld.dat.nc'
variable = 'cld'
model = 'cru_ts31'
scenario = 'observed'

# this read in will interpolate across NAs to fill in around the coastlines for later masking
historical = downscale.Dataset( historical_fn, variable, model, scenario, units=None, interp=True )

# DOWNSCALE
output_dir = './outputs'
clim_begin = '1961'
clim_end = '1990'
cru = downscale.DeltaDownscale( baseline, clim_begin, clim_end, historical, \
						metric='mean', ds_type='relative' )
cru.downscale( output_dir=output_dir )

```
