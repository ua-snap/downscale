# compare the DOT/DOF/LOGS with Matts::

import rasterio, os, glob
import numpy as np

new = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/derived_v2_parallel/GFDL-CM3/rcp85/dof'
old = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/test/dof_dot_logs_test_matt/dof/decadal_mean'

new = sorted(glob.glob(os.path.join(new,'*.tif')))
old = sorted(glob.glob(os.path.join(old,'*.tif')))

# drop not full decades...
new.pop(0); new.pop(-1)

new = np.array([rasterio.open(i).read(1) for i in new])
old = np.array([rasterio.open(i).read(1) for i in old])

# update the oob in the old series
old[old < -1] = -9999

diff = new - old

# /Volumes/Shared/Tech_Projects/DeltaDownscaling/project_data/derived_v2_parallel/GFDL-CM3/rcp85/dof/dof_GFDL-CM3_rcp85_2010s.tif
# /Volumes/Shared/Tech_Projects/DeltaDownscaling/project_data/test/dof_dot_logs_test_matt/dof/decadal_mean/dof_GFDL-CM3_rcp85_2010_2019.tif


rst_new = rasterio.open('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/derived_v2_parallel/GFDL-CM3/rcp85/dof/dof_GFDL-CM3_rcp85_2010s.tif').read(1)
rst_old = rasterio.open('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/test/dof_dot_logs_test_matt/dof/decadal_mean/dof_GFDL-CM3_rcp85_2010_2019.tif').read(1)