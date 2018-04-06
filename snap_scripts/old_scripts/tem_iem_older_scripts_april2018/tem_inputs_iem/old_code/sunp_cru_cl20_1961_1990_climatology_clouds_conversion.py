# script to convert the Sunshine variable from the CRU CL 2.0 Data to clouds
import rasterio, glob, os
import numpy as np

# list the files
input_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts20/sunp/akcan'
output_path = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/TEM_Data/cru_october_final/cru_ts20/cld/akcan'

if not os.path.exists( output_path ):
	os.makedirs( output_path )

l = sorted( glob.glob( os.path.join( input_path, '*.tif' ) ) )

# get metadata from the first in the input list:
meta = rasterio.open( l[0] ).meta

# stack the rasters with time as the first dimension
sunp = np.rollaxis( np.dstack([ rasterio.open( fn ).read( 1 ) for fn in l ]), -1 )

# mask it 
sunp = np.ma.masked_where( sunp < 0, sunp )

# calculate the clouds from the sunshine
cld = 100 - sunp

# and now write out the layers to GTiff
cld_split = [ cld[i, ...] for i in range( cld.shape[0] ) ]

for in_fn, arr in zip( l, cld_split ):
	output_filename = os.path.join( output_path, os.path.basename(in_fn).replace( 'sunp', 'cld_from_sunp' ) )
	with rasterio.open( output_filename, 'w', **meta ) as out:
		out.write( arr, 1 )


