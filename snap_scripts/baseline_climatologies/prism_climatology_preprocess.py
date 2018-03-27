# reproject, merge, and convert to GTiff PRISM 1961-1990 Climatology
# SNAP AKCanada Domain
# March 2018
# Author: Michael Lindgren (malindgren@alaska.edu)
#

def filelister( base_dir ):
    import itertools
    # list all files
    files = [ os.path.join(r,fn) for r,s,files in os.walk(base_dir) for fn in files if fn.endswith('.asc') or fn.endswith('.txt') ]
    # filter to the groups we want.
    regions = ['alaska', 'canada']
    variables = ['ppt','tmean','tmin','tmax']
    fd = { '{}-{}'.format(region, variable):
            _sort_files_chronological(list(filter(lambda x: region in x and variable in x and '14' not in x , files))) 
                for region,variable in itertools.product(regions,variables) }
    return fd

def _sort_files_chronological( files ):
    ''' specific function hardwired to work with the raw PRISM files'''
    months = [ int(os.path.basename(fn).split('.')[0].split('_')[-1]) for fn in files ]
    df = pd.DataFrame({'months':months,'files':files})
    return df.sort_values('months')['files'].tolist()

def reproject_canada_to3338( fn, template_fn, scalevals=False, dtype='float32' ):
    '''
    reproject and scale the temp values to floats from ints 
    '''
    from rasterio.warp import calculate_default_transform, reproject, Resampling
    
    # open up the SNAP template dataset  
    with rasterio.open(template_fn) as tmp:
        kwargs = tmp.meta.copy()
        tmp_arr = tmp.read(1)
        kwargs.update({
            'driver':'MEM', 
            'dtype':dtype
        })

    src_crs = {'init':'EPSG:4326'}
    dst_crs = {'init':'EPSG:3338'}
    nodata = -9999

    with rasterio.open( fn ) as src:
        types = { 'float32':np.float32, 'int32':np.int32 }
        out_arr = np.empty_like( tmp_arr.astype( types[ dtype ] ) )
        reproject(
            source=src.read(1),
            destination=out_arr,
            src_transform=src.transform,
            src_crs=src_crs,
            src_nodata=nodata,
            dst_transform=tmp.transform,
            dst_crs=tmp.crs,
            dst_nodata=nodata,
            resampling=Resampling.nearest)

    # open an in-memory dataset in the destination crs
    dst = rasterio.open( 'reproject', 'w', **kwargs )
    if scalevals: # scale the temperature PRISM values for Canada...
        out_arr[ out_arr != dst.nodata ] = out_arr[ out_arr != dst.nodata ] / 10.0
    dst.write( out_arr, 1 )
    return dst

def open_ak_precip( fn ):
    ''' 
    deal with diffs between dtypes between reproj canada and raw ak
    since float32 is our base dtype for ALF / IEM I will stick with 
    this as the dtype to store these data, despite being integers in 
    reality.
    '''
    with rasterio.open( fn ) as rst:
        meta = rst.meta
        arr = rst.read(1).astype( np.float32 )

        meta.update({
            'driver':'MEM',
            'dtype':'float32'
            })

    out = rasterio.open( '', 'w', **meta )
    out.write( arr, 1 )
    return out


def merge_ak_canada( ak_fn, can_fn, template_fn, variable, output_path ):
    from rasterio.merge import merge
    
    # reproject canada to 3338 and rescale temperature values from ints to floats
    dtype = 'float32'
    if variable == 'ppt':
        scalevals = False
        # open alaska AAI grid -- special for precip...
        ak = open_ak_precip( ak_fn )
    else: # temperature variables
        scalevals = True
        # open alaska AAI grid
        ak = rasterio.open( ak_fn ) # already in 3338
    
    can = reproject_canada_to3338( can_fn, template_fn, scalevals=scalevals, dtype=dtype )

    # get some output bounds information for the merge...
    with rasterio.open( template_fn ) as tmp:
        bounds = tmp.bounds

    # merge 'em
    merged, transform = merge( [ak,can], bounds=bounds, res=(2000,2000), nodata=-9999 )
    # update metadata
    meta = ak.meta
    count, height, width = merged.shape
    meta.update({
        'driver':'GTiff',
        'height':height,
        'width':width,
        'crs':{'init':'EPSG:3338'},
        'transform':transform,
        'compress':'lzw',
        'dtype':dtype
        })

    # output filenaming
    varname_lookup = {'tmean':'tas','tmin':'tasmin','tmax':'tasmax','ppt':'pr'}
    if variable == 'ppt':
        metric = 'total'
        units = 'mm'
    elif variable in ['tmean','tmax','tmin']:
        metric = 'mean'
        units = 'C'

    # write to disk
    try:
        if not os.path.exists( output_path ):
            os.makedirs( output_path )
    except:
        pass # deal with parallel processing in a clunky way

    month = os.path.basename(ak_fn).split('.')[0].split('_')[-1]
    if len(month) == 1:
        month = '0{}'.format(month)

    output_filename = os.path.join( output_path, '{}_{}_{}_akcan_prism_{}_1961_1990.tif'.format( varname_lookup[variable], metric, units, month ) )
    with rasterio.open( output_filename, 'w', **meta ) as out:
        out.write( merged )

    return output_filename


if __name__ == '__main__':
    import os, rasterio
    import pandas as pd
    import numpy as np

    base_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/raw/prism'
    output_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/prism'
    template_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_ar5_IPSL-CM5A-LR_rcp26_01_2006.tif'

    # use some hardwired functions to list and sort the data groups to a dict.
    files = filelister( base_dir )
    variables = ['ppt','tmean','tmin','tmax']

    for variable in variables:
        ak_files = files['{}-{}'.format('alaska',variable)]
        can_files = files['{}-{}'.format('canada',variable)]
        
        for ak, can in zip(ak_files, can_files):
            print( 'merging: {}'.format(variable) )
            
            # build output path
            varname_lookup = {'tmean':'tas','tmin':'tasmin','tmax':'tasmax','ppt':'pr'}
            output_path = os.path.join( output_dir, varname_lookup[variable] )
            
            # run the merge
            done = merge_ak_canada( ak, can, template_fn, variable, output_path )