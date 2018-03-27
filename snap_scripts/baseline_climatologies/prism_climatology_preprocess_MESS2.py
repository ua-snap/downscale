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

# def reproject_canada_to3338( fn, outbounds, scalevals=False, dtype='float32' ):
#     '''
#     reproject and scale the values to floats from ints 
#     '''
#     from rasterio.warp import calculate_default_transform, reproject, Resampling

#     src_crs = 'EPSG:4326'
#     dst_crs = 'EPSG:3338'

#     with rasterio.open( fn ) as src:
#         transform, width, height = calculate_default_transform(
#             src.crs, dst_crs, src.width, src.height, *src.bounds)
#         kwargs = src.meta.copy()
#         kwargs.update({
#             'crs': dst_crs,
#             'transform': transform,
#             'width': width,
#             'height': height,
#             'driver':'MEM',
#             'dtype':dtype
#         })
#         # open an in-memory dataset in the destination crs
#         dst = rasterio.open('', 'w', **kwargs)
#         # reproject the data to dst
#         reproject(
#             source=rasterio.band(src, 1),
#             destination=rasterio.band(dst, 1),
#             src_transform=src.transform,
#             src_crs=src_crs,
#             dst_transform=transform,
#             dst_crs=dst_crs,
#             resampling=Resampling.nearest)
        
#         if scalevals: # scale the temperature PRISM values for Canada...
#             arr = dst.read(1)
#             arr[ arr != dst.nodata ] = arr[ arr != dst.nodata ] / 10.0
#             dst.write( arr, 1 )

#         return dst

def reproject_canada_to3338( fn, template_fn, scalevals=False, dtype='float32' ):
    '''
    reproject and scale the values to floats from ints 
    '''
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    src_crs = {'init':'EPSG:4326'}
    dst_crs = {'init':'EPSG:3338'}

    with rasterio.open( fn ) as src:

        with rasterio.open( template_fn ) as tmp:
            tmp_meta = tmp.meta.copy()
        # open an in-memory dataset in the destination crs
        # tmp_meta.update( nodata=-9999, dtype=dtype, crs=dst_crs, driver='MEM' )
        tmp_meta.update( nodata=-9999, dtype=dtype, crs=dst_crs, driver='GTiff' )
        # dst = rasterio.open( '', 'w', **tmp_meta )
        dst = rasterio.open( '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/TEST_PRISM_CANADA.tif', 'w', **tmp_meta )
        # reproject the data to dst
        reproject(
            source=rasterio.band(src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=src.transform,
            src_crs=src_crs,
            dst_transform=dst.meta['transform'],
            dst_crs=dst_crs,
            resampling=Resampling.nearest)
        
        if scalevals: # scale the temperature PRISM values for Canada...
            arr = dst.read()
            arr[ arr != dst.nodata ] = arr[ arr != dst.nodata ] / 10.0
            dst.write( arr )

        return dst

# def open_ak_precip( fn ):
#     ''' 
#     deal with diffs between dtypes between reproj canada and raw ak
#     since float32 is our base dtype for ALF / IEM I will stick with 
#     this as the dtype to store these data, despite being integers in 
#     reality.
#     '''
#     with rasterio.open( fn ) as rst:
#         meta = rst.meta
#         arr = rst.read(1).astype( np.float32 )

#         meta.update({
#             'driver':'MEM',
#             'dtype':'float32'
#             })

#     with rasterio.open( '', 'w', **meta ) as out:
#         out.write( arr, 1 )
#     return out

def bounds_to_extent( bounds ):
    '''
    take input rasterio bounds object and return an extent
    '''
    l,b,r,t = bounds
    return [ (l,b), (r,b), (r,t), (l,t), (l,b) ]

def merge_ak_canada( ak_fn, can_fn, template_fn, variable, output_path ):
    from rasterio.merge import merge
    from shapely.geometry import Polygon
    from shapely.geometry.geo import mapping
    from rasterio.mask import mask
    
    # reproject canada to 3338 and rescale temperature values from ints to floats
    if variable == 'ppt':
        scalevals = False
        dtype = 'int32'
    else: # temperature variables
        scalevals = True
        dtype = 'float32'

    # get the template raster bounds...
    with rasterio.open( template_fn ) as tmp:
        outbounds = tmp.bounds
        geom = Polygon( bounds_to_extent( outbounds ) )

    can = reproject_canada_to3338( can_fn, template_fn, scalevals=scalevals, dtype=dtype )
    
    # open alaska AAI grid
    ak = rasterio.open( ak_fn ) # already in 3338

    # merge 'em
    merged, transform = merge( [ak,can], res=(2000,2000), nodata=-9999 )
    # update metadata
    meta = ak.meta
    count, height, width = merged.shape
    meta.update({
        'driver':'GTiff',
        'height':height,
        'width':width,
        'crs':{'init':'EPSG:3338'},
        'transform':transform,
        'compress':'lzw'
        })

    # output filenaming
    varname_lookup = {'tmean':'tas','tmin':'tasmin','tmax':'tasmax','ppt':'pr'}
    if variable == 'ppt':
        metric = 'total'
        units = 'mm'
    elif variable in ['tmean','tmax','tmin']:
        metric = 'mean'
        units = 'C'

    # crop to the extent that we want for AKCAN 2km SNAP data
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    tmp_meta = meta.copy()
    tmp_meta.update( driver='MEM' )
    with rasterio.open( '', 'w', **tmp_meta ) as src:
        out_image, out_transform = mask(src, [mapping(geom)], crop=True)
        out_meta = src.meta.copy()

    meta.update({'driver': 'GTiff',
                 'height': out_image.shape[1],
                 'width': out_image.shape[2],
                 'transform': out_transform})

    # with rasterio.open( "/tmp/masked.tif", "w", **out_meta ) as dest:
    #     dest.write(out_image)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

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
    output_dir = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/climatologies/prism_v4'

    # template raster used for output extent/crs/origin -- from version1 downscaling...
    # template_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_AR5_CCSM4_rcp26_01_2006.tif'
    template_fn = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/templates/akcan_2km/tas_mean_C_ar5_IPSL-CM5A-LR_rcp26_01_2006.tif'
    
    # use some hardwired functions to list and sort the data groups to a dict.
    files = filelister( base_dir )
    # variables = ['ppt','tmean','tmin','tmax']
    variables = ['tmin','tmax']
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
            break
        break