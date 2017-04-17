# compare precip OLD / NEW

def sort_files( files, split_on='_', elem_month=-2, elem_year=-1 ):
    '''
    sort a list of files properly using the month and year parsed
    from the filename.  This is useful with SNAP data since the standard
    is to name files like '<prefix>_MM_YYYY.tif'.  If sorted using base
    Pythons sort/sorted functions, things will be sorted by the first char
    of the month, which makes thing go 1, 11, ... which sucks for timeseries
    this sorts it properly following SNAP standards as the default settings.

    ARGUMENTS:
    ----------
    files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
    split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
    elem_month = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
        default:-2. For SNAP standard.
    elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
        default:-1. For SNAP standard.

    RETURNS:
    --------
    sorted `list` by month and year ascending. 

    '''
    import pandas as pd
    months = [ int(fn.split('.')[0].split( split_on )[elem_month]) for fn in files ]
    years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
    df = pd.DataFrame( {'fn':files, 'month':months, 'year':years} )
    df_sorted = df.sort_values( ['year', 'month' ] )
    return df_sorted.fn.tolist()

def only_years( files, begin=1901, end=2100, split_on='_', elem_year=-1 ):
    '''
    return new list of filenames where they are truncated to begin:end

    ARGUMENTS:
    ----------
    files = [list] list of `str` pathnames to be sorted by month and year. usually from glob.glob.
    begin = [int] four digit integer year of the begin time default:1901
    end = [int] four digit integer year of the end time default:2100
    split_on = [str] `str` character to split the filename on.  default:'_', SNAP standard.
    elem_year = [int] slice element from resultant split filename list.  Follows Python slicing syntax.
        default:-1. For SNAP standard.

    RETURNS:
    --------
    sliced `list` to begin and end year.
    '''
    import pandas as pd
    years = [ int(fn.split('.')[0].split( split_on )[elem_year]) for fn in files ]
    df = pd.DataFrame( { 'fn':files, 'year':years } )
    df_slice = df[ (df.year >= begin ) & (df.year <= end ) ]
    return df_slice.fn.tolist()

def list_files( base_path, ext='.tif' ):
    '''
    list files by walking a directory...

    ARGUMENTS:
    ----------
    base_path = [str] base directory where the data are stored
    ext = [str] extension to look for as a qualifier of data we want
                default: '.tif'
    RETURNS:
    --------
    list of filename str's in no particular order
    '''
    import os
    fn_list = []
    for root, subs, files in os.walk( base_path ):
        fn_list = fn_list + [ os.path.join( root, fn ) for fn in files if fn.endswith( ext ) ]
    return fn_list


if __name__ == '__main__':
    import rasterio, os, glob
    import numpy as np
    import pandas as pd
    from pathos.mp_map import mp_map

    variable = 'pr'
    # paths
    # mine = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_FINAL_OCT/GFDL-CM3/rcp60', variable )
    # matt = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/OLD_OCT_DELIVERY_KEEP_TMP/downscaled/GFDL-CM3/rcp60', variable )
    mine = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_FINAL_OCT/ts323/historical', variable )
    # matt = os.path.join( '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models/rcp60/GFDL-CM3', variable )
    matt = os.path.join( '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/ts323/historical', variable )

    # get the files
    globber = '*.tif'
    files = { group:sort_files( glob.glob( os.path.join( path, globber ) ) ) for group, path in zip(['mine', 'matt'],[mine, matt]) }

    # df it
    df = pd.DataFrame( files )
    
    # sliced = df.iloc[ slice( 0, 20 ) ]
    df = df.iloc[ slice( 0, 12 ) ]

    # open and stack
    def read_arr( fn ):
        with rasterio.open( fn ) as rst:
            arr = rst.read( 1 )
        return arr

    arr_list = [ mp_map( read_arr, df[i].tolist(), nproc=32 ) for i in ['mine', 'matt'] ]
    # arr_list = [ np.array([ rasterio.open(j).read(1) for j in df[i].tolist()]) for i in ['mine', 'matt'] ]

    # diff
    diff = np.array( arr_list[1] ) - np.array( arr_list[0] )

    def unique_mm( arr ):
        uniques = np.unique( arr )
        return uniques
        # return uniques.min(), uniques.max()

    uniques = [ unique_mm( i ) for i in diff ]


# # # # STACK em
    tmp_fn = df.iloc[0].tolist()[1]
    meta = rasterio.open( tmp_fn ).meta
    a = np.array(arr_list[0])
    b = np.array(arr_list[1])
    meta.update(compress='lzw', count=a.shape[0])
    out_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_FINAL_OCT_TESTING'
    for out_arr, out_fn in zip([a,b],[os.path.join( out_path, i+'.tif' ) for i in ['mine_3', 'matt_3']]):
        with rasterio.open( out_fn, 'w', **meta ) as out:
            out.write( out_arr )
    

dec_old = rasterio.open('/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models/rcp60/GFDL-CM3/derived/pr/decadal_mean/pr_decadal_mean_annual_total_mm_GFDL-CM3_rcp60_2010_2019.tif').read(1)
