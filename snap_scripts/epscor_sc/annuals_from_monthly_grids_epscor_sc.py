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

def get_monyear( fn ):
    '''
    specialized function to split the filename 
    pattern of SNAP data and extract 
    year and month information from it.
    '''
    fn, ext = os.path.splitext( fn )
    return fn.split( '_' )[-2:]

def get_month_seaon( fn ):
    ''' method to use a month number as a way to determine seasons '''
    # seasons
    seasonal_lookup = { 1:'DJF', 2:'DJF', 3:'MAM', 4:'MAM', 5:'MAM', \
                        6:'JJA', 7:'JJA', 8:'JJA',\
                         9:'SON', 10:'SON', 11:'SON', 12:'DJF' }

    fn = os.path.basename( fn )
    month, year = get_monyear( fn )
    return seasonal_lookup[ int(month) ]

def get_decade( fn ):
    month, year = get_monyear( fn )
    return year[:3]+'0s'

def make_grouper_df( files ):
    '''
    parse filenames to some of the common groupers

    ARGUMENTS:
    ----------
    files = [list] str filepaths to group

    RETURNS:
    --------
    pd.DataFrame with grouping variables as columns

    '''
    # some of the groupers we want to have
    months = [ get_monyear( fn )[0] for fn in files ]
    years = [ get_monyear( fn )[1] for fn in files ]
    seasons = [ get_month_seaon( fn ) for fn in files ]
    decades = [ get_decade( fn ) for fn in files ]
    columns = [ 'files', 'months', 'years', 'decades', 'seasons' ]
    dat = {'months':months, 'years':years, 'decades':decades, 
            'seasons':seasons, 'files':files}
    return pd.DataFrame( dat, columns=columns )

def _aggfun( group, out_fn, metric='mean', fn_col='files', decimals=1 ):
    '''
    aggregate the data in each group of filenames
    '''
    import rasterio
    metric_switch = { 'mean':np.mean, 'total':np.sum, 'min':np.min, 'max':np.max }

    name, df = group
    files = df[ fn_col ].tolist()
    
    with rasterio.open( files[0] ) as rst_tmp:
        meta = rst_tmp.meta
        mask = rst_tmp.read_masks( 1 ) == 0 # this can get you into trouble

    meta.update( compress='lzw' )
    if 'transform' in meta.keys():
        meta.pop( 'transform' )

    agg = metric_switch[ metric ]([ rasterio.open( fn ).read( 1 ) for fn in files ], axis=0 )

    agg = np.round( agg, decimals=decimals ) # round it
    if decimals < 1:
        agg = np.rint( agg )

    agg[ mask ] = meta[ 'nodata' ]

    # write it to disk
    dirname = os.path.dirname( out_fn )
    try:
        if not os.path.exists( dirname ):
            os.makedirs( dirname )
    except:
        pass

    with rasterio.open( out_fn, 'w', **meta ) as out:
        out.write( agg, 1 )
    return out_fn

def _aggfun_wrap( args ):
    return _aggfun( **args )


# def aggregate( groups, metric='mean', fn_col='files' ):
#   '''
#   groups = [list] of tuples where ( id, DataFrame )
#   '''
#   from pathos.mp_map import mp_map

#   done = mp_map( _aggfun, groups, nproc=32, **{'metric':metric, 'fn_col':fn_col} )


def generate_filenames( base_dir, output_dir, groups, freq ):
    names = [(name, group.iloc[0,:]['files']) for name, group in groups ]
    out_filenames = []
    for name, fn in names:
        fn, ext = os.path.splitext( fn )
        dirname, basename = os.path.split( fn )
        split = basename.split('_')
        year, month = [split.pop( len( split )-1 ) for i in range( 2 )]
        if 'seasons' in freq:
            newbase = '_'.join( split + [ name, year] ) + ext
        elif 'years':
            newbase = '_'.join( split + [ name ] ) + ext
        else:
            newbase = '_'.join( split+[ month, name ] ) + ext

        out_filenames = out_filenames + [os.path.join( dirname.replace( base_dir, output_dir ), newbase )]
    return out_filenames

if __name__ == '__main__':
    import os, rasterio, itertools
    import pandas as pd
    import numpy as np
    from pathos.mp_map import mp_map

    base_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled_FINAL_OCT'
    output_dir = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids_FINAL_OCT/annuals'
    # models = [ 'GFDL-CM3', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'GISS-E2-R', 'NCAR-CCSM4', '5ModelAvg' ]
    models = [ 'ts323' ]
    variables = [ 'tas','tasmin','tasmax', 'pr' ]
    # scenarios = [ 'historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85' ]
    scenarios = [ 'historical' ]

    for model, scenario, variable in itertools.product( models, scenarios, variables ):
        base_path = os.path.join( base_dir, model, scenario, variable )
        # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/downscaled/GFDL-CM3/rcp60/tas'
        if scenario == 'historical':
            begin = 1900
            end = 2005

            if model == 'ts323':
                begin = 1901
                end = 2014
        else:
            begin = 2006
            end = 2100

        if variable == 'pr':
            metric = 'total'
        else:
            metric = 'mean'

        files = sort_files( only_years( list_files( base_path ), begin, end ) )
        df = make_grouper_df( files )
        
        # the idea here is to use PANDAS groupby functions with the colums to produce groups of filenames
        # grouper for decadal seasonals
        # grouper = df.apply( lambda x: x['decades'] + x['seasons'], axis=1 )
        
        # # annuals
        grouper = df[ 'years' ]
        groups = [ group for group in df.groupby( grouper ) ] # unpack the groupby to a list
        
        # make some output filenames
        filenames = generate_filenames( base_dir, output_dir, groups, 'years' )

        if variable == 'pr':
            decimals = 0
        else:
            decimals = 1

        args_list = [ {'group':i,'out_fn':j, 'metric':metric, 'decimals':decimals } for i,j in zip( groups, filenames ) ]
        done = mp_map( _aggfun_wrap, args_list, nproc=32 )



