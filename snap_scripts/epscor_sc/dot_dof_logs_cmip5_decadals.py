# # # # 
# CALCULATE DOF/DOT/LOGS FROM TAS TIMESERIES
# # # # 
def tfg_days( x, err='off' ):
    ''' calculate DOF/DOT/LOGS for a vector of 12 chronological monthly values '''
    import itertools
    import numpy as np

    # filter the div by zero and comparison with np.nan warnings from numpy
    if err == 'off':
        np.warnings.filterwarnings( "ignore", category=RuntimeWarning )

    x[ x == 0 ] = -0.0001 # need to treat zero as freezing (working with signs)

    # positive or negative monthly temps
    s1 = np.sign( x )
    # products of consecutive months' signs: positive indicates no change; negative indicates a potential freeze or thaw transition
    s = s1[:11] * s1[1:]
    idx, = np.where( s < 0 )

    # may be length zero (no transitions)
    ind = np.sort( np.concatenate( [idx, idx+1] ) )

    if np.any( np.isnan( x ) == True ): # ignore cells with missing data
        dot, dof, grow = itertools.repeat( np.array([np.nan]), 3 )
        case = 1

    elif (len(ind) == 0) & (s1[0] > 0): # no transitions: all positive temps means no freeze day
        dot = np.array([0])
        dof, grow = itertools.repeat( np.array([365]), 2 )
        case = 2

    elif (len(ind) == 0) & (s1[0] < 0): # no transitions: all negative temps means no thaw day
        dot = np.array([365])
        dof, grow = itertools.repeat( np.array([0]), 2 )
        case = 3

    # [ML FIXED]
    elif len(ind) == 2: # only one transition during the year, thawing or freezing
        # places where we know the ground freezes and thaws, 
        #  but during a specific 12 months we just don't happen to witness both
        # only thaw occurs
        if x[ ind[0] ] < 0:
            # [ml] note:((ind[0]+1)-1) is ind[0]+1 is the month number and minus 1 is to get to previous month
            #      we could make that a call to a months array -- months = range(1, 12+1)
            dot = 15 + 30 * ((ind[0]+1)-1) - np.round( x[ ind[0] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
            dof = np.array([350]) # 350: we know the ground freezes so we use 350 rather than the special 365
            grow = dof - dot
            case = 4

        # only freeze occurs
        if x[ ind[0] ] > 0: 
            dof = 350 - 30 * (12-ind[1]-1) - np.round( x[ ind[1] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
            dot = np.array([15]) # 15: we know the ground thaws so we use 15 rather than the special 0
            grow = dof - dot
            case = 5

    # [ML FIXED]
    elif (len(ind) == 4 ) & (s1[0] < 0): # two transitions occur: thaw, then freeze (this is the ideal case; everything else is an idiosyncratic edge case)
        # [ml] note:((ind[0]+1)-1) is ind[0]+1 is the month number and minus 1 is to get to previous month
        #      we could make that a call to a months array -- months = range(1, 12+1)
        dot = 15 + 30 * ((ind[0]+1)-1) - np.round( x[ ind[0] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
        dof = 350 - 30 * (12-ind[3]-1) - np.round( x[ ind[3] ] / (np.diff( x[ ind[2:4] ] ) / 30.0), decimals=0 )
        grow = dof - dot
        case = 0
    
    # [ML FIXED]
    elif (len(ind) == 4) & (s1[0] > 0): # two transitions occur but backward to what is expected; freeze, then thaw
        if( ind[0] >= 7 ): # freeze occurs in second half of year as expected; late thaw is spurious
            # dof = 350 - 30 * (12-ind[1]-1) - np.round( x[ ind[1] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
            dof = 350 - 30 * (12-ind[1]-1) - np.round( x[ ind[1] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
            dot = np.array([15]) # ignore spurious post-freeze thaw; treat as early, unobserved thaw
            grow = dof - dot
            case = 6
        
        if ind[0] <= 6: # spurious freeze occurs in first half of year; thaw probably fine
            dot = 15 + 30 * ((ind[2]+1)-1) - np.round( x[ ind[2] ] / (np.diff( x[ ind[2:4] ]) / 30.0), decimals=0 )
            dof = np.array([350]) # ignore spurious early freeze; treat as late, unobserved freeze
            grow = dof - dot
            case = 7

    # [ML FIXED]    
    elif len(ind) > 4: # more than two transitions; at least one definitely spurious
        # [MATT Q]:
        #       what is the prepending 0 below? and what is its intention?
        #       what do u do if there is a use-case where idx-0 is already chosen? Py is ZERO-anchored...        
        ind2, = np.where( s < 0 )
        ind2 = ind2 + 1
        ind2 = np.insert( ind2, 0, np.array([0]) )

        # [ml] m1, m2 are month indexes
        m1, = np.where( np.diff( ind2 ) == np.max( np.diff( ind2 ) ) )
        m1 = m1[-1] + 1
        m2, = np.where( np.delete(np.diff( ind2 ), (m1-1)-1) == max( np.delete(np.diff( ind2 ), (m1-1)-1)) )
        m2, = np.where( np.delete(np.diff( ind2 ), (m1-1)) == max( np.delete(np.diff( ind2 ), (m1-1))) )
        m2 = m2[-1] + 1

        if m1 == m2:
            m2 = m2 - 1

        ind2 = ind2[ np.sort( np.append( m1, m2 ) ) ]
        ind = np.sort( np.append(ind2, ind2+1) ) - 1
       
        dot = 15 + 30 * (ind[1]-1) - np.round( x[ind[1]-1] / (np.diff( x[ ind[:2] ] ) / 30.0), 0) # [ml] SOME WEIRD -1's here...
        dof = 350 - 30 * (12-ind[3]-1) - np.round( x[ind[3]] / (np.diff( x[ ind[2:4] ] ) / 30.0), 0)
        grow = dof - dot
        case = 8

    else:
        dot, dof, grow = itertools.repeat( np.array([np.nan]), 3 )
    return np.concatenate([dof, dot, grow])

def read_arr( fn, band=1 ):
    ''' 
    read the array from a GTiff using rasterio without mem leaks 
    and return a tuple of (arr, meta)
    where `meta` is the file metadata 
    '''
    with rasterio.open( fn ) as rst:
        arr = rst.read( band )
        meta = rst.meta
    return (arr, meta)

def rasterize_shp( shp_fn, arr, affine, fill=0 ):
    '''
    convert a shapefile into a raster using a template arr and
    affine transform.
    '''
    from rasterio import features
    import geopandas as gpd
    shp = gpd.read_file( shp_fn )
    geoms = [ (geom, idx+1) for idx, geom in enumerate( shp.geometry ) ]

    mask = features.rasterize( geoms, out_shape=arr.shape,
                            fill=fill, transform=affine,
                            dtype='float32', all_touched=True )
    return mask
def get_monyear( fn ):
    '''
    specialized function to split the filename 
    pattern of SNAP data and extract 
    year and month information from it.
    '''
    fn, ext = os.path.splitext( fn )
    return fn.split( '_' )[-2:]

def get_decade( fn ):
    month, year = get_monyear( fn )
    return str(year[:-3])+ '0s'
def get_scenario( fn ):
    ''' 
    specialized function to return the rcp or historical for the 
    data files scenario from the filename components.  Its Quasi-smart
    but not perfect.
    '''
    basename = os.path.basename( fn )
    basename, ext = os.path.splitext( basename )
    bsplit = basename.split('_')
    scenario, = [ i for i in bsplit if i.startswith('rcp') or i == 'historical' ]
    return scenario

# def run( fn, mask, band=1 ):
#     ''' run the extractor '''
#     arr, meta = read_arr( fn, band=band )
#     month, year = get_monyear( fn )
#     return pd.DataFrame({ '-'.join([month, year]):{int(i):arr[ mask == i ][0] 
#                                 for i in np.unique( mask ) if i > 0} }).T

def run_decade( filenames ):
    import os, rasterio
    import numpy as np

    # load the data to an ndarray
    arr = np.array([ read_arr( fn )[0] for fn in filenames ])
    mask = rasterio.open( fn ).read_masks( 1 )
    
    # run the function
    out_arr = np.apply_along_axis( tfg_days, axis=0, arr=arr )

    # filename work
    fn = filenames[0]
    dirname, basename = os.path.split( fn )
    basename, ext = os.path.splitext( basename )
    new_base = basename.replace( 'tas_mean_monthly_mean_C', 'dot_dof_grow_mean_decadal' ).replace( '_01_', '_' )
    new_dir = dirname.replace( 'tas', 'dot_dof_grow' )

    # build output directory if needed.
    if not os.path.exists( new_dir ):
        os.makedirs( new_dir )

    # copy metadata
    with rasterio.open( fn ) as rst:
        meta = rst.meta
        meta.update( compress='lzw', count=3 )
        _ = meta.pop( 'transform' )

    # write out file
    output_filename = os.path.join( new_dir, new_base + ext )
    with rasterio.open( output_filename, 'w', **meta ) as out:
        # make it float32?
        out_arr = out_arr.astype( np.float32 )
        ind0, ind1 = np.where( mask == 0 )
        out_arr[ ..., ind0, ind1 ] = meta[ 'nodata' ]
        out.write( out_arr )

    return output_filename

if __name__ == '__main__':
    import os, rasterio, glob
    import numpy as np
    import pandas as pd
    from pathos.mp_map import mp_map
    from functools import partial
    from downscale.utils import sort_files, only_years
    import matplotlib
    matplotlib.use( 'agg' )
    from matplotlib import pyplot as plt
    import argparse

    # parse some args
    parser = argparse.ArgumentParser( description='downscale the AR5-CMIP5 data to the AKCAN extent required by SNAP' )
    parser.add_argument( "-b", "--base_path", action='store', dest='base_path', type=str, help="path to the directory where the decadal monthly downscaled data are stored" )
    parser.add_argument( "-m", "--model", action='store', dest='model', type=str, help="model name (exact)" )
    
    # parse the args and unpack
    args = parser.parse_args()
    base_path = args.base_path
    model = args.model

    # # TESTING STUFF 
    # base_path = '/workspace/Shared/Tech_Projects/EPSCoR_Southcentral/project_data/derived_grids/decadal_monthlies'
    # model = 'GFDL-CM3'
    # # END TESTING

    # # # # TESTING MATTS STUFF
    # base_path = '/Data/Base_Data/Climate/AK_CAN_2km/projected/AR5_CMIP5_models/rcp60/IPSL-CM5A-LR/derived/tas/decadal_mean'
    # model = 'IPSL-CM5A-LR'

    # RUN ALL SCENARIOS IN A SINGLE NODE FOR SINGLE MODEL
    # # list and sort the files from the directory
    files = [ os.path.join( r, f ) for r,s,files in os.walk( base_path ) for f in files if f.endswith('.tif') and 'tas_' in f and model in f ]
    # # TESTING ONE BELOW
    # files = [ os.path.join( r, f ) for r,s,files in os.walk( base_path ) for f in files if f.endswith('.tif') and 'tas_' in f and model in f and 'rcp60' in f and '_201' in f]
    # # END TESTING
    scenarios = [ get_scenario(fn) for fn in files ]
    file_groups = [ [ sorted( y.tolist() ) for x,y in j.groupby( [ fn.split('.')[0].split('_')[-1] for fn in j ] )] for i,j in pd.Series( files ).groupby( scenarios ) ]
    file_groups = [ j for i in file_groups for j in i ]
    
    # groupby decade
    # decades = [ fn.split('.')[0].split('_')[-1] for fn in files ]
    # file_groups = [ j.tolist() for i,j in pd.Series( files ).groupby( decades ) ]

    # NEW RUNNNER!
    done = mp_map( run_decade, file_groups, nproc=32 )

