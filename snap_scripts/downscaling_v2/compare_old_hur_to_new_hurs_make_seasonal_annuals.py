# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# DATA COMPARISONS OF NEW hurs AND OLD hur DOWNSCALED DATA
#  Michael Lindgren (malindgren@alaska.edu) September 2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def open_raster( fn, band=1 ):
    with rasterio.open( fn ) as rst:
        arr = rst.read( band )
    return arr

def seasonize( row ):
    if row['month'] in (12,1,2):
        season = 'DJF'
    if row['month'] in (3,4,5):
        season = 'MAM'
    if row['month'] in (6,7,8):
        season = 'JJA'
    if row['month'] in (9,10,11):
        season = 'SON'
    return season

def make_seasonal_average( files, model, scenario, season, output_path ):
    arr = np.round( np.mean([ open_raster( fn, band=1 ) for fn in files ], axis=0), 2 )
    basename = os.path.basename( files[0] )
    month_str = basename.split('_')[-2]
    basename = basename.replace('_'+month_str+'_', '_'+season+'_' )
    output_filename = os.path.join( output_path, model, scenario, basename )
    try:
        dirname = os.path.dirname( output_filename )
        if not os.path.exists( dirname ):
            _ = os.makedirs( dirname )
    except:
        pass

    with rasterio.open( files[0] ) as rst:
        meta = rst.meta.copy()
        mask = rst.read(1) == -9999
        meta.update( compress='lzw' )

    arr[ mask ] = -9999
    with rasterio.open( output_filename, 'w', **meta ) as out:
        out.write( arr.astype(np.float32), 1 )

    return output_filename

def season_group_ids( df ):
    # unpack
    groups = [ j for i,j in list(df.groupby(np.arange(0,df.shape[0])//3))]
    for count, df in enumerate(groups):
        df['season_id'] = count
    return pd.concat(groups)

def drop_firstyear( df ):
    if df.scenario.iloc[0] == 'historical':
        df.index = pd.date_range('1900-01-01', '2005-12-31', freq='M')
        df = df.loc['1900-03-01':'2005-11-30']
    else:
        df.index = pd.date_range('2006-01-01', '2100-12-31', freq='M')
        df = df.loc['2006-03-01':'2100-11-30']
    return df

def prep_filename_groups(group):
    return season_group_ids( drop_firstyear( group ) )

def make_args( row, output_path ):
    files = row['new_fn'].tolist()
    season = str(row['season'].iloc[0])
    model = row['model'].iloc[0]
    scenario = row['scenario'].iloc[0]
    return {'files':files, 'model':model, 'scenario':scenario, 'season':season, 'output_path':output_path}

def run( x ):
    return make_seasonal_average( **x )

if __name__ == '__main__':
    import os, glob
    import rasterio
    from rasterio.features import rasterize
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    import multiprocessing as mp
    from functools import partial
    import itertools

    # setup args
    ncpus = 32
    new_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary'
    old_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_old_hur_vap_Aug2018'
    output_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/hurs_seasonal_grids_testing'

    models = ['GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4']
    scenarios = ['historical','rcp26','rcp45','rcp60','rcp85']
    # seasons = {'DJF':(12,1,2),'MAM':(3,4,5),'JJA':(6,7,8),'SON':(9,10,11)}
    seasons = ['DJF','MAM','JJA','SON']

    # list the data
    new_files = [ os.path.join(r,fn) for r,s,files in os.walk(new_path) for fn in files if fn.endswith('.tif') ]
    old_files = [ os.path.join(os.path.dirname(fn).replace(new_path, old_path),'hur',os.path.basename(fn).replace('hurs','hur')).replace('/hurs/','/hur/') for fn in new_files ]

    # build into a dataframe
    colnames = ['variable','metric','units','group','model','scenario','month','year']
    files_df = pd.DataFrame([ os.path.basename(fn).split('.')[0].split('_') for fn in new_files ], columns=colnames)
    files_df['new_fn'] = new_files
    files_df['old_fn'] = old_files
    files_df['month'] = files_df['month'].astype(int)
    files_df['year'] = files_df['year'].astype(int)

    # make seasons
    files_df['season'] = files_df.apply( seasonize, axis=1 )

    # regroup and add a unique value for each season to get them chronologically
    files_df = files_df.sort_values(['model','scenario','year','month'])
    grouped = files_df.groupby(['model','scenario'])
    df_list = [ group for rowid, group in grouped ]

    pool = mp.Pool( 32 )
    out = pool.map( prep_filename_groups, df_list )
    pool.close()
    pool.join()

    for count, df in enumerate(out):
        print(count)
        groups = df.groupby('season_id')
        args = [ make_args(df,output_path) for rowid,df in groups ]

        pool = mp.Pool( 32 )
        final_out = pool.map( run, args )
        pool.close()
        pool.join()
