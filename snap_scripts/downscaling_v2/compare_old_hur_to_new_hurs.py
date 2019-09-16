# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# DATA COMPARISONS OF NEW hurs AND OLD hur DOWNSCALED DATA
#  Michael Lindgren (malindgren@alaska.edu) September 2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def open_raster( fn, band=1 ):
    with rasterio.open( fn ) as rst:
        arr = rst.read( band )
    return arr

def get_mean_area( arr, mask ):
    return np.mean(arr[ mask != 0 ])

def run( fn, mask ):
    arr = open_raster( fn )
    return get_mean_area( arr, mask )

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

    # area shapefile?
    shp_fn = '/Data/Base_Data/GIS/GIS_Data/Vector/AlaskaWatersheds/yukon_river_watershed_boundary.shp'

    models = ['GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4']
    scenarios = ['historical','rcp26','rcp45','rcp60','rcp85']

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

    # rasterize shapefile to mask or take the one from the file
    if shp_fn is not None:
        shp = gpd.read_file( shp_fn )
        shapes = list( shp.geometry )

        # mask
        with rasterio.open( new_files[0] ) as tmp:
            meta = tmp.meta.copy()
            tmp_arr = np.empty_like( tmp.read(1) )

        mask = rasterize( shapes, out_shape=tmp_arr.shape, fill=0, out=None, transform=meta['transform'], 
                            all_touched=False, default_value=1, dtype=None )
    else:
        with rasterio.open( new_files[0] ) as tmp:
            meta = tmp.meta.copy()
            mask = tmp.read_masks( 1 )

    # make a function to prepare for parallel procesing by group.
    f = partial( run, mask=mask )
    
    out = dict() # hold outputs
    for index_names, cur_df in files_df.groupby(['model','scenario']):
        print( index_names )
        
        # sort it
        new_df = cur_df.sort_values(['year','month']).copy()
        
        # new
        pool = mp.Pool( ncpus )
        new = pool.map( f, new_df['new_fn'].tolist() )
        pool.close()
        pool.join()
        
        # old
        pool = mp.Pool( ncpus )
        old = pool.map( f, new_df['old_fn'].tolist() )
        pool.close()
        pool.join()
        
        out[ '_'.join(index_names) ] = {'new':new, 'old':old }

    # futures
    years = (2006,2100)
    df_list = list()
    for key in out:
        if key.split('_')[1] != 'historical':
            df = pd.DataFrame(out[key])
            df.columns = ['_'.join([key,i]) for i in df.columns]
            df_list = df_list + [df]
    
    out_df_futures = pd.concat( df_list, axis=1 )
    out_df_futures.index = pd.date_range(str(years[0]), str(years[1]+1), freq='M')
    output_filename = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary/TESTING_NEW_OLD_hurs_hur/futures_yukonwatershed_areamean_allmodels.csv'
    out_df_futures.to_csv( output_filename, sep=',' )

    # historical
    years = (1900,2005)
    df_list = list()
    for key in out:
        if key.split('_')[1] == 'historical':
            df = pd.DataFrame(out[key])
            df.columns = ['_'.join([key,i]) for i in df.columns]
            df_list = df_list + [df]

    out_df_historicals = pd.concat( df_list, axis=1 )
    out_df_historicals.index = pd.date_range(str(years[0]), str(years[1]+1), freq='M')
    output_filename = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary/TESTING_NEW_OLD_hurs_hur/historicals_yukonwatershed_areamean_allmodels.csv'
    out_df_historicals.to_csv( output_filename, sep=',' )


    # # now lets do some plotting:
    # '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary/TESTING_NEW_OLD_hurs_hur/futures_akcan_areamean_allmodels.csv'
    # '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary/TESTING_NEW_OLD_hurs_hur/historicals_yukonwatershed_areamean_allmodels.csv'
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt

    out_df_futures = pd.read_csv('/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled_hurs_temporary/TESTING_NEW_OLD_hurs_hur/futures_akcan_areamean_allmodels.csv', index_col=0, parse_dates=True)
    futures_df = out_df_futures.loc[slice('2010','2099')]
    months = [ i.month for i in futures_df.index ]
    years = [ str(i.year)[:-1]+'0s' for i in futures_df.index ]
    
    futures_df['month'] = months
    futures_df['year'] = years

    l = futures_df.groupby(['year','month']).mean()
    scenarios = ['rcp26','rcp45','rcp60','rcp85']
    for model, scenario in itertools.product(models, scenarios):
        sub_df = l[[i for i in l.columns if scenario in i and model in i]]
        sub_df.index = [ '-'.join([str(j),str(i)]) for i,j in sub_df.index.tolist() ]
        sub_df.columns = [ i.split('_')[-1] for i in sub_df.columns ]

        sub_df.plot(kind='line', title='Comparison Of New and Old Relative Humidity Versions\n{} - {}\nMonthly Decadal Averages 2010s-2090s'.format(model.upper(), scenario.upper()))
        plt.savefig('compare_relative_humidity_decadals_{}_{}.png'.format(model,scenario))
        plt.close()

    # END


