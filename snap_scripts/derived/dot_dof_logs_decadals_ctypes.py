# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PYTHON VERSION OF MATTHEW LEONAWICZ's DOF/DOT/LOGS R Script
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Author: Michael Lindgren -- 2019 -- malindgren@alaska.edu
# LICENSE: MIT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

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
        dot = 0
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
        # print( "Condition unconsidered: {}".format( x.strip() ) )
    return np.array([dof, dot, grow])

def open_raster( fn ):
    ''' cleanly open the raster '''
    with rasterio.open( fn ) as rst:
        arr = rst.read(1).copy()
    return arr

def run_tfg_days( x ):
    '''wrapper to run the tfg days using the apply_along_axis mechanics'''
    out = np.full( 3, np.nan )
    if not np.isnan(x).any():
        vals = tfg_days( x ).squeeze()
        out[[0,1,2]] = vals
    return out

def run( x ):
    ''' run it '''
    # unpack
    x, out_base_fn = x
    year, sub_df = x
    
    nodata = -9999 # hardwired for SNAP data

    files = sub_df.fn.tolist()
    arr = np.array([ open_raster(fn) for fn in files ])
    arr[ arr == nodata ] = np.nan

    out_arr = np.apply_along_axis( run_tfg_days, axis=0, arr=arr )
    out_arr = out_arr[:3,...] # slice it back

    for idx,metric in enumerate(['dof','dot','logs']):
        out_fn = out_base_fn.format(metric, metric, str(year))
        dirname = os.path.dirname( out_fn )
        
        try:
            if not os.path.exists( dirname ):
                _ = os.makedirs( dirname )
        except:
            pass

        # get meta and mask for output arr and GTiff generation
        with rasterio.open(files[0]) as rst:
            meta = rst.meta.copy()
            meta.update( compress='lzw', dtype=np.int32 )
            mask = rst.read(1) == rst.nodata
        
        # write to GTiff
        with rasterio.open( out_fn, 'w', **meta ) as out:
            cur_arr = out_arr[ idx, ... ]
            cur_arr[ mask ] = -9999
            out.write( cur_arr.astype(np.int32) , 1 )

    return out_fn

def run_par_update_arr( idx ):
    tmp_out = np.ctypeslib.as_array(out_shared)
    i,j = idx
    tmp_out[:,i,j] = run_tfg_days( arr[:,i,j] )
    

if __name__ == '__main__':

    import os, glob, rasterio, itertools
    import numpy as np
    import pandas as pd
    import xarray as xr
    import multiprocessing as mp
    from multiprocessing import sharedctypes

    ncpus = 64
    models = ['GFDL-CM3','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NCAR-CCSM4','5ModelAvg',]
    # models = ['CRU-TS40',]
    scenarios = ['historical','rcp26','rcp45','rcp60','rcp85',]
    # scenarios = ['historical']
    
    for model, scenario in itertools.product( models, scenarios ):
        print(model, scenario)
        base_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/downscaled/{}/{}/tas'.format(model,scenario)
        out_path = '/workspace/Shared/Tech_Projects/DeltaDownscaling/project_data/derived_v2_parallel/{}/{}'.format(model,scenario)
        
        files = glob.glob(os.path.join(base_path, '*.tif'))
        colnames = ['variable', 'metric', 'units', 'project', 'model', 'scenario', 'month', 'year']
        df = pd.DataFrame([ os.path.basename(fn).split('.')[0].split('_')for fn in files ], columns=colnames)
        df['fn'] = files
        df.month = df.month.astype(int)
        df.year = df.year.astype(int)
        df = df.sort_values(['year', 'month']).reset_index( drop=True )
        df['decade'] = df.year.apply(lambda x: str(x)[:3]+'0s')
        args = list(df.groupby('decade'))
        args = [ ((year,sub_df),os.path.join(out_path,'{}','{}_'+model+'_'+scenario+'_{}.tif' )) for year,sub_df in args ]
        
        # get template file information
        with rasterio.open(files[0]) as tmp:
            meta = tmp.meta.copy()
            shape = tmp.shape
            rows,cols = shape
            arr = tmp.read(1)
            arr[arr == tmp.nodata] = np.nan

        for arg in args:
            x, out_base_fn = arg
            year, sub_df = x

            nodata = -9999 # hardwired for SNAP data

            files = sub_df.fn.tolist()
            times = pd.DatetimeIndex([ pd.Timestamp.strptime('{1}-{0}-15'.format(*os.path.basename(fn).split('.')[0].split('_')[-2:]), '%Y-%m-%d') for fn in files ])
            arr = np.array([ open_raster(fn) for fn in files ])
            arr[ arr == nodata ] = np.nan

            # make xarray dset
            ds = xr.Dataset({'dat':(['time','yc', 'xc'], arr.copy())},
                            coords={'xc': ('xc', np.arange(cols)),
                                    'yc': ('yc', np.arange(rows)),
                                    'time':times })

            # groupby and average to decadal series (12 values)
            da = ds.dat.groupby('time.month').mean(axis=0)
            arr = da.values.copy()
            del da, ds # cleanup

            indexes = list(zip(*np.where(~np.isnan(arr[0]))))
            # indexes = list(np.ndindex(arr.shape[-2:])) # old way

            # make the output arrays?
            count = 3
            out = np.ctypeslib.as_ctypes(np.zeros((count,rows,cols), dtype=np.float))
            out_shared = sharedctypes.RawArray(out._type_, out)
            p = mp.Pool(ncpus)
            p.map( run_par_update_arr, indexes )
            p.close()
            p.join()

            # bring these c-types arrays back to numpy arrays.
            out = np.ctypeslib.as_array(out_shared).astype(np.float32)

            # update the np.nans to something more useable and SNAP-ish
            for i in out:
                i[np.isnan(arr[0])] = -9999

            # dump to disk
            for idx,metric in enumerate(['dof','dot','logs']):
                out_fn = out_base_fn.format(metric, metric, str(year))
                dirname = os.path.dirname( out_fn )
                
                try:
                    if not os.path.exists( dirname ):
                        _ = os.makedirs( dirname )
                except:
                    pass

                # write to GTiff
                meta.update( compress='lzw', dtype=np.int32, count=1, nodata=-9999 )
                with rasterio.open( out_fn, 'w', **meta ) as out_rst:
                    cur_arr = out[ idx, ... ]
                    out_rst.write( cur_arr.astype(np.int32), 1 )

