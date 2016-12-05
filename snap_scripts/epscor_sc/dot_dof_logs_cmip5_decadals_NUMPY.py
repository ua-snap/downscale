# MAKE NUMPY RUNNER make x 3D

def tfg_days_3d( x, err='off' ):
    ''' calculate DOF/DOT/LOGS for a vector of 12 chronological monthly values '''
    import itertools
    import numpy as np

    # filter the div by zero and comparison with np.nan warnings from numpy
    if err == 'off':
        np.warnings.filterwarnings( "ignore", category=RuntimeWarning )

    # need to treat zero as freezing (working with signs)
    if np.any(x) == 0:
        x[ x == 0 ] = -0.0001

    nlayers, rows, cols = x.shape

    # make some 2-D dof / dot arrays
    dof = np.copy( x[0,...] )
    dot = np.copy( x[0,...] )

    # positive or negative monthly temps
    s1 = np.sign( x )
    # products of consecutive months' signs: positive indicates no change; negative indicates a potential freeze or thaw transition
    s = s1[:11, ...] * s1[1:, ...]

    # [ HARDWIRED ] pull out a mask using the known SNAP oob value
    oob_value = -3.39999995e+38
    oob_mask = (x[0,...] != oob_value).astype( int ) # pull out a mask for later
    
    # set anything that we dont want to look at to np.nan
    time_ind, lat_ind, lon_ind = np.where( x == oob_value )
    x[ ..., lat_ind, lon_ind ] = np.nan
    
    # set out-of-bounds values -- grabbing any other lurking ones
    time_ind, lat_ind, lon_ind = np.where( np.isnan( x ) )
    dof[ (lat_ind, lon_ind) ] = np.nan
    dot[ (lat_ind, lon_ind) ] = np.nan

    # FIND THE NEGATIVE VALUES IN THE FULL NDARRAY
    def where_less_zero( arr ):
        ''' where in 1d slices of 3d array are less than zero '''
        new_arr = np.copy( arr )
        new_arr[:] = np.nan
        ind, = np.where( arr < 0 )
        new_arr[ range( len(ind) ) ] = ind
        return new_arr

    new_ind = np.apply_along_axis( where_less_zero, axis=0, arr=s )
    new_ind_plus1 = new_ind + 1

    # get index lengths of values along axis=0
    def count_ind( arr ):
        ''' number of ~np.nan values '''
        return len( arr[ ~np.isnan( arr ) ] )

    # get the counts of values along the time dimension
    new_ind_counts = np.apply_along_axis( count_ind, axis=0, arr=new_ind )
    new_ind_plus1_counts = np.apply_along_axis( count_ind, axis=0, arr=new_ind_plus1 )
    ind_counts = new_ind_counts + new_ind_plus1_counts
    # above takes aboout 1.5 mins
    
    # # not sure if these are useful    
    # s1_above = np.sum( s1a, axis=0 ) > 0
    # s1_below = np.sum( s1a, axis=0 ) < 0
    # # # # 
    
    # # IS THE FIRST ELEM OF THE NEW ARRAY > 0
    def s1_n0_greater0( arr ):
        return arr[ 0 ] > 0

    s1_n0_greater = np.apply_along_axis( s1_n0_greater0, axis=0, arr=s1 )

    # no transitions: all positive temps means no freeze day
    cur_ind = np.where( (ind_counts == 0) & (s1_n0_greater == True) )
    dot[ cur_ind ] = 0 # set dot to zero
    dof[ cur_ind ] = 365 # set dof to 365
    # case = 2

    # no transitions: all negative temps means no thaw day
    cur_ind = np.where( (ind_counts == 0) & (s1_n0_greater == False) )
    dot[ cur_ind ] = 365 # set dot to 365
    dof[ cur_ind ] = 0 # set dof to zero    
    # case = 3

    # # # #END CURRENT

    # [ML FIXED]
    # only one transition during the year, thawing or freezing
    cur_ind = np.where( (ind_counts == 2) & (x[0,...] < 0) )
    
# # # # # # # # # # #     
# ONCE WE FIGURE OUT THIS PART THE REST IS REALLY THE SAME THING WITH SLIGHT MODS.
    # THIS PART IS CONFUSING ME in ROUND 3 here...
    dot = 15 + 30 * ((ind[0]+1)-1) - np.round( x[ ind[0] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
    dot = 15 + 30 * ((ind[0, cur_ind[0], cur_ind[1]]+1)-1) - np.round( x[ ind[0] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
    dof = dof[ cur_ind ] = 350 # 350: we know the ground freezes so we use 350 rather than the special 365
    # case = 4
# # # # # # # # # # # 

        # places where we know the ground freezes and thaws, 
        #  but during a specific 12 months we just don't happen to witness both
        # only thaw occurs
        # if x[ ind[0] ] < 0:
            # [ml] note:((ind[0]+1)-1) is ind[0]+1 is the month number and minus 1 is to get to previous month
            #      we could make that a call to a months array -- months = range(1, 12+1)

        # [ round3 ] THIS GUY IS THE SAME AS THE ONE WE NEED TO FIGURE OUT ABOVE JUST FLIPPED!
        # only freeze occurs
        if x[ ind[0] ] > 0: 
            dof = 350 - 30 * (12-ind[1]-1) - np.round( x[ ind[1] ] / (np.diff( x[ ind[:2] ] ) / 30.0), decimals=0 )
            dot = np.array([15]) # 15: we know the ground thaws so we use 15 rather than the special 0
            grow = dof - dot
            case = 5

    # [ ! ON THIS HERE ! ]
    # two transitions occur: thaw, then freeze (this is the ideal case; everything else is an idiosyncratic edge case)
    cur_lat, cur_lon = np.where( (ind_counts == 4) & (x[0,...] < 0) )
    # THIS IS NON-WORKING!

    # [ml] note:((ind[0]+1)-1) is ind[0]+1 is the month number and minus 1 is to get to previous month
    #      we could make that a call to a months array -- months = range(1, 12+1)
    new_ind_cur = new_ind[ ..., cur_lat, cur_lon ]
    dot[ cur_ind ] = 15 + 30 * ((new_ind[0, cur_ind[0], cur_ind[1]]+1)-1) - np.round( x[ new_ind[0]:, cur_ind[0], cur_ind[1] ] / (np.diff( x[ ind[:2, cur_ind[0], cur_ind[1]] ] ) / 30.0), decimals=0 )
    
    dof[ cur_ind ] = 350 - 30 * (12-ind[3]-1) - np.round( x[ ind[3] ] / (np.diff( x[ ind[2:4] ] ) / 30.0), decimals=0 )
    grow = dof - dot
    case = 0

    
    # [ML FIXED]
    cur_ind = np.where( (ind_counts == 4) & (x[0,...] > 0) )

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
        ind2 = np.sort( np.concatenate( (ind2, np.array([0])) ) )

        diffed = np.diff( ind2 )
        # [ml] m1, m2 are month indexes
        m1, = np.where( diffed == np.max( diffed ) )#[0][-1]
        m1 = m1 + 1
        m2 = np.where( np.delete(np.diff( ind2 ), (m1-1)-1) == max( np.delete(np.diff( ind2 ), (m1-1)-1)) )[-1] + 1

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
    return 'dof:{} dot:{} logs:{}'.format( dof,dot,grow )


if __name__ == '__main__':
    # get the GD filenames here 
    # filenames = []
    # set up the array up here...
    arr = np.array([ read_arr( fn )[0] for fn in filenames ])
    x = arr.copy()

    # INSERT THIS TEST DATA TO PLACES IN THE ARRAY
    test_vals = [np.array([-16, -5, -1, 3, 5, 10, 12, 16, 11, -3, -15, -16]),
        np.array([-16, -5, -1, 3, 5, 10, 12, 16, 11, np.nan, -15, -16]),
        np.array([1, 3, 4, 6, 7, 12, 15, 12, 8, 9, 4, 2]),
        np.array([-16, -15, -13, -11, -10, -5, 0, -2, -4, -12, -13, -16]),
        np.array([-16, -13, -8, -6, 1, 4, 7, 11, 8, 4, 2, 1]),
        np.array([1, 3, 1, 5, 8, 10, 14, 11, 7, -2, -5, -2]),
        np.array([1, 3, 1, 5, 8, 10, 14, 10, 4, -1, -4, 1]),
        np.array([1, -5, -4, -2, 3, 5, 10, 8, 6, 4, 4, 1]),
        np.array([-11, 1, -7, -3, 2, 6, 11, 10, 8, -1, -5, -10]) ]

    ind = [ (800, 1000+idx) for idx in range( len( test_vals ) ) ]

    for idx, a in zip( ind, test_vals ):
        x[ ..., idx[0], idx[1] ] = a

    # how to get to the data we added for testing later on
    x[ ..., 800, 1000:np.array(ind)[:,1].max() ]

    




