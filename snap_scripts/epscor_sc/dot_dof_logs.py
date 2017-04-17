# # # # 
# PYTHON VERSION OF MATT's DOF/DOT/LOGS
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
    return 'dof:{} dot:{} logs:{}'.format( dof,dot,grow )


if __name__ == '__main__':
    import numpy as np
    # test data
    x_list = [ np.array([-16, -5, -1, 3, 5, 10, 12, 16, 11, -3, -15, -16]),
          np.array([-16, -5, -1, 3, 5, 10, 12, 16, 11, np.nan, -15, -16]),
          np.array([1, 3, 4, 6, 7, 12, 15, 12, 8, 9, 4, 2]),
          np.array([-16, -15, -13, -11, -10, -5, 0, -2, -4, -12, -13, -16]),
          np.array([-16, -13, -8, -6, 1, 4, 7, 11, 8, 4, 2, 1]),
          np.array([1, 3, 1, 5, 8, 10, 14, 11, 7, -2, -5, -2]),
          np.array([1, 3, 1, 5, 8, 10, 14, 10, 4, -1, -4, 1]),
          np.array([1, -5, -4, -2, 3, 5, 10, 8, 6, 4, 4, 1]),
          np.array([-11, 1, -7, -3, 2, 6, 11, 10, 8, -1, -5, -10]) ]

    for ct, x in enumerate( x_list ):
        print( '{} - {}'.format( ct, tfg_days( x ) ) )


# # # # OUTPUTS:
# 0 - dof:[ 284.] dot:[ 83.] logs:[ 201.]
# 1 - dof:[ nan] dot:[ nan] logs:[ nan]
# 2 - dof:[365] dot:0 logs:[365]
# 3 - dof:[0] dot:[365] logs:[0]
# 4 - dof:[350] dot:[ 131.] logs:[ 219.]
# 5 - dof:[ 283.] dot:[15] logs:[ 268.]
# 6 - dof:[ 284.] dot:[15] logs:[ 269.]
# 7 - dof:[350] dot:[ 117.] logs:[ 233.]
# 8 - dof:[ 287.] dot:[ 123.] logs:[ 164.]


