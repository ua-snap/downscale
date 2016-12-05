# TESTING SCRIPT FOR COMPARISON WITH PYTHON VERSION OF THIS FUNCTION.

# values: valid day of freeze or thaw values are 15:350.
# 0 and 365 are special, indicating always frozen (freeze: 0, thaw: 365, grow: 0) or always thawed (freeze: 365, thaw: 0, grow: 365)

# This function operates on a vector x of 12 monthly average temperatures
tfg_days <- function(x){
    if(any(x==0, na.rm=T)) x[x == 0] <- -0.0001 # need to treat zero as freezing (working with signs)
    s1 <- sign(x) # positive or negative monthly temps
    s <- s1[1:11]*s1[2:12] # products of consecutive months' signs: positive indicates no change; negative indicates a potential freeze or thaw transition
    ind <- sort(c(which(s < 0), which(s < 0) + 1)) # may be length zero (no transitions)

    if(any(is.na(x))){ # ignore cells with missing data
        dot <- dof <- grow <- NA
        case <- 1
    } else if(length(ind)==0 & s1[1] > 0){ # no transitions: all positive temps means no freeze day
        dot <- 0
        dof <- grow <- 365
        case <- 2
        print(case)
    } else if(length(ind)==0 & s1[1] < 0) { # no transitions: all negative temps means no thaw day
        dot <- 365
        dof <- grow <- 0
        case <- 3
        print(case)
    } else if(length(ind)==2) { # only one transition during the year, thawing or freezing
        # places where we know the ground freezes and thaws, but during a specific 12 months we just don't happen to witness both
        if(x[ind[1]] < 0){ # only thaw occurs
            dot <- 15 + 30*(ind[1]-1) - round(x[ind[1]] / (diff(x[ind[1:2]]) / 30), 0)
            dof <- 350 # 350: we know the ground freezes so we use 350 rather than the special 365
            grow <- dof - dot
            case <- 4
            print(case)
        }
        if(x[ind[1]] > 0){ # only freeze occurs
            dof <- 350 - 30*(12-ind[2]) - round(x[ind[2]] / (diff(x[ind[1:2]]) / 30), 0)
            dot <- 15 # 15: we know the ground thaws so we use 15 rather than the special 0
            grow <- dof - dot
            case <- 5
            print(case)
        }
    } else if(length(ind)==4 & s1[1] < 0) { # two transitions occur: thaw, then freeze (this is the ideal case; everything else is an idiosyncratic edge case)
        # ind is the index where 
        dot <- 15 + 30*(ind[1]-1) - round( x[ind[1]] / (diff(x[ind[1:2]]) / 30), 0)
        dof <- 350 - 30*(12-ind[4]) - round(x[ind[4]] / (diff(x[ind[3:4]]) / 30), 0)
        grow <- dof - dot
        case <- 0
        print(case)
    } else if(length(ind)==4 & s1[1] > 0){  # two transitions occur but backward to what is expected; freeze, then thaw
        if(ind[1] >= 7) { # freeze occurs in second half of year as expected; late thaw is spurious
            dof <- 350 - 30*(12-ind[2]) - round(x[ind[2]] / (diff(x[ind[1:2]]) / 30), 0)
            dot <- 15 # ignore spurious post-freeze thaw; treat as early, unobserved thaw
            grow <- dof - dot
            case <- 6
            print(case)
        }
        if(ind[1] <= 6){ # spurious freeze occurs in first half of year; thaw probably fine
            dot <- 15 + 30*(ind[3]-1) - round(x[ind[3]] / (diff(x[ind[3:4]]) / 30), 0)
            dof <- 350 # ignore spurious early freeze; treat as late, unobserved freeze
            grow <- dof - dot
            case <- 7
            print(case)
        }
    } else if(length(ind) > 4){ # more than two transitions; at least one definitely spurious
        ind2 <- c(0, which(s < 0))
        m1 <- tail(which(diff(ind2) == max(diff(ind2))), 1) + 1
        m2 <- tail(which(diff(ind2)[-(m1-1)] == max(diff(ind2)[-(m1-1)])), 1) + 1
        if(m1==m2) m2 <- m2-1
        ind2 <- ind2[sort(c(m1, m2))]
        ind <- sort(c(ind2, ind2+1))
        dot <- 15 + 30*(ind[1]-1) - round(x[ind[1]] / (diff(x[ind[1:2]]) / 30), 0)
        dof <- 350 - 30*(12-ind[4]) - round(x[ind[4]] / (diff(x[ind[3:4]]) / 30), 0)
        grow <- dof - dot
        case <- 8
        print(case)
    } else {
        dot <- dof <- grow <- NA
        # stop(print(paste("Condition unconsidered:", paste(x, collapse=" "))))
    }
    return( paste(' dof',dof, ' dot',dot, ' grow',grow, sep=':' ) )
}

x <- c(-6.1999998, -9.1999998, -3.0999999, 1.4, 7.3000002, 11.0, 13.6, 13.6, 10.5, 1.3, -4.0999999, -7.0999999)

# test data
x <- list(
  c(-16, -5, -1, 3, 5, 10, 12, 16, 11, -3, -15, -16),
  c(-16, -5, -1, 3, 5, 10, 12, 16, 11, NA, -15, -16),
  c(1, 3, 4, 6, 7, 12, 15, 12, 8, 9, 4, 2),
  c(-16, -15, -13, -11, -10, -5, 0, -2, -4, -12, -13, -16),
  c(-16, -13, -8, -6, 1, 4, 7, 11, 8, 4, 2, 1),
  c(1, 3, 1, 5, 8, 10, 14, 11, 7, -2, -5, -2),
  c(1, 3, 1, 5, 8, 10, 14, 10, 4, -1, -4, 1),
  c(1, -5, -4, -2, 3, 5, 10, 8, 6, 4, 4, 1),
  c(-11, 1, -7, -3, 2, 6, 11, 10, 8, -1, -5, -10)
)

# run it
for( i in data ){ 
    print( tfg_days( i ) ) 
}


# # # # OUTPUTS:
# [1] " dof:284: dot:83: grow:201"
# [1] " dof:NA: dot:NA: grow:NA"
# [1] " dof:365: dot:0: grow:365"
# [1] " dof:0: dot:365: grow:0"
# [1] " dof:350: dot:131: grow:219"
# [1] " dof:283: dot:15: grow:268"
# [1] " dof:284: dot:15: grow:269"
# [1] " dof:350: dot:117: grow:233"
# [1] " dof:287: dot:123: grow:164"


