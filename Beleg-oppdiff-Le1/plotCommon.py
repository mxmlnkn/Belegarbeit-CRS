#!/usr/bin/python
# -*- coding: utf-8 -*-

Zstanal   = 0.0545
vSelected = [ 0.05, 0.1, 0.15, 0.2, 0.25, 0.5 ]
colors    = [ 'r', 'g', 'b', 'orange', 'm', 'c', 'k', 'y', 'dimgray' ]
figTiny   = [4,3]   # used especially for Y_i plots, because 3 per line in pdf
figSmall  = [5,3]   # plots which are 2 per line in pdf or which just don't have much details
#figNormal = [6,4]   # for 1 per line plot
#figLarge  = [8,5]   # One per line with many details

def bisectionExtrema( f,a,b,nIterations=16,debug=False ):
    """
    " finds the extremum of f(x) in the interval (a,b)
    " assumes that b > a and that only one extremum exists in [a,b]
    """
    """
    " ^
    " |           .´
    " |  .'',   .'
    " |.'    ´'`
    " +--------------------
    "  a  c  b
    """
    extremaWasInside = True
    for i in range(nIterations):
        assert( b > a )
        c  = 0.5 *(a+b)
        # everything smaller than interval width / 6 should basically be enough
        # if the factor is too small it may result in problems like the
        # bisection quitting too early, thereby giving back wrong maxima!
        dx = 1e-2*(b-a)
        # if floating point precision exhausted then these differences will be 0
        # for a and b do only onesided, because else we would leave the given interval!
        left   = f(a+dx) > f(a   )
        middle = f(c+dx) > f(c-dx)
        right  = f(b   ) > f(b-dx)
        if left == middle and middle == right and i == 0:
            extremaWasInside = False

        if debug:
            print "f at x=",a,"going up?",left  ," ( f(x+dx)=",f(a+dx),", f(x   =",f(a   )
            print "f at x=",c,"going up?",middle," ( f(x+dx)=",f(c+dx),", f(x-dx)=",f(c-dx)
            print "f at x=",b,"going up?",right ," ( f(x   )=",f(b   ),", f(x-dx)=",f(b-dx)

        # If the sign of the derivative is the same for all, then
        # the maximum is not inside the specified interval!
        if ( left == middle and middle == right ):
            if extremaWasInside:
                break   # this can also happen if dx is too small to resolve, so we break the search
            else:
                return None, None, None # unreasonable result, therefore error code
        elif left == middle:
            a = c
        elif middle == right:
            b = c
        else:
            # This happens if there are two extrema inside interval
            return None, None, None

    c = 0.5*(a+b)
    return f(c), c, 0.5*(b-a)

def axisIsLog( ax, axis ):
    if axis == 'x':
        return ax.get_xscale() == 'log'
    elif axis == 'y':
        return ax.get_yscale() == 'log'
    else:
        assert( False, "axis neither 'x' nor 'y'!" )


def axisMin( ax, axis ):
    xmin=float('+inf')
    isLog = axisIsLog( ax, axis )
    for line in ax.get_lines():
        if axis == 'x':
            x = line.get_xdata()
        else:
            x = line.get_ydata()
        if isLog:
            x = x[ x>0 ]
        xmin = min( xmin, min(x) )
    return xmin

def axisMax( ax, axis ):
    xmax=float('-inf')
    isLog = axisIsLog( ax, axis )
    for line in ax.get_lines():
        if axis == 'x':
            x = line.get_xdata()
        else:
            x = line.get_ydata()
        if isLog:
            x = x[ x>0 ]
        xmax = max( xmax, max(x) )
    return xmax

def autoRange( ax, axis, lb, rb = None ):
    if rb == None:
        rb = lb
    isLog = axisIsLog( ax, axis )

    xmin=axisMin( ax, axis )
    xmax=axisMax( ax, axis )

    from math import log,exp

    if isLog:
        dx   = log(xmax) - log(xmin)
        xmin /= exp( lb*( dx ) )
        xmax *= exp( rb*( dx ) )
    else:
        dx = xmax - xmin
        xmin -= lb*dx
        xmax += rb*dx

    if axis == 'x':
        ax.set_xlim( [xmin,xmax] )
    else:
        ax.set_ylim( [xmin,xmax] )

def autoRangeXY( ax, lb = 0.1, rb = None, bb = None, tb = None ):
    if rb == None:
        rb = lb
    if tb == None:
        tb = lb
    if bb == None:
        bb = lb

    autoRange( ax, 'x', lb, rb )
    autoRange( ax, 'y', bb, tb )

from math import log10,ceil,floor
def autoLabel( ax, axis, nbins=5, roundFunc=ceil ):
    #xmin  = axisMin( ax, axis )
    #xmax  = axisMax( ax, axis )
    #isLog = axisIsLog( ax, axis )
    #if isLog:
    #    dx   = roundFunc( ( log10(xmax) - log10(xmin) ) / nbins )
    #else:
    #    assert( False, "Not yet implemented" )
    #
    #from numpy import arange
    #n0 = int( floor( log10( xmin ) ) )
    #n1 = int( ceil ( log10( xmax ) ) )
    ##print "n0 =",n0,", n1 =",n1,", dx =",dx
    #xpos = 10.**( n0 + arange(nbins+2)*dx )
    #ax.set_xticks( xpos )
    ##print "set xlabels at : ", xpos
    print ""


def relErr( x, y ):
    from numpy import zeros
    assert( len(x) == len(y) )
    non0 = abs(y) > 1e-16
    #non0 = abs(y) != 0
    tmp = ( x[non0] - y[non0] ) / y[non0]
    res = zeros( len(y) )
    res[non0] = tmp
    return y[non0], abs(res[non0])

def finishPlot( fig, ax, fname, loc='best' ):
    l = ax.legend( loc=loc, prop={'size':10}, labelspacing=0.2,
                   fancybox=True, framealpha=0.5 )
    #if l != None:
    #    l.set_zorder(0)  # alternative to transparency
    fig.tight_layout()
    fig.savefig( fname+".pdf" )
    print "[Saved '"+fname+".pdf']"
    from matplotlib.pyplot import close
    close( fig )

def p(x0,y0,x):
    """
    This function approximates f so that f(x0)=y0 and returns f(x)
    """
    assert( len(x0) == len(y0) )
    # number of values. Approximating polynomial is of degree n-1
    n = len(x0)
    assert( n > 1 )
    res = 0
    for k in range(n):
        prod = y0[k]
        for j in range(n):
            if j != k:
                prod *= (x0[j]-x)/(x0[j]-x0[k])
        res += prod
    return res

def readUlfFile( name ): # name without extension
    # read data from file (some times we want 0.10, but have v=0.1, that's what this try is for
    try:
        path = name+'.ulf'
        f    = open( path,'r')
    except IOError:
        try:
            path = name+'0.ulf'
            f    = open( path,'r')
        except IOError:
            try:
                path = name+'00.ulf'
                f    = open( path,'r')
            except IOError:
                path = name
                f    = open( path,'r')
    from numpy import array,loadtxt,arange
    header = array( f.readline().split() )
    data   = loadtxt( path,skiprows=1 )
    # create dictionary, which returns column number, if we input column header string
    hdict  = dict( zip( header, arange(len(header)) ) )

    return data,hdict

def calcZ( YCH4, YO2 ):
    YO20    = 0.23
    YCH40   = 1.0
    nu      = 3.99
    return ( nu*YCH4 - YO2 + YO20 )/( nu*YCH40 + YO20 )

def calcZBilger( data, hdict, debug=False ):
    from numpy import array,zeros,ones,abs,sum
    ########################### Elementmassenbrueche ###########################

    # .ulf 1st line needs to be in the same order !!!:
    # CH4 O2 H2O H2O2 CO2 CH3 HO2 CH2O HCO CH3O CO OH H O H2 AR N2
    stoffnamen = [ "CH4", "O2", "H2O", "H2O2", "CO2", "CH3", "HO2", "CH2O", "HCO",
                   "CH3O", "CO", "OH", "H", "O", "H2", "AR", "N2" ]
    Nstoffe    = len(stoffnamen)
    elemnamen  = [ "C", "H", "O", "Ar", "N" ]
    Nelem      = len(elemnamen)
    i0         = hdict["CH4"]
    Ystoffe    = data[ :, i0:(i0+len(stoffnamen)) ].transpose()
    Nx         = len(Ystoffe[0,:])
    Melem      = array( [ 12.0107, 1.00794, 15.9994, 39.948, 14.0067 ] )
    a = array( [
        1,4,0,0,0,   # CH4
        0,0,2,0,0,   # O2
        0,2,1,0,0,   # H20
        0,2,2,0,0,   # H2O2
        1,0,2,0,0,   # CO2
        1,3,0,0,0,   # CH3
        0,1,2,0,0,   # HO2
        1,2,1,0,0,   # CH2O
        1,1,1,0,0,   # HCO
        1,3,1,0,0,   # CH3O
        1,0,1,0,0,   # CO
        0,1,1,0,0,   # OH
        0,1,0,0,0,   # H
        0,0,1,0,0,   # O
        0,2,0,0,0,   # H2
        0,0,0,1,0,   # Ar
        0,0,0,0,2    # N2
    ] )
    a = a.reshape( [Nstoffe,Nelem] )
    # reshape((3, 2) -> [ [0, 1],[2, 3],[4, 5] ]
    # Calculate molar mass of compounds from the molar masses of the elements they consist of
    Mstoffe = zeros( len(stoffnamen) )
    for i in range(Nstoffe):
        Mstoffe[i] = sum(a[i]*Melem)

    if debug:
        print "Stoffnamen: ",len(stoffnamen)
        print "elemnamen : ",len(elemnamen)
        print "a         : ",len(stoffnamen)
        itest = 2
        print stoffnamen[itest]," consists of ",
        for j in range(Nelem):
            print a[itest,j]," ",elemnamen[j],"-, "
        print "atoms and has a molar mass of ",Mstoffe[itest]," g/mol"

    # Calculate element mixture fractions
    Z = zeros( [Nelem,Nx] )
    for j in range(Nelem):
        for i in range(Nstoffe):
            Z[j] += a[i,j]*Melem[j] / Mstoffe[i] * Ystoffe[i,:]
        if debug:
            print "Z_",elemnamen[j]," = ",Z[j]
    # Note: access is Z[ielem,ix] (for that the data was transposed when writing to Ystoffe!)
    #       sum only over all elements (for every position ix) -> should get array
    #       of length Nx with every element being ~1.0

    if debug:
        print "Nx = ",Nx
        print "sum_j Z_j = ",sum( Z,axis=0 )
        print len(Z),len(Z[0])
    meanerr = abs( sum( Z,axis=0 )-ones(Nx) ).sum() / (1.0*Nx)
    if debug or meanerr >= 1e-2:
        print "mean error = ",meanerr
    assert( meanerr < 1e-2 )

    ############## Elementmischungsbruch / Bilger-Mischungsbruch ###############

    YCH4B   = 1.0
    YO2B    = 0.0
    YCH4Ox  = 0.0
    YO2Ox   = 0.23
    nuCH4   = 1.0
    nuO2    =-2.0
    MCH4    = Mstoffe[ 0 ]
    MO2     = Mstoffe[ hdict["O2"]-hdict["CH4"] ]
    betaB   = 2.0*( YCH4B  / ( nuCH4 * MCH4 ) + YO2B  / ( nuO2 * MO2 ) )
    betaOx  = 2.0*( YCH4Ox / ( nuCH4 * MCH4 ) + YO2Ox / ( nuO2 * MO2 ) )
    acoeff  = 1.0
    bcoeff  = 4.0
    jC      = 0
    jH      = 1
    jO      = 2
    beta    = Z[jC] / ( acoeff*nuCH4*Melem[jC] ) + Z[jH] / ( bcoeff*nuCH4*Melem[jH] ) + Z[jO] / ( nuO2*Melem[jO] )
    ZBilger = ( beta-betaOx )/( betaB-betaOx )

    if debug:
        print "M_O2 = ",MO2,", M_CH4 = ",MCH4,", M_O = ",Melem[jO],", M_C = ",Melem[jC]
        print "beta_B  = ",betaB  ,", beta_Ox  = ",betaOx
        print "beta[0] = ",beta[0],", beta[-1] = ",beta[-1]
        print len(x),len(ZBilger)

    return ZBilger

def calcChi(data,hdict):
    from numpy import max,zeros
    if "chi" in hdict and max( data[:,hdict["chi"]] ) != 0:
        return data[:,hdict["chi"]]
    x       = data[:,hdict["X"]]
    Z       = data[:,hdict["Z"]]
    Nx      = len(x)
    # Calculate chi = - nabla Z_Bilger
    dZdx    = zeros(Nx)
    for i in range(Nx-1):
        dZdx[i] = (Z[i+1]-Z[i])/(x[i+1]-x[i])  # FDS
    dZdx[-1] = (Z[-1]-Z[-2])/(x[-1]-x[-2])     # BDS on last point
    D        = data[:,hdict["lambda"]] / ( data[:,hdict["cpMean"]] * data[:,hdict["rho"]] )
    return 2*D*(dZdx)**2

def precalcOppdiffValues( basedir ):
    v_list          = []
    Tmax_list       = []
    Tst_list        = []
    chist_list      = []
    chist_sr        = []
    chist_intp_list = []  # chi_stoic calculated using interpolation, therefore slightly more correct!
    strainrate_list = []
    PVst_list       = []
    PVmax_list      = []

    #basedir = './results'
    # analyse all files in 'results/v_*.ulf' and append T_max and chi_st to lists
    import os
    for filename in os.listdir(basedir):
        path = basedir+'/'+filename
        if not os.path.isfile(path):  # ignore directories
            continue

        # extract number between s0 and s1, meaning extract 0.30 from v_0.30.ulf
        s0 = 'v_'
        s1 = '.ulf'
        v  = float( path[ path.find(s0)+len(s0) : path.find(s1) ] )
        v_list.append( v )

        data, hdict = readUlfFile( path )
        T    = data[:,hdict["T"]]
        imax = T.argmax()
        Tmax_list.append( T[imax] )

        ZUlf    = data[:,hdict["Z"]]
        ist     = abs( ZUlf - Zstanal ).argmin()
        chi     = calcChi(data,hdict)

        chist   = chi[ist] # i-stoch calculated above with Z_st = 0.0545
        chist_list.append( chist )
        Tst_list.append( T[ist] )
        #chist_sr.append( r"$\chi_\mathrm{st}=" + format( chist, '.3f' ) + r"\,s^{-1}$" )
        print ""
        print "v =",v
        print "Chi around >chi_st<:", chi[ist-1],",>", chi[ist],"<,", chi[ist+1]
        print "Z   around >Z_st<  :",ZUlf[ist-1],",>",ZUlf[ist],"<,",ZUlf[ist+1]

        Zl = data[ist-1,hdict["Z"]]
        Zc = data[ist  ,hdict["Z"]]
        Zr = data[ist+1,hdict["Z"]]
        chist_intp = chi[ist-1] * (Zc-Zstanal)*(Zr-Zstanal)/( (Zc-Zl)*(Zr-Zl) ) + \
                     chi[ist  ] * (Zl-Zstanal)*(Zr-Zstanal)/( (Zl-Zc)*(Zr-Zc) ) + \
                     chi[ist+1] * (Zl-Zstanal)*(Zc-Zstanal)/( (Zl-Zr)*(Zc-Zr) )
        chist_intp_list.append( chist_intp )
        print "Chist interpolated :",chist_intp
        chist_sr.append( r"$\chi_\mathrm{st}=" + format( chist_intp, '.3f' ) + r"\,\mathrm{s}^{-1}$" )

        from math import sqrt
        # Calculate strain rate
        # k = 2/L*( v1*sqrt(rho1/rho2) + v2 )
        # we used v1 = -v2 for our setup, simplifying things drastically:
        # k = 2/L*( sqrt(rho1/rho2) - 1)*v ~ v
        L     = 0.02  # m <- can be seen with AXIS_MIN and AXIS_MAX in ULF settings
        rhoB  = 0.656 # kg/m^3 for CH4
        rhoOx = 1.16  # kg/m^3 for air
        k     = 2./L*abs( 1.0-sqrt(rhoB/rhoOx) )*v
        strainrate_list.append( k )
        print "strain rate k     :",k

        # Fortschrittsvariable
        PV = data[:,hdict["CO"]] + data[:,hdict["CO2"]]
        PVst_list.append( PV[ist] )
        PVmax_list.append( PV[imax] )

    from numpy import array
    return array( v_list          ), \
           array( Tmax_list       ), \
           array( Tst_list        ), \
           array( chist_list      ), \
           array( chist_sr        ), \
           array( chist_intp_list ), \
           array( strainrate_list ), \
           array( PVst_list       ), \
           array( PVmax_list      )

def precalcFlameletValues( basedir ):
    fname_list      = []
    chist_list      = []
    chist_sr        = []
    chist_calc_list = []
    chist_intp_list = []
    Tmax_list       = []
    Tst_list        = []
    PVst_list       = []
    PVmax_list      = []

    #basedir = './results'
    # analyse all files in 'results/chist_*.ulf' and append T_max and chi_st to lists
    import os
    for filename in os.listdir(basedir):
        path = basedir+'/'+filename
        if not os.path.isfile(path):  # ignore directories
            continue

        fname_list.append(path)

        # extract number between s0 and s1, meaning extract 0.30 from chist_0.30.ulf
        s0 = 'chist_'
        s1 = '.ulf'
        chist  = float( path[ path.find(s0)+len(s0) : path.find(s1) ] )
        chist_list.append( chist )
        chist_sr.append( r"$\chi_\mathrm{st}=" + format( chist, '.3f' ) + r"\,\mathrm{s}^{-1}$" )

        data, hdict = readUlfFile( path )
        T    = data[:,hdict["T"]]
        imax = T.argmax()
        Tmax_list.append( T[imax] )

        # Calculate Chi_stoich and compare to those input into the simulation
        ZUlf       = data[:,hdict["Z"]]
        ist        = abs( ZUlf - Zstanal ).argmin()
        chi        = calcChi(data,hdict)
        chist_calc = chi[ist]
        chist_calc_list.append( chist_calc )
        Tst_list.append( T[ist] )

        Zl = ZUlf[ist-1]
        Zc = ZUlf[ist]
        Zr = ZUlf[ist+1]
        chist_intp = chi[ist-1] * (Zc-Zstanal)*(Zr-Zstanal)/( (Zc-Zl)*(Zr-Zl) ) + \
                     chi[ist  ] * (Zl-Zstanal)*(Zr-Zstanal)/( (Zl-Zc)*(Zr-Zc) ) + \
                     chi[ist+1] * (Zl-Zstanal)*(Zc-Zstanal)/( (Zl-Zr)*(Zc-Zr) )
        print "Test lagrange: ",chist_intp," ?= ",p( data[ist-1:ist+2,hdict["Z"]], chi[ist-1:ist+2], Zstanal )
        chist_intp2 = p( data[ist-2:ist+2,hdict["Z"]], chi[ist-2:ist+2], Zstanal )

        chist_intp_list.append( chist_intp )

        print "Chist configured              : ",chist
        print "Chist calculated naively      : ",chist_calc
        print "Chist interpolated (3 points) : ",chist_intp
        print "Chist interpolated (5 points) : ",chist_intp2
        print "Chist def p                   : ",p( ZUlf[ ist-1:ist+2 ], chi[ ist-1:ist+2 ], Zstanal )
        print "x: ZUlf = ", ZUlf[ ist-1:ist+2 ]
        print "y: T    = ", T   [ ist-1:ist+2 ]
        Tintp = lambda Z: p( ZUlf[ ist-1:ist+2 ], T[ ist-1:ist+2 ], Z )
        print "polynomial : ",Tintp( ZUlf[ist-1] ), Tintp( ZUlf[ist+1] )
        Tmax, Zmax, dZ = bisectionExtrema( Tintp , ZUlf[ist-2], ZUlf[ist+2], debug=True )
        print "Z borders                     : ",ZUlf[ist-2], ZUlf[ist+12]
        print "Zmax def p                    : ",Zmax," +- ",dZ
        print "Tmax def p                    : ",Tmax
        print "Tst def p                     : ",p( ZUlf[ ist-1:ist+2 ], T[ ist-1:ist+2 ], Zstanal )

        from matplotlib.pyplot import figure, plot, show
        from numpy import linspace
        figure()
        Z = linspace( ZUlf[ist-2], ZUlf[ist+2], 100, endpoint=True )
        plot( Z, Tintp(Z) )
        plot( ZUlf[ ist-2:ist+3 ], T[ ist-2:ist+3 ], 'ro' )
        show()

        exit()

        # Fortschrittsvariable
        PV = data[:,hdict["CO"]] + data[:,hdict["CO2"]]
        PVst_list.append( PV[ist] )
        PVmax_list.append( PV[imax] )

    from numpy import array
    return array( fname_list      ), \
           array( chist_list      ), \
           array( chist_sr        ), \
           array( chist_calc_list ), \
           array( chist_intp_list ), \
           array( Tmax_list       ), \
           array( Tst_list        ), \
           array( PVst_list       ), \
           array( PVmax_list      )
