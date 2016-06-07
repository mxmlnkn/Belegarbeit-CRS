#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os

v = []; Tmax = []; Chi_st = []; T_st = []

basedir = './results'
for filename in os.listdir(basedir):
    path = basedir+'/'+filename
    if not os.path.isfile(path):
        continue

    s0='v_'; s1='.ulf'
    v.append( float( path[ path.find(s0)+len(s0) : path.find(s1) ] ) )

    f = open(path,'r')
    header = np.array( f.readline().split() )
    data   = np.loadtxt( path,skiprows=1 )

    # create dictionary, which returns column number, if we input column header string
    hdict  = dict( zip( header, np.arange(len(header)) ) )

    Tmax.append  ( data[:,hdict["T"]].max() )
    
    x       = data[:,hdict["X"]]
    Zm      = data[:,hdict["Z"]]
    Nx      = len(x)
    Zstanal = 0.0545
    ist     = abs(Zm-Zstanal).argmin()  
    dZdx    = np.zeros(Nx)
    for i in range(Nx-1):
        dZdx[i] = (Zm[i+1]-Zm[i])/(x[i+1]-x[i])  # FDS
    dZdx[-1] = (Zm[-1]-Zm[-2])/(x[-1]-x[-2])  # BDS on last point
    D        = data[:,hdict["lambda"]] / ( data[:,hdict["cpMean"]] * data[:,hdict["rho"]] )
    chi      = 2*D*(dZdx)**2
    chist    = chi[ist] # i-stoch calculated above with Z_st = 0.0545

    Chi_st.append( chist )
    T_st.append( data[ ist, hdict["T"] ] )

    # k = 2/L*( v1*sqrt(rho1/rho2) + v2 )
    # we used v1 = -v2 for our setup, simplifying things drastically:
    # k = 2/L*sqrt(rho1/rho2) * 2*v ~ v
    # for gases the density is almost the same
    # L = 0.02 m <- can be seen with AXIS_MIN and AXIS_MAX in ulf settings
    #k.append( 4.*L )

outdata = []
outdata.append( [v, Tmax, Chi_st] )
# numpy version too old on node (1.4.2)
#np.savetxt( fname=path, X=outdata, header="#  v/(m/s)   Tmax/K    Chi_st    PV" )

plt.figure()
plt.xlabel(r"$v$")
plt.ylabel(r"$T_\mathrm{max}$")
plt.plot( v, Tmax, 'o' )

plt.figure()
plt.xlabel(r"$v$")
plt.ylabel(r"$\chi_\mathrm{stoch.}$")
plt.plot( v, Chi_st, 'o' )

plt.figure()
plt.xlabel(r"Chi")
plt.ylabel(r"T_max")
plt.plot( Chi_st, Tmax,'o' )

plt.figure()
plt.xlabel(r"Chi")
plt.ylabel(r"T_st")
plt.plot( Chi_st, T_st,'o' )

plt.show()

