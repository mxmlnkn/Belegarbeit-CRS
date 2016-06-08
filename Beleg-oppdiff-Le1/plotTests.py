#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os
from plotCommon import *

def plotSCurves():
    # needs: chist_list Tmax_list strainrate_list PVst_list chist_list PVmax_list v_list iSelected finishPlot
    fig, ax = [], []
    for i in range(7):
        fig.append( plt.figure( figsize=figSmall ) )
    l_chist = r"$\chi_\mathrm{st} / s^{-1}$"
    l_k     = r"Inverse Streckungsrate $k^{-1} / \mathrm{s}$"
    l_pvst  = r"$\mathrm{PV}_\mathrm{st} := Y_\mathrm{CO}( Z_\mathrm{st} ) + Y_{\mathrm{CO}_2}( Z_\mathrm{st} ) $"
    l_pvmax = r"$\mathrm{PV}_{T_\mathrm{max}} := \left[Y_\mathrm{CO}\left( Z \right) + Y_{\mathrm{CO}_2}\left( Z \right) \right]_{T(Z)=T_\mathrm{max}}$"
    xlabels = [ l_chist, l_chist, l_k, l_k, l_pvst, l_pvst, l_pvmax ]

    for i in range(len(fig)):
        ax.append( fig[i].add_subplot( 111, xlabel = xlabels[i],
                ylabel = r"$T_\mathrm{max} / ^\circ\mathrm{K}$",
                xmargin = 0.1, ymargin = 0.1 ) )
        #ax[-1].tick_params(axis='both', which='major', labelsize=10)
        #ax[-1].tick_params(axis='both', which='minor', labelsize=8 )

    # S-Curves
    burning = Tmax_list > 400
    ax[0].plot( chist_list                     , Tmax_list             , 'bo'  )
    ax[0].plot( chist_list        [ burning   ], Tmax_list[ burning   ], 'b-'  )
    ax[2].plot( 1./strainrate_list             , Tmax_list             , 'bo'  )
    ax[2].plot( 1./strainrate_list[ burning   ], Tmax_list[ burning   ], 'b-'  )
    ax[4].plot( PVst_list         [ burning   ], Tmax_list[ burning   ], 'bo'  )
    # zooms / selected configurations
    ax[1].plot( chist_list        [ iSelected ], Tmax_list[ iSelected ], 'bo-' )
    ax[3].plot( 1./strainrate_list[ iSelected ], Tmax_list[ iSelected ], 'bo-' )
    ax[5].plot( PVst_list         [ iSelected ], Tmax_list[ iSelected ], 'bo'  )
    ax[6].plot( PVmax_list        [ iSelected ], Tmax_list[ iSelected ], 'bo'  )

    ax[0].set_xlim( [ 0, 1.1*np.max( chist_list ) ] )
    labels = [ str(v) for v in v_list ]

    for i in range( len(v_list) ):
        # S curve Tmax over chi_st
        if v_list[i] > 0.2:
            ax[0].annotate( labels[i], xycoords='data', textcoords='offset points',
                size=9, xytext=(-8,5), xy=( chist_list[i]    , Tmax_list[i] ) )
        # S curve Tmax over strain rate
        if ( v_list[i] <= 0.1 or v_list[i] == 1.0 ) and not v_list[i] == 0.025:
            ax[2].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(-8,5), xy=( 1./strainrate_list[i], Tmax_list[i] ) )
        if v_list[i] == 0.025:
            ax[2].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(-12,5), xy=( 1./strainrate_list[i], Tmax_list[i] ) )

    for i in iSelected:
        ax[1].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(5,-3), xy=( chist_list[i]        , Tmax_list[i] ) )
        ax[3].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(5,-3), xy=( 1./strainrate_list[i], Tmax_list[i] ) )
        ax[5].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(5,-3), xy=( PVst_list[i]         , Tmax_list[i] ) )
        ax[6].annotate( labels[i], xycoords='data', textcoords='offset points',
            size=9, xytext=(5,-3), xy=( PVmax_list[i]        , Tmax_list[i] ) )

    autoRange( ax[5], 'x', 0.1, 0.2 )

    fnames = [
        "chist-Tmax"     , # 0
        "chist-Tmax-zoom", # 1
        "k-Tmax"         , # 2
        "k-Tmax-zoom"    , # 3
        "PV-Tmax"        , # 4
        "PV-Tmax-zoom"   , # 5
        "PVmax-Tmax-zoom"  # 6
    ]
    for i in range(len(fnames)):
        finishPlot( fig[i], ax[i], fnames[i] )

def plotMischungsbruch():
    ############ Z_Bilger, Z_ULF and Z_calc comparisons for 5 chi_st ###########
    fig, ax = [], []
    for i in range(6):
        fig.append( plt.figure( figsize=figTiny ) )
    sr_ZULF = r"Mischungsbruch $Z_\mathrm{ULF}$"
    # Z_Bilfer-Z
    ax.append( fig[0].add_subplot( 111,
        xlabel = sr_ZULF,
        ylabel = r"$Z_\mathrm{Bilger}$",
        xlim   = [0,1],
        ylim   = [0,1]
    ) )
    # Z_Bilfer-Z_ULF-diff
    ax.append( fig[1].add_subplot( 111,
        xlabel = sr_ZULF,
        ylabel = r"$Z_\mathrm{Bilger}-Z_\mathrm{ULF}$",
        xlim   = [0,1]
    ) )
    # Z_Bilger-Z_calc
    ax.append( fig[2].add_subplot( 111,
        xlabel = r"Mischungsbruch $Z$",
        ylabel = r"$Z_\mathrm{Bilger}$",
        xlim   = [0,1],
        ylim   = [0,1]
    ) )
    # Z_Bilger-Z_calc-diff
    ax.append( fig[3].add_subplot( 111,
        xlabel = r"Mischungsbruch $Z$",
        ylabel = r"$Z_\mathrm{Bilger}-Z$",
        xlim   = [0,1]
    ) )
    # Z_Bilger-Z_calc-relerr
    ax.append( fig[4].add_subplot( 111,
        xlabel = r"Mischungsbruch $Z$",
        ylabel = r"$\left|\frac{Z_\mathrm{Bilger}-Z}{Z}\right|$",
        xscale = 'log',
        yscale = 'log'
    ) )
    # Z_Bilger-Z_ULF-relerr
    ax.append( fig[5].add_subplot( 111,
        xlabel = sr_ZULF,
        ylabel = r"$\left|\frac{Z_\mathrm{Bilger}-Z_\mathrm{ULF}}{Z_\mathrm{ULF}}\right|$",
        xscale = 'log',
        yscale = 'log'
    ) )

    iC = -1
    for i in iSelected:
        iC += 1
        data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )
        Z       = calcZ( data[:,hdict["CH4"]], data[:,hdict["O2"]] )
        ZBilger = calcZBilger( data, hdict )
        ZUlf    = data[:,hdict["Z"]]

        #print "Z       : ",Z
        #print "ZUlf    : ",ZUlf
        #print "ZBilger : ",ZBilger
        print "min(ZBilger)  = ", np.min(ZBilger), ", max(ZBilger) = ",np.max(ZBilger)
        print "min(Z_ULF)    = ", np.min(ZUlf)   , ", max(Z_ULF)   = ",np.max(ZUlf)
        print "ZBilger > 0.8 : ", np.sum(ZBilger > 0.8)

        ax[0].plot( ZUlf, ZBilger       , '.', color=colors[iC], label=chist_sr[i] )
        ax[1].plot( ZUlf, ZBilger-ZUlf  , '.', color=colors[iC], label=chist_sr[i] )
        ax[2].plot( Z   , ZBilger       , '.', color=colors[iC], label=chist_sr[i] )
        ax[3].plot( Z   , ZBilger-Z     , '.', color=colors[iC], label=chist_sr[i] )
        x,y = relErr(ZBilger,Z)
        ax[4].plot( x,y, '.', color=colors[iC], label=chist_sr[i] )
        x,y = relErr(ZBilger,ZUlf)
        ax[5].plot( x,y, '.', color=colors[iC], label=chist_sr[i] )
        #print "ZUlf (almost zero sorted out) : ",x
        #print "ZBilger (almost zero sorted out) : ",y

    ax[0].plot( [0,1], [0,1], 'k--', label=u"Identität" )
    ax[2].plot( [0,1], [0,1], 'k--', label=u"Identität" )

    autoLabel( ax[4], 'x', nbins=6 )
    autoLabel( ax[5], 'x', nbins=6 )
    autoRangeXY( ax[4] )
    autoRangeXY( ax[5] )
    # fuck doesn't fucking work
    #ax[4].locator_params( axis='x', nbins=5 )  # reduce number of ticks because of overalpping
    #ax[5].locator_params( axis='x', nbins=5 )  # reduce number of ticks because of overalpping
    filenames=[
        "ZULF-ZBilger-for-5-chist"        , # 0
        "ZULF-ZBilger-diff-for-5-chist"   , # 1
        "Zcalc-ZBilger-for-5-chist"       , # 2
        "Zcalc-ZBilger-diff-for-5-chist"  , # 3
        "Zcalc-ZBilger-relerr-for-5-chist", # 4
        "ZULF-ZBilger-relerr-for-5-chist"   # 5
    ]
    for i in range(len(fig)):
        finishPlot( fig[i], ax[i], filenames[i] )

def plotErfcProfiles():
    ############## chi vs. erfc profile (using interpolated chist) #############
    from scipy.special import erfcinv
    fig, ax = [], []
    for i in range(4):
        fig.append( plt.figure( figsize=figTiny ) )
        ax.append( fig[-1].add_subplot( 111,
            xlabel = r"Mischungsbruch $Z_\mathrm{ULF}$",
            ylabel = r"Skalare Dissipationsrate $\chi / s^{-1}$",
            xlim   = [0,1]
        ) )
    ax[0].plot( [0,0], 'k--', label=r'analytisch' )
    ax[1].plot( [0,0], 'k--', label=r'analytisch$\cdot 1.2$' )
    ax[2].plot( [0,0], 'k--', label=r'analytisch' )
    ax[3].plot( [0,0], 'r-' , label=r'analytisch$\cdot 1.2$' )
    ax[3].plot( [0,0], 'b-' , label=r'analytisch' )
    iComp = abs( chist_intp_list - 1.5 ).argmin()

    counter = -1
    for i in iSelected:
        counter += 1
        data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )
        ZUlf = data[:,hdict["Z"]]

        # simulated / experimental chi
        for j in range(len(ax)):
            if j==3:
                continue
            ax[j].plot( ZUlf, calcChi(data,hdict), '.', color=colors[counter],
                        label=chist_sr[i] )
        # chi-erfc-profile
        chianal = chist_list[i]*np.exp(2.*( erfcinv(2.*Zstanal)**2 -
                                            erfcinv(2.*ZUlf   )**2 ))
        ax[0].plot( ZUlf, chianal, '--', color=colors[counter] )
        # chi-1.2*erfc-profile
        ax[1].plot( ZUlf, 1.2*chianal, '--', color=colors[counter] )
        # chi_interpolated-erfc-profile
        ax[2].plot( ZUlf, chianal / chist_list[i] * chist_intp_list[i], '--', color=colors[counter] )

    for i in [iComp]:
        counter += 1
        data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )
        ZUlf = data[:,hdict["Z"]]

        # chi-erfc-profile
        chianal = chist_list[i]*np.exp(2.*( erfcinv(2.*Zstanal)**2 -
                                            erfcinv(2.*ZUlf   )**2 ))

        ax[3].plot( ZUlf, 1.2*chianal, 'r-' )
        ax[3].plot( ZUlf, chianal / chist_list[i] * chist_intp_list[i], 'b-' )

        # simulated / experimental chi
        ax[3].plot( ZUlf, calcChi(data,hdict), 'k.', markersize=2.0,
                    label=chist_sr[i] )

    #ax[3].set_title( chist_sr[iComp] )
    ax[3].set_xlim([0,0.1]) # roughly 2*Zstanal as upper limit
    ax[3].set_ylim( [0,5] )
    ax[3].plot( [Zstanal,Zstanal], ax[3].get_ylim(), 'k--' , label=r"$Z_\mathrm{stoch}$" )

        #chianal = strainrate_list[i]/np.pi*np.exp(-2.*erfcinv(2.*ZUlf)**2)
    filenames = [
        "chianal"     , # 1
        "chianal1.2"  , # 2
        "chianal-intp", # 3
        "chianal-zoom"  # 4
    ]
    for i in range(len(fig)):
        finishPlot( fig[i], ax[i], filenames[i] )

v_list         , \
Tmax_list      , \
Tst_list       , \
chist_list     , \
chist_sr       , \
chist_intp_list, \
strainrate_list, \
PVst_list      , \
PVmax_list     = precalcOppdiffValues( './results' )

iSelected = [ i for i in range( len(v_list) ) if v_list[i] in vSelected ]

######################### y limits for species profiles ########################

# find a suitable and comparable range for Y_i and T (maximum  / minimum of all files to read and compare)
yminT  = 1e7
ymaxT  = 0
yminYi = 1.0
ymaxYi = 0
for i in range(len(v_list)):
    data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )
    yminT  = min( yminT , min( data[:,hdict["T"]] ) )
    ymaxT  = max( ymaxT , max( data[:,hdict["T"]] ) )
    for name in ["CO","CO2","H2","H2O","OH"]:
        yminYi = min( yminYi, min( data[:,hdict[name]] ) )
        ymaxYi = max( ymaxYi, max( data[:,hdict[name]] ) )
# HARDCODED! (ALL IS FUTILE)
yminT  = 300
ymaxT  = 2300
yminYi = 0.0
ymaxYi = 0.13

################# Plots #################

plotSCurves()
plotErfcProfiles()
plotMischungsbruch()

################################### x over ix ##################################

data, hdict = readUlfFile( "./results/v_0.5" ) # "./oppdifJet_ct_CH4Air_smooke_final.ulf"
# access to column labeld T: data[:,hdict["T"]]
Zcalc   = calcZ( data[:,hdict["CH4"]], data[:,hdict["O2"]] )
x       = data[:,hdict["X"]]
ZUlf    = data[:,hdict["Z"]]
Nx      = len(x)
ist     = abs(ZUlf-Zstanal).argmin()
xst     = x[ist]
Zst     = ZUlf[ist]

chist   = calcChi( data, hdict )[ist] # i-stoch calculated above with Z_st = 0.0545
Tmax    = np.max( data[:,hdict["lambda"]] )
print "i_st = ",ist,", x_st = ",xst,", Z_st = ",Zst,", chi_st = ",chist,", T_max = ",Tmax

fig = plt.figure( figsize=figSmall )
ax = fig.add_subplot( 111,
    xlabel = u"Zähler",
    ylabel = "Ort $x / m$",
    xlim   = [0,len(x)],
    ylim   = [x.min(),x.max()]
)
ax.plot( x, 'b.' )
finishPlot( fig, ax, "ix" )


################################### Z over x ###################################

#ZBilger = calcZBilger(data, hdict, True)
#   M_O2 =  31.9988 , M_CH4 =  16.04246 , M_O =  15.9994 , M_C =  12.0107
#   beta_B  =  0.124669159219 , beta_Ox  =  -0.00718776954136
#   beta[0] =  0.124669159219 , beta[-1] =  -0.00718776954136
ZBilger = calcZBilger( data, hdict )
fig = plt.figure( figsize=figSmall )
ax  = fig.add_subplot( 111 )
ax.set_xlabel( r"Ort $x / m$" )
ax.set_ylabel( r"Mischungsbruch $Z$" )
ax.set_xlim( [x.min(),x.max()] )
ax.set_ylim( [-0.2,1.05] )
ax.plot( x, Zcalc            , 'bx', label=r"$Z(x)$" )
ax.plot( x, ZBilger          , 'm.', label=r"$Z_\mathrm{Bilger}(x)$" )
ax.plot( x, ZUlf             , 'k-', label=r"$Z_\mathrm{ULF}(x)$" )
ax.plot( x, 50*(ZUlf-Zcalc)  , 'g+', label=r"$50 \cdot \left( Z_\mathrm{ULF}(x) - Z(x) \right)$" )
ax.plot( x, 50*(ZUlf-ZBilger), 'y+', label=r"$50 \cdot \left( Z_\mathrm{ULF}(x) - Z_\mathrm{Bilger}(x) \right)$" )
ax.plot( ax.get_xlim(), [Zstanal,Zstanal], 'k--', label=r"$Z_\mathrm{stoch}$" )
ax.plot( [xst,xst], ax.get_ylim(), 'k:' , label=r"$x_\mathrm{stoch}$" )
#ax.plot( xst, Zstanal        , 'ro', label=r"$Z_\mathrm{stoch}$" )
#ax.plot( x  , Zcalc          , 'b.', label=r"$Z(x) = \frac{ \nu Y_{\mathrm{CH}_4}(x) - Y_{\mathrm{O}_2}(x) + Y_{\mathrm{O}_2,\mathrm{Ox}} }{ \nu Y_{\mathrm{CH}_4,\mathrm{B}} + Y_{\mathrm{O}_2,\mathrm{Ox}} }$" )
#ax.plot( relErr( ZUlf, Zcalc   ), 'g+', label=r"$\frac{ \left( Z_\mathrm{ULF}(x) - Z(x) \right) }{ Z(x) }$" )
#ax.plot( relErr( ZUlf, ZBilger ), 'y+', label=r"$\frac{\left( Z_\mathrm{ULF}(x) - Z_\mathrm{Bilger}(x) \right) }{ Z_\mathrm{Bilger}(x) }$" )
finishPlot( fig, ax, "Mischungsbruchvergleich" )
non0,err = relErr( ZUlf, Zcalc )
print "[Mischungsbruchvergleich.pdf] max. rel. error (Z_ULF-Z)/Z =", max( abs( err ) )
non0,err = relErr( ZUlf, ZBilger )
print "[Mischungsbruchvergleich.pdf] max. rel. error (Z_ULF-Z_Bilger)/Z_Bilger =", max( abs( err ) )
plt.close( fig )


############################### Z_Bilger over Z ################################

fig = plt.figure( figsize=figSmall )
ax = fig.add_subplot( 111,
    xlabel = r"Mischungsbruch $Z$",
    ylabel = r"Elementmischungsbruch $Z_\mathrm{Bilger}$",
    xlim   = [0,1],
    ylim   = [0,1]
)
ax.plot  ( ZUlf, ZBilger, 'b.' )
finishPlot( fig, ax, "Z-ZBilger" )


########################### species profiles over x ############################

for i in range(len(v_list)):
    data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )

    fig = plt.figure( figsize=figTiny )
    plt.title ( r"$v="+str(v_list[i]) + r"\,\mathrm{m}/\mathrm{s}$, " + chist_sr[i], fontsize=12 )

    x = data[:,hdict["X"]]
    axL = fig.add_subplot(111)
    # putting these options inside add_subplot will weirdly result in weird double tick names ... bug ???
    axL.set_xlabel( r"Ort x / m" )
    axL.set_ylabel( u"Massenbrüche "+r"$Y_i$" )
    axL.set_xlim( np.min(x),np.max(x) )
    axL.set_ylim( yminYi, ymaxYi + 0.1*(ymaxYi-yminYi) )
    for name in ["CO","CO2","H2","H2O","OH"]:
        axL.plot( x, data[:,hdict[name]], '-', label=r"$Y_\mathrm{"+name+r"}$" )
    axL.plot( [xst,xst], axL.get_ylim(), 'k:', label=r"$x_\mathrm{stoch}$" )
    axR = axL.twinx()
    axR.plot( x, data[:,hdict["T"]], 'k--', label="T" )
    axL.plot( 0,0, 'k--', label="T" ) # only plot in order to include it in legend!
    axR.set_ylabel( r"Temperatur $T / ^\circ\mathrm{K}$" )
    axR.set_xlim( np.min(x),np.max(x) )
    axR.set_ylim( yminT, ymaxT + 0.1*(ymaxT-yminT) )
    finishPlot( fig, axL, "Yi-over-x-v-"+str(v_list[i]) )


########################### species profiles over Z ############################

for i in range(len(v_list)): # only use 5 lowest velocities
    data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )

    fig = plt.figure( figsize=figTiny )
    plt.title ( r"$v="+str(v_list[i]) + r"\,\mathrm{m}/\mathrm{s}$, " + chist_sr[i], fontsize=12 )

    ZUlf = data[:,hdict["Z"]]
    axL = plt.subplot(111)
    axL.set_xlabel( r"Mischungsbruch $Z_\mathrm{ULF}$" )
    axL.set_ylabel( u"Massenbrüche "+r"$Y_i$" )
    axL.set_xlim( 0,1 )
    axL.set_ylim( yminYi, ymaxYi + 0.1*(ymaxYi-yminYi) )
    for name in ["CO","CO2","H2","H2O","OH"]:
        axL.plot( ZUlf, data[:,hdict[name]], '-', label=r"$Y_\mathrm{"+name+r"}$" )
    axL.plot( [Zstanal,Zstanal], axL.get_ylim(), 'k:', label=r"$Z_\mathrm{stoch}$" )
    axR = axL.twinx()
    axR.plot( ZUlf, data[:,hdict["T"]], 'k--', label="T" )
    axL.plot(0,0,'k--',label="T") # ghost plot for legend
    axR.set_ylabel( r"Temperatur $T / ^\circ\mathrm{K}$" )
    axR.set_xlim( 0,1 )
    axR.set_ylim( yminT, ymaxT + 0.1*(ymaxT-yminT) )
    finishPlot( fig, axL, "Yi-over-Z-v-"+str(v_list[i]) )


################################# Y_i over PV ##################################

for i in range(len(v_list)):
    data, hdict = readUlfFile( 'results/v_'+str(v_list[i]) )
    PV  = data[:,hdict["CO"]] + data[:,hdict["CO2"]]
    fig = plt.figure( figsize=figTiny )
    ax  = fig.add_subplot( 111,
        title  = r"$v="+str(v_list[i]) + r"\,\mathrm{m}/\mathrm{s}$, " + chist_sr[i],
        xlabel = r"Fortschrittsvariable PV / %",
        ylabel = u"Massenbrüche "+r"$Y_i$ / %",
        xlim   = [ np.min(PV)*100, 1.1*np.max(PV)*100 ],
        ylim   = [ yminYi    *100, 1.1*ymaxYi    *100 ]
    )
    ax.plot( [ PVst_list[i]*100, PVst_list[i]*100 ], ax.get_ylim(), 'k:', label=r"$\mathrm{PV}_\mathrm{stoch}$" )
    ax.plot( [ PVmax_list[i]*100, PVmax_list[i]*100 ], ax.get_ylim(), 'k--', label=r"$\mathrm{PV}_{T_\mathrm{max}}$" )
    for name in ["CO","CO2","H2","H2O","OH"]:
        Yi = data[:,hdict[name]]
        ax.plot( PV[ :len(Yi)/2 ]*100, Yi[ :len(Yi)/2 ]*100, '.', label=r"$Y_\mathrm{"+name+r"}$" )
        #ax.plot( PV*100, Yi*100, '.', label=r"$Y_\mathrm{"+name+r"}$" )
    finishPlot( fig, ax, "Yi-over-PV-v-"+str(v_list[i])+"-half" )

    ax.lines = []
    ax.set_color_cycle(None)
    ax.plot( [ PVst_list[i]*100, PVst_list[i]*100 ], ax.get_ylim(), 'k:', label=r"$\mathrm{PV}_\mathrm{stoch}$" )
    ax.plot( [ PVmax_list[i]*100, PVmax_list[i]*100 ], ax.get_ylim(), 'k--', label=r"$\mathrm{PV}_{T_\mathrm{max}}$" )
    for name in ["CO","CO2","H2","H2O","OH"]:
        Yi = data[:,hdict[name]]
        ax.plot( PV*100, Yi*100, '.', label=r"$Y_\mathrm{"+name+r"}$" )
    finishPlot( fig, ax, "Yi-over-PV-v-"+str(v_list[i]) )

#plt.show()

