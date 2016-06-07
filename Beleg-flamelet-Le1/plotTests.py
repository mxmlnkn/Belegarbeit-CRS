#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os
from plotCommon import *

opp_vSelected = vSelected
sr_chist  = [ '0.01','0.1','0.25','0.4918','0.8800','1.2686','1.6383','1.9816','3.7859' ]

def plotSCurves():
    # needs: chist_list Tmax_list strainrate_list PVst_list chist_list PVmax_list v_list iSelected finishPlot
    fig, ax = [], []
    for i in range(7):
        fig.append( plt.figure( figsize=figSmall ) )
    l_chist = r"$\chi_\mathrm{st} / s^{-1}$"
    l_k     = r"Inverse Streckungsrate $k^{-1} / \mathrm{s}$"
    l_pvst  = r"$\mathrm{PV}_\mathrm{st} := Y_\mathrm{CO}( Z_\mathrm{st} ) + Y_{\mathrm{CO}_2}( Z_\mathrm{st} ) $"
    xlabels = [ l_chist, l_chist, l_pvst, l_pvst, l_chist, l_pvst, l_chist ]
    l_Tmax  = r"$T_\mathrm{max} / ^\circ\mathrm{K}$"
    l_Tst   = r"$T_\mathrm{st} / ^\circ\mathrm{K}$"
    ylabels = [ l_Tmax, l_Tmax, l_Tmax, l_Tmax, l_Tmax, l_Tmax, l_Tst ]

    for i in range(len(fig)):
        ax.append( fig[i].add_subplot( 111, xlabel = xlabels[i],
                   ylabel = ylabels[i], xmargin = 0.1, ymargin = 0.1 ) )

    # Tmax over chi_st (for Le=1 and Le != 1)
    ax[4].plot( fLe1chist, fLe1Tmax, 'bo-', label=fLe1label )
    ax[4].plot( oLe1chist, oLe1Tmax, 'ro-', label=oLe1label )
    ax[4].plot( fLeVchist, fLeVTmax, 'co-', label=fLeVlabel )
    ax[4].plot( oLeVchist, oLeVTmax, 'o-' , label=oLeVlabel, color='orange' )
    autoRange( ax[4], 'x', 0, 0.2 )
    chistLim = ax[4].get_xlim()
    chistLim = ( 0, chistLim[1] )
    ax[4].set_xlim( chistLim )

    # Tmax over PV (zoom) | PV-Tmax-comparison.pdf
    ax[5].plot( fLe1PVst, fLe1Tmax, 'bo-', label=fLe1label )
    ax[5].plot( oLe1PVst, oLe1Tmax, 'ro-', label=oLe1label )
    ax[5].plot( fLeVPVst, fLeVTmax, 'co-', label=fLeVlabel )
    ax[5].plot( oLeVPVst, oLeVTmax, 'o-' , label=oLeVlabel, color='orange' )
    #ax[5].tick_params(axis='both', which='major', labelsize=10)
    #ax[5].tick_params(axis='both', which='minor', labelsize=8 )
    autoRange( ax[5], 'y', 0.15, 0.86 )
    autoRange( ax[5], 'x', 0.1 , 0.1 )
    PVLim = ax[4].get_xlim()

    # S-Curves
    burning = Tmax_list > 500
    ax[0].plot( chist_list         , Tmax_list         , 'bo'  )
    ax[0].plot( chist_list[burning], Tmax_list[burning], 'bo-', label=fLe1label )
    #ax[0].plot( opp_chist_list, opp_Tmax_list, 'ro'  )
    #ax[0].plot( opp_chist_list[ opp_Tmax_list > 500 ],
    #            opp_Tmax_list [ opp_Tmax_list > 500 ], 'ro-', label=oLe1label )
    ax[0].set_xlim( chistLim )

    ax[2].plot( PVst_list[burning], Tmax_list[burning], 'bo-' )
    autoRange( ax[2], 'x', 0.1, 0.15 )

    # Tmax over chi_st (zoom) (for Le=1)
    ax[1].plot( fLe1chist, fLe1Tmax, 'bo-', label=fLe1label )
    ax[1].plot( oLe1chist, oLe1Tmax, 'ro-', label=oLe1label )
    ax[1].set_xlim( chistLim )

    # Tmax over PV (zoom) (for Le=1)
    ax[3].plot( fLe1PVst, fLe1Tmax, 'bo-', label=fLe1label )
    ax[3].plot( oLe1PVst, oLe1Tmax, 'ro-', label=oLe1label )
    #ax[3].tick_params(axis='both', which='major', labelsize=10)
    #ax[3].tick_params(axis='both', which='minor', labelsize=8 )
    autoRange( ax[3], 'x', 0.1, 0.2 )

    # Comparison: T(chi_st) over chi_st
    ax[6].plot( fLe1chist, fLe1Tst, 'bo-', label=fLe1label )
    ax[6].plot( oLe1chist, oLe1Tst, 'ro-', label=oLe1label )
    ax[6].plot( fLeVchist, fLeVTst, 'co-', label=fLeVlabel )
    ax[6].plot( oLeVchist, oLeVTst, 'o-' , label=oLeVlabel, color='orange' )
    ax[6].set_xlim( chistLim )

    # make annotations on labels where needed
    for i in range( len(chist_list) ):
        # S curve Tmax over chi_st
        #if chist_list[i] >= 0.4:
        #    ax[0].annotate( labels[i], xycoords='data', textcoords='offset points',
        #        size=9, xytext=(-8,5), xy=( chist_list[i], Tmax_list[i] ) )
        ax[2].annotate( format(chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(5,-3), xy=( PVst_list[i] , Tmax_list[i] ) )

    for i in iSelected:
        #ax[1].annotate( labels[i], xycoords='data', textcoords='offset points',
        #    size=9, xytext=(5,0), xy=( chist_list[i], Tmax_list[i] ) )
        ax[3].annotate( format(chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(5,-3), xy=( PVst_list[i] , Tmax_list[i] ) )
        ax[5].annotate( format(chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(-26,-3), xy=( PVst_list[i] , Tmax_list[i] ) )

    for i in opp_iSelected:
        ax[3].annotate( format(opp_chist_list[i],'.2f'), xycoords='data', textcoords='offset points',size=9, xytext=(5,-3), xy=( opp_PVst_list[i] , opp_Tmax_list[i] ) )
        ax[5].annotate( format(opp_chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(5,-3), xy=( opp_PVst_list[i] , opp_Tmax_list[i] ) )

    for i in flameletLeVar_iSelected:
        ax[5].annotate( format(flameletLeVar_chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(5,-3), xy=( flameletLeVar_PVst_list[i] , flameletLeVar_Tmax_list[i] ) )

    for i in oppLeVar_iSelected:
        ax[5].annotate( format(oppLeVar_chist_list[i],'.2f'), xycoords='data', textcoords='offset points', size=9, xytext=(-26,-3), xy=( oppLeVar_PVst_list[i] , oppLeVar_Tmax_list[i] ) )

    fnames = [
        "chist-Tmax"           ,  # 0
        "chist-Tmax-zoom"      ,  # 1
        "PV-Tmax"              ,  # 2
        "PV-Tmax-zoom"         ,  # 3
        "chist-Tmax-comparison",  # 4
        "PV-Tmax-comparison"   ,  # 5
        "chist-Tst-comparison"    # 6
    ]
    for i in range(len(fnames)):
        finishPlot( fig[i], ax[i], fnames[i] ) # last is 'outside'

def plotErfcProfiles():
    # chi vs. erfc profile (using interpolated chist)
    # needs: chist_list Zstanal ZUlf chist_intp_list
    from scipy.special import erfcinv
    fig, ax = [], []
    for i in range(2):
        fig.append( plt.figure( figsize=figSmall ) )
        ax.append( fig[-1].add_subplot( 111,
            xlabel = r"Mischungsbruch $Z_\mathrm{ULF}$",
            ylabel = r"Skalare Dissipationsrate $\chi / s^{-1}$",
            xlim   = [0,1]
        ) )
    ax[0].plot( [0,0], 'k--', label=r'analytisch' )
    ax[1].plot( [0,0], 'k--', label=r'analytisch' )
    counter = -1
    for i in iSelected:
        counter += 1
        data, hdict = readUlfFile( 'results/chist_'+str(chist_list[i]) )
        ZUlf = data[:,hdict["Z"]]

        # simulated / experimental chi
        for j in range(len(ax)):
            ax[j].plot( ZUlf[::4], calcChi(data,hdict)[::4], '.',
                        color=colors[counter], label=chist_sr[i] )
        # chi-erfc-profile
        chianal = chist_list[i]*np.exp(2.*( erfcinv(2.*Zstanal)**2 -
                                            erfcinv(2.*ZUlf   )**2 ))
        ax[0].plot( ZUlf, chianal, '--', color=colors[counter] )
        # chi_interpolated-erfc-profile
        ax[1].plot( ZUlf, chianal / chist_list[i] * chist_calc_list[i], '--', color=colors[counter] ) #chist_intp_list[i]
        #chianal = strainrate_list[i]/np.pi*np.exp(-2.*erfcinv(2.*ZUlf)**2)
    filenames = [ "chianal", "chianal-calc" ]
    for i in range(len(fig)):
        finishPlot( fig[i], ax[i], filenames[i] )

opp_v_list         , \
opp_Tmax_list      , \
opp_Tst_list       , \
opp_chist_list     , \
opp_chist_sr       , \
opp_chist_intp_list, \
opp_strainrate_list, \
opp_PVst_list      , \
opp_PVmax_list     = precalcOppdiffValues( '../Beleg-oppdiff-Le1/results' )
opp_iSelected = [ i for i in range( len(opp_v_list) ) if opp_v_list[i] in opp_vSelected ]

oppLeVar_v_list         , \
oppLeVar_Tmax_list      , \
oppLeVar_Tst_list       , \
oppLeVar_chist_list     , \
oppLeVar_chist_sr       , \
oppLeVar_chist_intp_list, \
oppLeVar_strainrate_list, \
oppLeVar_PVst_list      , \
oppLeVar_PVmax_list     = precalcOppdiffValues( '../Beleg-oppdiff-Levar/results' )
oppLeVar_iSelected = [ i for i in range( len(oppLeVar_v_list) ) if oppLeVar_v_list[i] in opp_vSelected ]

fname_list     , \
chist_list     , \
chist_sr       , \
chist_calc_list, \
chist_intp_list, \
Tmax_list      , \
Tst_list       , \
PVst_list      , \
PVmax_list = precalcFlameletValues( './results' )
iSelected = np.arange(len(chist_list))[ np.logical_and( chist_list >= 0.1, chist_list <= 2.0 ) ]

flameletLeVar_fname_list     , \
flameletLeVar_chist_list     , \
flameletLeVar_chist_sr       , \
flameletLeVar_chist_calc_list, \
flameletLeVar_chist_intp_list, \
flameletLeVar_Tmax_list      , \
flameletLeVar_Tst_list       , \
flameletLeVar_PVst_list      , \
flameletLeVar_PVmax_list = precalcFlameletValues( '../Beleg-flamelet-Levar/results' )
flameletLeVar_iSelected = np.arange(len( flameletLeVar_chist_list ))[ np.logical_and( flameletLeVar_chist_list >= 0.1, flameletLeVar_chist_list <= 2.0 ) ]


fLe1chist =               chist_list[               iSelected ]
oLe1chist =           opp_chist_list[           opp_iSelected ]
fLeVchist = flameletLeVar_chist_list[ flameletLeVar_iSelected ]
oLeVchist =      oppLeVar_chist_list[      oppLeVar_iSelected ]

fLe1Tmax  =                Tmax_list[               iSelected ]
oLe1Tmax  =            opp_Tmax_list[           opp_iSelected ]
fLeVTmax  =  flameletLeVar_Tmax_list[ flameletLeVar_iSelected ]
oLeVTmax  =       oppLeVar_Tmax_list[      oppLeVar_iSelected ]

fLe1Tst   =                Tmax_list[               iSelected ]
oLe1Tst   =            opp_Tmax_list[           opp_iSelected ]
fLeVTst   =  flameletLeVar_Tmax_list[ flameletLeVar_iSelected ]
oLeVTst   =       oppLeVar_Tmax_list[      oppLeVar_iSelected ]

fLe1PVst  =                PVst_list[               iSelected ]
oLe1PVst  =            opp_PVst_list[           opp_iSelected ]
fLeVPVst  =  flameletLeVar_PVst_list[ flameletLeVar_iSelected ]
oLeVPVst  =       oppLeVar_PVst_list[      oppLeVar_iSelected ]

fLe1label =r"flamelets $\mathrm{Le}=1$"
oLe1label =r"oppdiffJet $\mathrm{Le}=1$"
fLeVlabel =r"flamelets $\mathrm{Le}\neq 1$"
oLeVlabel =r"oppdiffJet $\mathrm{Le}\neq 1$"


######################### y limits for species profiles ########################

# find a suitable and comparable range for Y_i and T (maximum  / minimum of all files to read and compare)
yminT  = 1e7
ymaxT  = 0
yminYi = 1.0
ymaxYi = 0
for i in iSelected:
    data, hdict = readUlfFile( './results/chist_'+str(chist_list[i]) )
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

##################### Fast Plot to comapre T_st with T_max #####################

fig = plt.figure( figsize=([6,4]) )
ax  = fig.add_subplot( 111,
    xlabel = r"$T_\mathrm{st}$",
    ylabel = r"$T_\mathrm{max} - T_\mathrm{st}$"
)
ax.plot( fLe1Tst, fLe1Tmax - fLe1Tst, label=fLe1label )
ax.plot( oLe1Tst, oLe1Tmax - oLe1Tst, label=oLe1label )
ax.plot( fLeVTst, fLeVTmax - fLeVTst, label=fLeVlabel )
ax.plot( oLeVTst, oLeVTmax - oLeVTst, label=oLeVlabel )
finishPlot( fig, ax, "Tmax-Tst" )


####################### Z_Bilger vs. Zcalc for 5 chi_st ########################
# nur aus Interesse, keine Pflicht

fig = plt.figure( figsize=([6,4]) )
ax  = fig.add_subplot( 111,
    xlabel = r"Mischungsbruch $Z_\mathrm{calc}$",
    ylabel = r"$Z_\mathrm{Bilger}$",
    xlim   = [0,1],
    ylim   = [0,1]
)
for i in range(len(chist_list)):
    data, hdict = readUlfFile( fname_list[i] )
    Z       = calcZ( data[:,hdict["CH4"]], data[:,hdict["O2"]] )
    ZBilger = calcZBilger( data, hdict )
    print "ZBilger over Z_calc: min(ZBilger) = ",np.min(ZBilger),", max(ZBilger) = ",np.max(ZBilger)
    print "ZBilger over Z_calc: min(Zcalc) = ",np.min(Z),", max(Zcalc) = ",np.max(Z)
    print "ZBilger > 0.8: ",np.sum(ZBilger > 0.8)

    ax.plot( Z, ZBilger, '.', label=chist_sr[i] )
    ax.plot( Z, 10*(ZBilger-Z), '.', label=chist_sr[i] )
finishPlot( fig, ax, "Zcalc-ZBilger-5-chist" )


############################### species profiles ###############################

for i in range(len(chist_list)):
    data, hdict = readUlfFile( 'results/chist_'+str(chist_list[i]) )
    fig = plt.figure( figsize=figTiny )
    ZUlf = data[:,hdict["Z"]]
    axL = fig.add_subplot( 111,
        title  = chist_sr[i],
        xlabel = r"Mischungsbruch $Z_\mathrm{ULF}$",
        ylabel = r"Massenbrueche $Y_i$"
    )
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
    finishPlot( fig, axL, "Yi-over-Z-chist-"+sr_chist[i] )


######################### species profiles (comparison) ########################

species=["CO","CO2","H2O","H2","OH","CH4","O2","T"]
titles=[ r"$\mathrm{CO}$",
         r"$\mathrm{CO}_2$",
         r"$\mathrm{H}_2\mathrm{O}$",
         r"$\mathrm{H}_2$",
         r"$\mathrm{OH}$",
         r"$\mathrm{CH}_4$",
         r"$\mathrm{O}_2$",
         "Temperature" ]
fig,ax = [],[]
for i in range(len(species)):
    fig.append( plt.figure( figsize=figTiny ) )
    ax.append( fig[-1].add_subplot( 111,
        title  = titles[i],
        xlabel = r"Mischungsbruch $Z_\mathrm{ULF}$",
        ylabel = r"Massenbrueche $Y_\mathrm{"+species[i]+r"}$",
        xlim   = [ 0, 1 ]
    ) )
ax[-1].set_ylabel( r"Temperatur $T / ^\circ\mathrm{K}$" )
ax[-1].set_ylim( yminT, ymaxT + 0.1*(ymaxT-yminT) )
ax[-1].set_title( r"" )

configurations=[ 'results/chist_1.9816',
                 '../Beleg-flamelet-Levar/results/chist_1.9816',
                 '../Beleg-oppdiff-Le1/results/v_0.25',
                 '../Beleg-oppdiff-Levar/results/v_0.25' ]
labels=[ r'Flamelets $\mathrm{Le}=1$',
         r'Flamelets $\mathrm{Le}\neq1$',
         r'OppDiffJet $\mathrm{Le}=1$',
         r'OppDiffJet $\mathrm{Le}\neq1$' ]
for ifname in range(len(configurations)):
    data, hdict = readUlfFile( configurations[ifname] )
    for i in range(len(species)-1):
        ax[i].plot( data[:,hdict["Z"]], data[:,hdict[species[i]]], '-', label=labels[ifname] )
    ax[-1].plot( data[:,hdict["Z"]], data[:,hdict["T"]], '-', label=labels[ifname] )
for i in range(len(species)):
    ax[i].plot( [Zstanal,Zstanal], ax[i].get_ylim(), 'k:', label=r"$Z_\mathrm{stoch}$" )

for i in range(len(species)):
    finishPlot( fig[i], ax[i], "Yi-over-Z-chist-1.98-comp-"+species[i] )


################################# Y_i over PV ##################################

for i in range(len(chist_list)):
    data, hdict = readUlfFile( fname_list[i] )
    PV  = data[:,hdict["CO"]] + data[:,hdict["CO2"]]
    fig = plt.figure( figsize=figTiny )
    ax  = fig.add_subplot( 111,
        title  = chist_sr[i],
        xlabel = r"Fortschrittsvariable PV / %",
        ylabel = r"Massenbrueche $Y_i$ / %",
        xlim   = [ np.min(PV)*100, 1.1*np.max(PV)*100 ],
        ylim   = [ yminYi    *100, 1.1*ymaxYi    *100 ]
    )
    for name in ["CO","CO2","H2","H2O","OH"]:
        Yi = data[:,hdict[name]]
        #ax.plot( PV[ :len(Yi)/2 ]*100, Yi[ :len(Yi)/2 ]*100, '.', label=r"$Y_\mathrm{"+name+r"}$" )
        ax.plot( PV*100, Yi*100, '.', label=r"$Y_\mathrm{"+name+r"}$" )
    ax.plot( [ PVst_list[i]*100, PVst_list[i]*100 ], ax.get_ylim(), 'k:', label=r"$\mathrm{PV}_\mathrm{stoch}$" )
    ax.plot( [ PVmax_list[i]*100, PVmax_list[i]*100 ], ax.get_ylim(), 'k--', label=r"$\mathrm{PV}_{T_\mathrm{max}}$" )
    finishPlot( fig, ax, "Yi-over-PV-chist-"+sr_chist[i] )

#plt.show()
