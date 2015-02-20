#!/usr/env python
#
# Calculate the cross correlation of two time series
# Estimates the uncertainties using Monte Carlo Simulation

import matplotlib.pyplot as plt  # plot library
import numpy as np  # array manipulation
from astropy import table  # handle data tables
from scipy import signal, stats  # signal processing tools, statistic tools


def crosscorr(x, y, t):
    ''' Caluclate the cross correlation function and timeshifts for a
        pair of time series x,y
    '''

    # normalize input series
    x -= x.mean()
    y -= y.mean()
    x /= x.std()
    y /= y.std()

    # calculate cross-correlation function
    corr = signal.correlate(x,y)/float(len(x))

    # transform time axis in offset units
    #lags = np.arange(corr.size) - (t.size - 1)
    #tstep = (t[-1] - t[0])/float(t.size)
    #offset = lags*tstep
    tstep = t[1]-t[0]
    inc = tstep/3.0
    offset = np.arange(-max(t), max(t)+inc, tstep)

    # time shift is found for the maximum of the correlation function
    shift = offset[np.argmax(corr)]

    # new time axis to plot shifted time series
    newt = t + shift

    # correct time intervals if shift bigger than half the interval
    if min(newt) > (max(t)/2):
         newt = newt - max(t)
         shift = shift - max(t)
    elif max(newt) < (min(t)/2):
         newt = newt + min(t)
         shift = shift + min(t)

    return corr, offset, newt, shift


if __name__ == "__main__":
    '''
    MAIN RROGRAM

    Correlation of 2 time series in an ascii or fits file

    '''

#   Welcome message and general instructions ==================================
    print '''

    Time delay calculation between 2 time series

    NOTE FOR USING 2 FILES:
    The first input file must contains the time column and the first series

        For ASCII FILES is assumed that the first line contains the columns
        headers, and the columns are separated by blank spaces.
        Blank lines and lines starting with the carachter '#' are
        ignored.

        For FITS FILES are assumed that the first extension contains the data
        of the columns of interest.

        The time series are assumed to have equally spaced data (binned data).

        ERROR columns are assumed to be the 1-sigma deviation of the data.

        If indentifying columns by position, start counting from ZERO
    '''
#==============================================================================

#=== User Interface to Read Files =============================================
    fmt = str(raw_input("Input files format (ascii/fits): "))

    inpt1 = str(raw_input("File 1: "))
    inpt2 = str(raw_input("File 2: "))

#    tb1 = table.Table.read(inpt1, format=fmt, hdu=1)
#    tb2 = table.Table.read(inpt2, format=fmt, hdu=1)
    
    if fmt=='fits':
        tb1 = table.Table.read(inpt1, format=fmt, hdu=1)
        tb2 = table.Table.read(inpt2, format=fmt, hdu=1)
    if fmt=='ascii':
        tb1 = table.Table.read(inpt1, format=fmt)
        tb2 = table.Table.read(inpt2, format=fmt)

    timecol = raw_input("Time column name or position (starting from zero): ")
    xcol = raw_input("series 1 column name or position (starting from zero): ")
    xerr = raw_input("series 1 error column name or position (starting from zero): ")
    ycol = raw_input("series 2 column name or position (starting from zero): ")
    yerr = raw_input("series 2 error column name or position (starting from zero): ")
#==============================================================================

#==== Interpret user input  ==================================================
#===  if columns is an integer use as column position
#===  if not assume its a string (column name)
#===  other cases will raise an ERROR and stop execution
    try:
        timecol = int(timecol)
    except ValueError:
        timecol = str(timecol)

    try:
        xcol = int(xcol)
    except ValueError:
        xcol = str(xcol)

    try:
        xerr = int(xerr)
    except ValueError:
        xerr = str(xerr)

    try:
        ycol = int(ycol)
    except ValueError:
        ycol = str(ycol)

    try:
        yerr = int(yerr)
    except ValueError:
        yerr = str(yerr)
#=============================================================================

#   Store columns in variables x(t), y(t), t
    t = tb1.field(timecol)
    x = tb1.field(xcol)
    xe = tb1.field(xerr)
    y = tb2.field(ycol)
    ye = tb2.field(yerr)

#-----------------------------------------------------------------------------
#   Change NaN and Negative values to zero (as xronos crosscorr task
#   see: https://heasarc.gsfc.nasa.gov/docs/xanadu/xronos/help/crosscor.html

    excludedx = []
    excludedy = []
    numexx = 0
    numexy = 0
    print ''
    print 'Nan and negative values = 0'
    print ''
    for i in xrange(len(x)):
        if (x[i] >=0):
            pass
        else:
            excludedx.append((i,x[i]))
            numexx += 1
            x[i] = 0
            xe[i] = 0

        if (y[i] >=0):
            pass
        else:
            excludedy.append((i,y[i]))
            numexy += 1
            y[i]=0
            ye[i] = 0
#------------------------------------------------------------------------------

##  start time from Zero
    tstart = min(t)
    tend = max(t)
    t -= min(t)

### === MONTE CARLO SIMULATIONS ===============================================
    nsimulations = int(raw_input("How many simulations?: "))

#   for each point in x and y, generates 'nsimulation' new points
#   new value = (uncertainty*random.values + mean)
#   the randon points follow a normal distribution with 1-sigma equals to the
#   1-sigma error bars from the original series

    aux1 = []
    aux2 = []
    for i, meanx in enumerate(x):
        newx = xe[i]*np.random.randn(nsimulations) + meanx
        aux1.append(newx)
    for j, meany in enumerate(y):
        newy = ye[j]*np.random.randn(nsimulations) + meany
        aux2.append(newy)

#   rearange itens, newxses will contain one time series in each element
    newxses = []
    newyses = []
    for n in xrange(nsimulations):
        newxses.append(np.array([aux1[m][n] for m in xrange(len(aux1))]))
    for n in xrange(nsimulations):
        newyses.append(np.array([aux2[m][n] for m in xrange(len(aux2))]))

#======= DEBUG OPTIONS, comment if not necessary ==============================
#   plot the distribution of a single point to check if it follows a normal
#   distribution
    aheader = 'dispersao dos pontos falsos para um ponto da curva de luz 1 \n'
    np.savetxt('pointspread1.dat.gz', aux1[len(t)/2],
                delimiter=' ', header=aheader, comments='#')

    aheader = 'dispersao dos pontos falsos para um ponto da curva de luz 2 \n'
    np.savetxt('pointspread2.dat.gz', aux2[len(t)/2],
                delimiter=' ', header=aheader, comments='#')

    plt.hist(aux1[len(t)/2], label='mid point of curve 1', alpha=0.5)
    plt.hist(aux2[len(t)/2], label='mid point of curve 2', alpha=0.5)
    plt.title('Point spread distribution')
    plt.legend(loc='best')
    plt.savefig('pointspread.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

    aheader = 'Todas as curvas de luz falsas do grupo A \n'
    np.savetxt('fakecurves1.dat.gz', np.array(newxses).T,
                delimiter=' ', header=aheader, comments='#')

    aheader = 'Todas as curvas de luz falsas do grupo B \n'
    np.savetxt('fakecurves2.dat.gz', np.array(newxses).T,
                delimiter=' ', header=aheader, comments='#')

    aheader = 'coluna de tempo para plotar com as curvas dos arquivos "fakecurves" \n'
    np.savetxt('time.dat.gz', t,
                delimiter=' ', header=aheader, comments='#')

#   plot all fake points and original curve on top to check
    for simulated in newxses:
       plt.plot(t, simulated, '.', alpha=0.6)
    plt.errorbar(t, x, yerr=xe, fmt='k+-', linewidth=2.0)
    plt.title("Colored: randomized points, Black: Original lightcurve")
    plt.savefig('fake_x_lightcurves.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

#   plot 10 new x lightcurves and original on top to check
#    for simulated in newxses[:10]:
#       plt.plot(t, simulated, '.-', alpha=0.5)
#    plt.errorbar(t, x, yerr=xe, fmt='k+-', linewidth=2.0, alpha=0.7)
#    plt.title("Colored: 10 randomized lightcurves, Black: Original lightcurve")
#    plt.savefig('10_x_lightcurves.pdf', bbox_inches='tight', format='pdf',
#            papertype='a4', orientation='landscape')
#    plt.show()
#    plt.cla()

#   plot new y lightcurves and original on top to check
    for simulated in newyses:
        plt.plot(t, simulated, '.', alpha=0.6)
    plt.errorbar(t, y, yerr=ye, fmt='k+-', linewidth=2.0)
    plt.title("Colored: randomized lightcurves, Black: Original lightcurve")
    plt.savefig('fake_y_lightcurves.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

#   plot 10 new y lightcurves and original on top to check
#    for simulated in newyses[:10]:
#        plt.plot(t, simulated, '.-', alpha=0.5)
#    plt.errorbar(t, y, yerr=ye, fmt='k+-', linewidth=2.0, alpha=0.7)
#    plt.title("Colored: 10 randomized lightcurves, Black: Original lightcurve")
#    plt.savefig('10_y_lightcurves.pdf', bbox_inches='tight', format='pdf',
#            papertype='a4', orientation='landscape')
#    plt.show()
#    plt.cla()
#==============================================================================

#====== Statistical adjust of the simulations results =========================

#   calcule the various correlation functions
#   and stores the time-shifts to build a distribution
    shiftes = []
    corrs = []
    offsets = []
    for newx, newy in zip(newxses, newyses):
        newcorr, newoffset, nnewt, newshift = crosscorr(newx, newy, t)
        shiftes.append(newshift)
        corrs.append(newcorr)
        offsets.append(newoffset)

    corrs = corrs[::-1]
    corrs.append(newoffset)
    corrs = corrs[::-1]

    aheader = '''Todas as correlacoes entre os grupos A e B \n
                 - primeira coluna = delays, demais colunas = correlacao\n'''
    np.savetxt('allcorrelations.dat.gz', np.array(corrs).T,
                delimiter=' ', header=aheader, comments='#')

    for correlation, offset in zip(corrs[1:], offsets):
        plt.plot(offset, correlation, alpha=0.6)
    plt.xlabel('Offset [s]')
    plt.ylabel('Correlation')
    plt.savefig('correlations.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    query = str(raw_input('Plotar todas as correlacoes? (y/n): '))
    if query == 'y':
        plt.show()
    plt.cla()

#   calculates the natural binning of the lightcurve
    stp = abs(t[1] - t[0])

#   histogram bin size equals of the lightcurves bins
    binlist = np.arange(-max(t), max(t), step=stp)

#   a smaller bin size to check small flutuations of the distribution
    binlist2 = np.arange(-max(t), max(t), step=stp/10)

    aheader = 'Delays de todas as correlacoes entre os grups A e B \n'
    np.savetxt('alldelays.dat.gz', np.array(shiftes).T,
                delimiter=' ', header=aheader, comments='#')

#   plot original time shift distribution
    plt.hist(shiftes, bins=binlist, alpha=0.7)
    plt.hist(shiftes, bins=binlist2, alpha=0.5)
    plt.title('Distribution Function')
    plt.savefig('distribution_full.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

#   selected time shift limits for physical reasons if necessary
#   use min(shiftes) and max(shiftes) if not
    minoffset = float(raw_input('Enter Low limit for offset: '))
    maxoffset = float(raw_input('Enter High limit for offset: '))

#   newshifts = shiftes
    newshifts = [shiftes[i] for i in xrange(len(shiftes))
            if ((shiftes[i]>minoffset) and (shiftes[i]<maxoffset))]

#   fit normal distribution to the selected distribution
    mean, sigma = stats.norm.fit(newshifts)

    binlist = np.arange(minoffset, maxoffset, step=stp)
    binlist2 = np.arange(minoffset, maxoffset, step=stp/10)
#   plot selected time shift distribution
    plt.hist(newshifts, bins=binlist, normed=True, alpha=0.7)
    plt.hist(newshifts, bins=binlist2, normed=True, alpha=0.5)
#   create a x-axis for the gaussian funtion with 1000 points
    xgaus = np.linspace(minoffset, maxoffset, 10000)
#   generates the gaussian curve with mean and sigma
    gauss = stats.norm.pdf(xgaus, mean, sigma)
#   plot gaussian curve over histogram, with values on legend
    plt.plot(xgaus, gauss, color='k', linewidth=2.0,
            label='mean={0:.2f}, sigma={1:.2f}'.format(mean,sigma))
    plt.title('Selected Distribution Function')
    plt.legend(loc='best')
    plt.savefig('distribution_adjusted.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

# =============================================================================
#   Calculates correlation of x and y time series
    corr, offset, newt, shift = crosscorr(x, y, t)

    corrshift = shift

# === BEGIN of BLOCK ==========================================================
# == Comment this block to use results
#    free of monte-carlo statistics

#   time shift given by the maximum of the distribution
    shift = mean

#   new time axis to plot shifted time series
    newt = t + shift

#   correct time intervals if shift bigger than half the interval
    if min(newt) > (max(t)/2):
         newt = newt - max(t)
         shift = shift - max(t)
    elif max(newt) < (min(t)/2):
         newt = newt + min(t)
         shift = shift + min(t)
#=============================================== END of BLOCK =================

#   Print Info
    print ""
    print ""
    print "==================================================================="
    print "                                  INFO                             "
    print "==================================================================="
    print "Initial Time: {0}".format(tstart)
    print "Final Time:   {0}".format(tend)
    print "Total observed time: {0}".format(max(t))
    print "Temporal bin size: {0}".format(int(t[1]-t[0]))
    print "Number of bins:    {0}\n".format(len(x))
    print 'founded {0} negative or NaN values in the first lightcurve'.format(numexx)
    print 'list of index and values (swaped for zeros)'
    print excludedx
    print '\nfounded {0} negative or NaN values in the second lightcurve'.format(numexy)
    print 'list of index and values (swaped for zeros)'
    print excludedy
    print "\nResult from the direct cross-correlation function:"
    print "time shift = {0:.2f}\n".format(corrshift)
    print 'Results from the simulated distribution:\n'
    print '{0} Simulations'.format(nsimulations)
    print 'time shift = {0:.2f} +- {1:.2f}'.format(shift, sigma)
    print "\n=========================== END ================================="


    aheader = 'Correlacao entre as curvas 1 e 2 \n'
    np.savetxt('crosscorr.dat.gz', np.transpose([offset, corr]),
                delimiter=' ', header=aheader, comments='#')

#   plot correlation function
    plt.plot(offset, corr, 'o-')
    # position of chosen value
    plt.vlines(shift, min(corr), max(corr), 'k', 'dashed',
            'mean offset = {0:1f}'.format(shift))
    plt.xlabel('Offset [time units]', fontsize=12)
    plt.ylabel('Correlation coeficient', fontsize=12)
    plt.title('Correlation Function')
    plt.legend(loc='best')
    plt.savefig('crosscorr.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

#   plot original time series
    plt.errorbar(t, x, yerr=xe, label='series 1')
    plt.errorbar(t, y, yerr=ye, label='series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.legend(loc='best')
    plt.savefig('lightcurves.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()

#   plot original time series 1 plus shifted time series 2
    plt.errorbar(t, x, yerr=xe, label='series 1')
    plt.plot(newt, y, 'r', label='shifted series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.legend(loc='best')
    plt.savefig('new_lightcurves.pdf', bbox_inches='tight', format='pdf',
            papertype='a4', orientation='landscape')
    plt.show()
    plt.cla()
