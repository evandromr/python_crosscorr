#!/usr/env python

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
    lags = np.arange(corr.size) - (t.size - 1)
    tstep = (t[-1] - t[0])/float(t.size)
    offset = lags*tstep

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
    print ''
    print 'Nan and negative values = 0'
    print ''
    for i in xrange(len(x)):
        if (x[i] >=0):
            pass
        else:
            x[i] = 0
            xe[i] = 0
        
        if (y[i] >=0):
            pass
        else:
            y[i]=0
            ye[i] = 0
#------------------------------------------------------------------------------

##  start time from Zero
    t -= min(t)

### === MONTE CARLO SIMULATIONS ===============================================
    nsimulations = int(raw_input("How many simulations?: "))

#   for each point in x and y, generates 'nsimulation' new points
#   new value = (uncertainty*random.values + mean)
#   the randon files follow a normal distribution with 1-sigma equals the
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

#======= DEBUG OPTIONS, may comment this part =================================
#   plot the distribution of a single point to check if it follows a normal
#   distribution
    plt.hist(aux1[len(t)/2], label='mid point of curve 1', alpha=0.5)
    plt.hist(aux2[len(t)/2], label='mid point of curve 2', alpha=0.5)
    plt.title('Point spread distribution')
    plt.legend(loc='best')
    plt.savefig('pointspread.ps')
    plt.show()
    plt.cla()

#   plot 10 new x lightcurves and original on top to check
    for simulated in newxses[:10]:
       plt.plot(t, simulated, '.-', alpha=0.5)
    plt.errorbar(t, x, yerr=xe, fmt='k+-', linewidth=2.0, alpha=0.7)
    plt.title("Colored: 10 randomized lightcurves, Black: Original lightcurve")
    plt.savefig('10_x_lightcurves.ps')
    plt.show()
    plt.cla()

#   plot new y lightcurves and original on top to check
    for simulated in newyses[:10]:
        plt.plot(t, simulated, '.-', alpha=0.5)
    plt.errorbar(t, y, yerr=ye, fmt='k+-', linewidth=2.0, alpha=0.7)
    plt.title("Colored: 10 randomized lightcurves, Black: Original lightcurve")
    plt.savefig('10_y_lightcurves.ps')
    plt.show()
    plt.cla()
#==============================================================================

#====== Statistical adjust of the simulations results =========================

#   calcule the various correlation functions
#   and stores the time-shifts to build a distribution
    shiftes = []
    for newx, newy in zip(newxses, newyses):
        newcorr, newoffset, nnewt, newshift = crosscorr(newx, newy, t)
        shiftes.append(newshift)

#   calculates the natural binning of the lightcurve
    stp = abs(t[1] - t[0])

#   histogram bin size equals of the lightcurves bins
    binlist = np.arange(-max(t), max(t), step=stp)

#   a smaller bin size to check small flutuations of the distribution
    binlist2 = np.arange(-max(t), max(t), step=stp/10)

#   plot original time shift distribution
    plt.hist(shiftes, bins=binlist, alpha=0.7)
    plt.hist(shiftes, bins=binlist2, alpha=0.5)
    plt.title('Distribution Function')
    plt.savefig('distribution_full.ps')
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
    plt.savefig('distribution_adjusted.ps')
    plt.show()
    plt.cla()

# =============================================================================
#   Calculates correlation of x and y time series
    corr, offset, newt, shift = crosscorr(x, y, t)

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

#   visualize calculated time shift
    print 'time shift = {0:.2f} +- {1:.2f}'.format(shift, sigma)

#   open file to write results of the correlation function
    out = open('crosscorr.dat', 'w')
#   write offset and corr to file 'crosscorr.dat' in 2 columns
    for i in xrange(len(corr)):
        out.write('{0} {1} \n'.format(offset[i], corr[i]))
    out.close()

#   plot correlation function
    plt.plot(offset, corr, 'o-')
    # position of chosen value
    plt.vlines(shift, min(corr), max(corr), 'k', 'dashed',
            'mean offset = {0:1f}'.format(shift))
    plt.xlabel('Offset [time units]', fontsize=12)
    plt.ylabel('Correlation coeficient', fontsize=12)
    plt.title('Correlation Function')
    plt.legend(loc='best')
    plt.savefig('crosscorr.ps')
    plt.show()
    plt.cla()

#   plot original time series
    plt.errorbar(t, x, yerr=xe, label='series 1')
    plt.errorbar(t, y, yerr=ye, label='series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.savefig('lightcurves.ps')
    plt.show()
    plt.cla()

#   plot original time series 1 plus shifted time series 2
    plt.errorbar(t, x, yerr=xe, label='series 1')
    #plt.plot(t, y, label='series 2')
    plt.plot(newt, y, 'r', label='shifted series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.legend(loc='best')
    plt.savefig('new_lightcurves.ps')
    plt.show()
    plt.cla()
