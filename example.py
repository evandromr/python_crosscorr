#!/usr/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal, stats
import pyfits as fits


def corrfunc(x, y, t):
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
    Creates 2 fake time series, and calculates the cross correlation and time
    delay between then
    '''

#   Time Series
    T = float(raw_input('Insert a period for the time series: '))
    P = float(raw_input('Insert the total sampled time: '))
    stp = float(raw_input('Insert the time step (sample rate): '))

#   time delay
    delay = float(raw_input('Insert the time delay between the 2 series: '))

#   1-sigma errors
    sigma = float(raw_input('Insert the fake 1-sig deviation for the data: '))

    # time data
    t = np.arange(T, P, step=stp)

    t -= min(t)
    # sinusoidal time series
    x = np.sin((2.0*np.pi*t)/T)
    y = np.sin(((2.0*np.pi*(t - delay))/T))

    x = x + sigma*np.random.randn(len(x))
    y = y + sigma*np.random.randn(len(y))

#   number of simulations
    nsimulations = int(raw_input('Insert the number of simulations :'))

#   generates 'nsimulations' fake time series
    aux1 = []
    aux2 = []
    for i, meanx in enumerate(x):
        newx = 0.3*np.random.randn(nsimulations) + meanx
        aux1.append(newx)
    for j, meany in enumerate(y):
        newy = 0.3*np.random.randn(nsimulations) + meany
        aux2.append(newy)

    newxses = []
    newyses = []
    for n in xrange(nsimulations):
        newxses.append(np.array([aux1[m][n] for m in xrange(len(aux1))]))
    for n in xrange(nsimulations):
        newyses.append(np.array([aux2[m][n] for m in xrange(len(aux2))]))


#======= DEBUG OPTION ==================================================
#   plot new x lightcurves and original on top to check
    for simulated in newxses:
       plt.plot(t, simulated, '.')
    plt.errorbar(t, x, yerr=0.2, fmt='k+-', linewidth='2.0')
    plt.show()
    plt.cla()

#   plot new y lightcurves and original on top to check
    for simulated in newyses:
        plt.plot(t, simulated, '.')
    plt.errorbar(t, y, yerr=0.2, fmt='k+-', linewidth='2.0')
    plt.show()
    plt.cla()
#=======================================================================

#   store calculated time shift for each simulated curve
    shiftes = []
    for newx, newy in zip(newxses, newyses):
        newcorr, newoffset, nnewt, newshift = corrfunc(newx, newy, t)
        shiftes.append(newshift)

#   histogram binning equal of time step
    binlist = np.arange(-max(t), max(t), step=stp)

#   plot original time shift distribution
    plt.hist(shiftes, bins=binlist, normed=True, alpha=0.6)
    plt.title('Distribution Function')
    plt.show()
    plt.cla()

#   histogram binnin manually defined (step)
    binlist2 = np.arange(-max(t), max(t), step=5)

#   plot original time shift distribution
#    plt.hist(shiftes, bins=binlist2, normed=True, alpha=0.6)
#    plt.title('Distribution Function')
#    plt.show()
#    plt.cla()

#   calculates the mean and sigma of original distribution (without selection)
    mean, sigma = stats.norm.fit(shiftes)
    print 'Results from the total distribution (without selection)'
    print 'time shift = {0:.2f} +- {1:.2f}'.format(mean, sigma)
    print ' '

#   selected time shift limits for physical reasons
#   use min(shiftes) and max(shiftes) if not
    minoffset = float(raw_input('Enter Low limit for offset: '))
    maxoffset = float(raw_input('Enter High limit for offset: '))

#   newshifts = shiftes
    newshifts = [shiftes[i] for i in xrange(len(shiftes))
            if ((shiftes[i]>minoffset) and (shiftes[i]<maxoffset))]

#   fit normal distribution
    mean, sigma = stats.norm.fit(newshifts)

#   histogram binning equals of time step
    binlist = np.arange(minoffset, maxoffset, step=stp)

#   smaller histogram bin, set mannually
    binlist2 = np.arange(minoffset, maxoffset, step=1)

#   plot selected time shift distribution
    plt.hist(newshifts, bins=binlist, normed=True, alpha=0.6)
    plt.hist(newshifts, bins=binlist2, normed=True, alpha=0.6)

#   create a x-axis for the gaussian funtion with 1000 points
    xgaus = np.linspace(minoffset, maxoffset, 10000)

#   generates the gaussian curve with mean and sigma
    gauss = stats.norm.pdf(xgaus, mean, sigma)

#   plot gaussian curve over histogram, with values on legend
    plt.plot(xgaus, gauss, color='k', linewidth=2.0,
            label='mean={0:.2f}, sigma={1:.2f}'.format(mean,sigma))
    plt.title('Selected Distribution Function')
    plt.legend(loc='best')
    plt.show()
    plt.cla()

# =========================================================================
#   Calculates correlation of x and y time series
    corr, offset, newt, shift = corrfunc(x, y, t)

# === BEGIN of BLOCK ======================================================
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
#=============================================== END of BLOCK ==============

#   visualize calculated time shift
    print 'results from the selected distribution'
    print 'time shift = {0:.2f} +- {1:.2f}'.format(shift, sigma)

#   open file to write the cross correlation function
    out = open('crosscorr.dat', 'w')
#   write offset and corr to file 'crosscorr.dat' in 2 columns
    for i in xrange(len(corr)):
        out.write('{0} {1} \n'.format(offset[i], corr[i]))
    out.close()

#   plot correlation function
    plt.plot(offset, corr, 'o-')
    # position of maximum chosen value
    plt.vlines(shift, min(corr), max(corr), 'k', 'dashed',
            'mean offset = {0:1f}'.format(shift))
    plt.xlabel('Offset [time units]', fontsize=12)
    plt.ylabel('Correlation coeficient', fontsize=12)
    plt.title('Correlation Function')
    plt.legend(loc='best')
    plt.show()
    plt.cla()

#   plot original time series
    plt.plot(t, x, label='series 1')
    plt.plot(t, y, label='series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.show()
    plt.cla()

#   plot original time series plus shifted time series
    plt.plot(t, x, label='series 1')
    #plt.plot(t, y, label='series 2')
    plt.plot(newt, y, 'r', label='shifted series 2')
    plt.xlabel('Time [s]', fontsize=12)
    plt.ylabel('Normalized Count Rate [counts/s]', fontsize=12)
    plt.legend(loc='best')
    plt.show()
    plt.cla()
