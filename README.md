cross-correlation
=================

A python implementation of a cross-correlation task that finds time delays 
between two time series, with monte-carlo simulations to estimate the
uncertainties.

Assuming that the error of the time series are 1-sigma deviations of the value,
generates several fake curves with random points that follow a normal
distribution with the same 1-sigma deviation. Calculates the cross-correlation
function and time-delay for each fake curve. Find the mean and standard
deviation of the distribution of time-delays.

**External packages needed:**
  - matplotlib
  - numpy
  - astropy
  - scipy
