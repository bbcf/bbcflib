from bbcflib.bFlatMajor import common
from bbcflib import btrack as track
import numpy
from scipy.fftpack import fft, ifft
from numpy import conjugate
from math import log

def normalize(x):
    """Substracts the average and divides by the standard deviation."""
    mu = numpy.mean(x)
    isigma = 1.0/numpy.sqrt((x*x).mean()-mu*mu)
    return (x-mu)*isigma

def correlation(trackList, start, end, limits=(-1000,1000)):
    """
    Calculates the cross-correlation between two tracks Q1,Q2. Returns a vector
    containing the correlation at each lag in this order (L for max lag value):
    [-L,-L+1,...-1,0,1,...,L-1,L].

    CHECK THAT THE FOLLOWING STILL HOLDS:

    A negative lag indicates that Q2 is shifted to the right w.r.t Q1,
    a positive lag - to the left.
    So to get the correlation at lag +4, one has to look at the -L+5 th
    element of the array.

    Example:

    |_____ /^\ _________|         lag 0
    |______________/^\__|


    |_____ /^\ _________|         lag -8
       ->   |______________/^\__|


            |_____ /^\ _________| lag +8
    |______________/^\__|  <-

    :param trackList: list of **two** FeatureStream objects.
    :param start,end: (int) bounds of the region to consider, in bp.
    :param limits: (tuple (int,int)) maximum lag to consider. [(-1000,1000)]
    :rtype: list of floats
    """
    # One could profit from numpy to reduce the memory space used for
    # storing these - long - arrays ('dtype')
    x = [numpy.array([s[0] for s in common.unroll(t,start,end,'score')])
         for t in trackList]
    x = [normalize(t) for t in x]
    N = len(x[0])+limits[1]-limits[0]-1
    ##### convert to nearest power of 2, fft gets orders of magnitude faster...
    N = 2**int(log(2+N,2)+.5)
    def _corr(x1,x2,N):
        corr = ifft(conjugate(fft(x1,N))*fft(x2,N))/len(x1)
        corr = numpy.concatenate((corr[N+limits[0]:], corr[:limits[1]+1]))
        return numpy.real(corr)
    return [_corr(x1,x2,N) for n,x1 in enumerate(x) for x2 in x[(n+1):]]

################################################################################
