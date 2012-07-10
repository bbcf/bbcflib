from bbcflib.bFlatMajor import common
import numpy
#from scipy.fftpack import fft, ifft
from numpy.fft import fft, ifft
from numpy import conjugate
from math import log

def normalize(x):
    """Substracts the average and divides by the standard deviation."""
    mu = numpy.mean(x)
    isigma = 1.0/numpy.sqrt((x*x).mean()-mu*mu)
    return (x-mu)*isigma

def correlation(trackList, start, end, limits=None):
    """
    Calculates the cross-correlation between two tracks Q1,Q2. Returns a vector
    containing the correlation at each lag in this order (L/R for resp. min/max lag value):
    [-L,-L+1,...-1,0,1,...,R-1,R]. If more than two tracks are given in *trackList*,
    returns a list of correlation vectors, one for every distinct pair of tracks.

    A positive lag indicates that Q2 is shifted to the right w.r.t Q1,
    a negative lag - to the left.
    So to get the correlation at lag +4, one has to look at the L+4 th
    element of the array (in Python notation - or L+5 counting from 1);
    for lag -4, at the L-4 rd element (in Python notation - or L-3 counting from 1).

    Example::

        |_____ /^\ _________|         lag 0
        |______________/^\__|


        |_____ /^\ _________|         lag +8
           ->   |______________/^\__|


                |_____ /^\ _________| lag -8
        |______________/^\__|  <-

    :param trackList: list of FeatureStream objects of the same size.
    :param start,end: (int) bounds of the region to consider, in bp.
    :param limits: (tuple (int,int)) maximum lag to consider. [(-len(Q1)+1,len(Q1)-1)]
    :rtype: list of floats, or list of lists of floats.
    """
    ##### One could profit from numpy to reduce the memory space used for
    ##### storing these - long - arrays ('dtype' is float64 by default).
    x = [numpy.array([s[0] for s in common.unroll(t,start,end)]) for t in trackList]
    x = [normalize(t) for t in x]
    if not limits:
        limits = (-len(x[0])+1 , len(x[0])-1)
    N = len(x[0])+limits[1]-limits[0]-1
    ##### convert to nearest power of 2, fft gets orders of magnitude faster...
    N = 2**int(log(2+N,2)+.5)
    def _corr(x1,x2,N):
        corr = ifft(fft(x1,N)*conjugate(fft(x2,N)))#/len(x1)
        corr = numpy.concatenate((corr[N+limits[0]:], corr[:limits[1]+1]))
        return numpy.real(corr)
    if len(trackList) == 2:
        return _corr(x[0],x[1],N)
    else:
        return [_corr(x1,x2,N) for n,x1 in enumerate(x) for x2 in x[n+1:]]

################################################################################
