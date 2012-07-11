from bbcflib.bFlatMajor import common
import numpy
try:
    from scipy.fftpack import fft, ifft
except ImportError:
    from numpy.fft import fft, ifft
from numpy import conjugate
from math import log

def normalize(x):
    """Substracts the average and divides by the standard deviation."""
    x = numpy.asarray(x)
    mu = numpy.mean(x)
    isigma = 1.0/numpy.sqrt((x*x).mean()-mu*mu)
    return (x-mu)*isigma

def correlation(trackList, start, end, limits=(-1000,1000), with_acf=False):
    """
    Calculates the cross-correlation between two streams and
    returns a vector containing the correlation at each lag in this order 
    (L/R for resp. limits[0], limits[1]):
    [L,L+1,...,R-1,R]. 
    If more than two tracks are given in *trackList*,
    returns a list of correlation vectors, one for every distinct pair of tracks.
    If `with_acf` is True, self-correlations will also be included in the list.

    A negative lag indicates that track 2 is shifted to the right w.r.t track 1,
    a positive lag - to the left.
    So to get the correlation at lag +4, one has to look at the 4-L th
    element of the array;
    for lag -4, at the -4-L th element.

    Example::

        |_____ /^\ _________|         lag 0
        |______________/^\__|


        |_____ /^\ _________|         lag -8
           ->   |______________/^\__|


                |_____ /^\ _________| lag +8
        |______________/^\__|  <-

    :param trackList: list of FeatureStream objects
    :param start,end: (int) bounds of the region to consider, in bp.
    :param limits: (tuple (int,int)) maximum lag to consider. [-1000,1000]
    :param with_acf: (bool) include auto-correlations. [False]
    :rtype: list of floats, or list of lists of floats.
    """
    ##### One could profit from numpy to reduce the memory space used for
    ##### storing these - long - arrays ('dtype' is float64 by default).
    x = [numpy.array([s[0] for s in common.unroll(t,start,end)]) for t in trackList]
    x = [normalize(t) for t in x]
    if limits[1]-limits[0] > 2*len(x[0]):
        limits = (-len(x[0])+1,len(x[0])-1)
    N = len(x[0])+limits[1]-limits[0]-1
    ##### convert to nearest power of 2, fft gets orders of magnitude faster...
    N = 2**int(log(2+N,2)+.5)
    def _corr(x1,x2,N):
        corr = ifft(fft(x1,N)*conjugate(fft(x2,N)))/len(x1)
        corr = numpy.concatenate((corr[N+limits[0]:], corr[:limits[1]+1]))
        return numpy.real(corr)[::-1]
    if with_acf:
        return [[_corr(x1,x2,N) for x2 in x[n:]] for n,x1 in enumerate(x)]
    elif len(trackList) == 2:
        return _corr(x[0],x[1],N)
    else:
        return [[_corr(x1,x2,N) for x2 in x[n+1:]] for n,x1 in enumerate(x[:-1])]

################################################################################
