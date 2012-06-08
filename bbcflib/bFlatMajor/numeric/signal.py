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
    """
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
