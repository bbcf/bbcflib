from bbcflib.gfminer.common import unroll
from bbcflib.track import FeatureStream
from numpy.fft import fft, ifft
from numpy import conjugate,array,asarray,mean,sqrt,real,hstack
from numpy import concatenate as ncat
from math import log

def score_array(trackList,fields=['score']):
    """Returns a numeric array with the *fields* columns from each input track
    and a vector of row labels, taken from the *name* field which must match in all tracks."""
    if not(isinstance(trackList,(list,tuple))): trackList=[trackList]
    sidx = [[t.fields.index(f) for f in fields] for t in trackList]
    nidx = [t.fields.index('name') for t in trackList]
    dico = dict((k[nidx[0]],[k[i] for i in sidx[0]]) for k in trackList[0])
    nums = asarray(dico.values())
    labs = asarray(dico.keys())
    for n,tn in enumerate(trackList[1:]):
        dico = dict((k[nidx[n+1]],[k[i] for i in sidx[n+1]]) for k in tn)
        nums = hstack((nums,[dico[k] for k in labs]))
    return (nums,labs)

def _normalize(x):
    """Substracts the average and divides by the standard deviation."""
    x = asarray(x)
    if any(abs(x-x[0])>1e-6):
        mu = mean(x)
        isigma = 1.0/sqrt((x*x).mean()-mu*mu)
        return (x-mu)*isigma
    else:
        return x-x[0]

def correlation(trackList, regions, limits=(-1000,1000), with_acf=False):
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
    :param regions: a tuple (start,end) or a FeatureStream with the bounds of the regions to consider (see `unroll`).
        In the latter case, all regions will be concatenated.
    :param limits: (tuple (int,int)) maximum lag to consider. [-1000,1000]
    :param with_acf: (bool) include auto-correlations. [False]
    :rtype: list of floats, or list of lists of floats.
    """
    ##### One could profit from numpy to reduce the memory space used for
    ##### storing these - long - arrays ('dtype' is float64 by default).
    if isinstance(regions,FeatureStream):
        _reg = list(regions)
        x = [array([s[0] for s in unroll(t,FeatureStream(_reg,fields=regions.fields))])
             for t in trackList]
        _reg = []
    else:
        x = [array([s[0] for s in unroll(t,regions)]) for t in trackList]
    x = [_normalize(t) for t in x]
    if limits[1]-limits[0] > 2*len(x[0]):
        limits = (-len(x[0])+1,len(x[0])-1)
    N = len(x[0])+limits[1]-limits[0]-1
    ##### convert to nearest power of 2, fft gets orders of magnitude faster...
    N = 2**int(log(2+N,2)+.5)
    def _corr(x1,x2,N):
        corr = ifft(conjugate(fft(x1,N))*fft(x2,N))/len(x1)
        corr = ncat((corr[N+limits[0]:], corr[:limits[1]+1]))
        return real(corr)
    if with_acf:
        return [[_corr(x1,x2,N) for x2 in x[n:]] for n,x1 in enumerate(x)]
    elif len(trackList) == 2:
        return _corr(x[0],x[1],N)
    else:
        return [[_corr(x1,x2,N) for x2 in x[n+1:]] for n,x1 in enumerate(x[:-1])]

################################################################################
