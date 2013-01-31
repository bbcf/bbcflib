from bbcflib.bFlatMajor.common import unroll
from bbcflib.btrack import FeatureStream
from numpy.fft import fft, ifft
from numpy import conjugate,array,asarray,mean,median,sqrt,exp,real,hstack,nonzero,prod,around,argsort
from numpy import log as nlog,concatenate as ncat, float as nfloat
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

def normalize(trackList,method='total',field='score'):
    """Normalizes the scores in every stream from *trackList* using the given *method*.
    It assumes that each of the streams represents the same features, i.e. the n-th element
    of one stream corresponds to the n-th element of another

    [!] This function will temporarily store everything in memory.

    :param method: normalization method:
        * ``'total'`` divides every score vector by its sum (total number of reads) x 10^7 .
        * ``'deseq'`` applies DESeq's normalization ("size factors") - considering every track
            as belonging to a different group.
        * ``'quantile'`` applies quantile normalization.
    :param field: (str) name of the field containing the scores (must be the same for all streams).
    """
    def _total(scores):
        for n,col in enumerate(scores):
            scores[n] = scores[n]/sum(scores[n])
        return scores
    def _deseq(scores):
        cnts = scores[:,nonzero(prod(scores,axis=0))[0]] # none of the counts is zero
        loggeomeans = mean(nlog(cnts),axis=0) # -inf if division by 0
        size_factors = exp(median(nlog(cnts)-loggeomeans,axis=1))
        scores = around((scores.T / size_factors).T,2)
        return scores
    def _quantile(scores):
        ordering = argsort(scores)
        for n in range(len(scores)):
            scores[n] = scores[n][ordering[n]]
        means = mean(scores,0)
        for n in range(len(scores)):
            scores[n] = around(means[argsort(ordering[n])],2)
        return scores

    if method == 'total': f = _total
    elif method == 'deseq': f = _deseq
    elif method == 'quantile': f = _quantile
    elif hasattr(method,'__call__'): f = method # custom function
    else: raise ValueError("Unknown normalization method (got %s)" % method)
    if isinstance(trackList,(list,tuple)):
        allcontents = [None]*len(trackList)
        allscores = [None]*len(trackList)
        for n,t in enumerate(trackList):
            # make a copy of the whole content, and scores separately
            score_idx = t.fields.index(field)
            allcontents[n] = []
            allscores[n] = []
            for x in t:
                allcontents[n].append(list(x))
                allscores[n].append(x[score_idx])
        allscores = f(asarray(allscores))
        for n,content in enumerate(allcontents):
            score_idx = trackList[n].fields.index(field)
            for k,x in enumerate(content):
                x[score_idx] = allscores[n][k]
                content[k] = tuple(x)
        return [FeatureStream(c,fields=trackList[n].fields) for n,c in enumerate(allcontents)]
    elif method == 'total' or method == f:
        content = []
        scores = []
        score_idx = trackList.fields.index(field)
        for x in trackList:
            content.append(list(x))
            scores.append(x[score_idx])
        scores = f(asarray([scores,],dtype=nfloat))
        for k,x in enumerate(content):
            x[score_idx] = scores[0][k]
            content[k] = tuple(x)
        return FeatureStream(content,fields=trackList.fields)
    else:
        return trackList

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
