# coding: utf-8

import sys
from bbcflib.bFlatMajor import common
from bbcflib.bFlatMajor.stream import concatenate
from bbcflib.btrack import FeatureStream
from numpy import float as nfloat, asarray, zeros

def _sum(scores,denom=None):
    return sum(scores)

def _arithmetic_mean(scores,denom):
    return sum(scores)*denom

def _geometric_mean(scores,denom):
## more precise/less efficient: exp(sum([log(x) for x in scores])*denom)
    return (reduce(lambda x, y: x*y, scores))**denom

def _min(scores,denom=None):
    return min(scores) if scores else 0

def _max(scores,denom=None):
    return max(scores) if scores else 0

def _qnth(vec, n):
    pivot = vec[0]
    below = [s for s in vec if s < pivot]
    above = [s for s in vec if s > pivot]
    i, j = len(below), len(vec)-len(above)
    if n < i:      return _qnth(below, n)
    elif n >= j:   return _qnth(above, n-j)
    else:          return pivot

def _median(scores,denom=None):
    if len(scores) % 2:
        return _qnth(scores,(len(scores)-1)/2)
    else:
        return (_qnth(scores,len(scores)/2-1)+_qnth(scores,len(scores)/2))*.5

_score_functions = {'arithmetic':_arithmetic_mean, 'geometric':_geometric_mean, 'sum':_sum,
                    'mean': _arithmetic_mean, 'min':_min, 'max':_max, 'median':_median}

@common.ordered
def merge_scores(trackList, method='arithmetic'):
    """
    Creates a stream with per-base average of several score tracks::

        X1: __________666666666______
        X2: _____2222222222__________
        R:  _____11111444443333______

    :param trackList: list of FeatureStream objects.
    :param method: (str) type of average: one of 'arithmetic','geometric', or 'sum' (no average).
        If one of the scores is 0, 1 is used instead for calculating the mean.
    :rtype: FeatureStream
    """
    tracks = [FeatureStream(common.sentinelize(x,[sys.maxint]*len(x.fields)), x.fields)
              for x in trackList]
    tracks = [common.reorder(t,['start','end','score']) for t in tracks]
    fields = [f for f in tracks[0].fields if all([f in t.fields for t in tracks])] # common fields
    elements = [list(x.next()) for x in tracks]
    track_denom = 1.0/len(trackList)

    mean_fn = _score_functions.get(method,_arithmetic_mean)
    for i in xrange(len(tracks)-1, -1, -1):
        if elements[i][0] == sys.maxint:
            tracks.pop(i)
            elements.pop(i)

    def _stream(tracks):
        while tracks:
            start = min([x[0] for x in elements])
            end = min([x[0] for x in elements if x[0]>start]+[x[1] for x in elements])
            scores = [x[2] for x in elements if x[1]>start and x[0]<end]
            if len(fields) > 3:
                rest = []
                for i in range(len(fields[3:])):
                    r = [str(x[3+i]) for x in elements if not(x[3+i] is None) and x[1]>start and x[0]<end]
                    if all([x == r[0] for x in r]):
                        rest.append(r[0])
                    else:
                        rest.append( "|".join(r) )
                yield (start, end, mean_fn(scores,track_denom)) + tuple(rest)
            else:
                yield (start, end, mean_fn(scores,track_denom))
            for i in xrange(len(tracks)-1, -1, -1):
                if elements[i][0] < end:
                    elements[i][0] = end
                if elements[i][1] <= end:
                    elements[i] = list(tracks[i].next())
                if elements[i][0] == sys.maxint:
                    tracks.pop(i)
                    elements.pop(i)
    return FeatureStream(_stream(tracks),fields)

###############################################################################
def filter_scores(trackScores,trackFeatures,method='sum',strict=False,annotate=False,flatten=common.cobble):
    """
    Extract from *trackScores* only the regions overlapping *trackFeatures*'s regions.
    Warning: both score and features streams must be sorted! (use `common.sorted_stream` if necessary).
    Example::

        X: _____#########__________#############_______
        Y: __________666666666___2222776_444___________
        R: __________6666__________22776_444___________

    Note: *trackFeatures* is :func:`cobbled <bbcflib.bFlatMajor.common.cobble>` by default (to avoid
    score duplications). An alternative is :func:`fusion <bbcflib.bFlatMajor.common.fusion>`, or nothing.
    If strand information is present in both *trackScores* and *trackFeatures*, only scores inside
    a region of the same strand are kept.

    :param trackScores: (FeatureStream) one -sorted- score track.
        If a list fo streams is provided, they will be merged (using `merge_scores`).
    :param trackFeatures: (FeatureStream) one -sorted- feature track.
        If a list fo streams is provided, they will be merged (using `concatenate`).
    :param method: (str) `merge_scores` *method* argument, in case *trackScores* is a list. ['sum']
    :param strict: (bool) if True, only score regions from *trackScores* that are
        strictly contained in a feature region of *trackFeatures* will be returned. [False]
    :param annotate: (bool) if True, supplementary annotation (and the corresponding fields)
        from *trackFeatures* will be added to the result. [False]
    :param flatten: (func) one of None, `common.fusion` or `common.cobble`.
        Function to be applied to *trackFeatures* before all. [common.cobble]
    :rtype: FeatureStream
    """
    def _stream(ts,tf,stranded):
        info_idx = [k for k,f in enumerate(tf.fields) if f not in ts.fields]
        tf = common.sentinelize(tf,[sys.maxint]*len(tf.fields))
        Y = [(-sys.maxint,-sys.maxint,0.0)]
        for x in ts:
            xstart = x[ts.fields.index('start')]
            xend = x[ts.fields.index('end')]
            ynext = Y[-1]
            # Load into Y all feature items which intersect score x
            while ynext[0] < xend:
                ynext = tf.next()
                if ynext[1] > xstart: Y.append(ynext)
            n = 0
            while Y[n][1] <= xstart: n+=1
            Y = Y[n:]
            for y in Y:
                if stranded and (x[ts.fields.index('strand')] != y[tf.fields.index('strand')]):
                    continue
                info = tuple([y[k] for k in info_idx]) if annotate else ()
                if strict and (y[0] > xstart or y[1] < xend): continue
                if y[0] >= xend : continue    # keep for next iteration
                start = xstart if y[0] < xstart else y[0]
                end   = xend   if y[1] > xend   else y[1]
                yield (start,end)+tuple(x[2:])+info

    if isinstance(trackFeatures,(list,tuple)): trackFeatures = concatenate(trackFeatures)
    if isinstance(trackScores,(list,tuple)): trackScores = merge_scores(trackScores,method)
    _info_fields = [f for f in trackFeatures.fields if f not in trackScores.fields] if annotate else []
    stranded = 'strand' in (set(trackScores.fields) & set(trackFeatures.fields))
    if flatten is None:
        _tf = trackFeatures
    else:
        _tf = flatten(trackFeatures,stranded=stranded)
    _ts = common.reorder(trackScores,['start','end'])
    return FeatureStream(_stream(_ts,_tf,stranded), _ts.fields+_info_fields)

###############################################################################
def score_by_feature(trackScores,trackFeatures,fn='mean'):
    """
    For every feature from *trackFeatures*, get the list of all scores it contains
    and apply an operation *fn* on this list (by default, scores are averaged).
    Warning: both score and feature streams must be sorted! (use `common.sorted_stream` is necessary).
    The output is a stream similar to *trackFeatures* but with an additional `score` field
    for each stream in *trackScores*::
        fn = 'mean':

        X: ------##########--------------##########------
        Y: ___________666666666__________6666666666______
        R: ______[   3.   ]______________[   6.   ]______


        fn = 'sum':

        X : ------##########--------------##########------
        Y1: ___________666666666__________6666666666______
        Y2: ___222222_____________________333_____________
        R : ______[  30,6  ]______________[  60,9  ]______

    :param trackScores: (list of) one or several -sorted- score track(s) (FeatureStream).
    :param trackFeatures: (FeatureStream) one -sorted- feature track.
    :param fn: (str of function): operation applied to the list of scores from one feature.
        Can be one of 'sum','mean','median','min','max', or a custom function.
    :rtype: FeatureStream
    """
    def _stream(ts,tf):
        X = [common.sentinelize(x, [sys.maxint]*len(x.fields)) for x in ts]
        S = [[(-sys.maxint,-sys.maxint,0.0)] for t in ts]
        if hasattr(fn,'__call__'):
            _fn = fn
        else:
            _fn = _score_functions.get(fn,_arithmetic_mean)
        for y in tf:
            ystart = y[tf.fields.index('start')]
            yend = y[tf.fields.index('end')]
            scores = ()
            for i in range(len(ts)):
                xnext = S[i][-1]
                # Load into S all score items which intersect feature y
                while xnext[0] < yend:
                    xnext = X[i].next()
                    if xnext[1] > ystart: S[i].append(xnext)
                n = 0
                while S[i][n][1] <= ystart: n+=1
                S[i] = S[i][n:]
                scores_y = []
                for s in S[i]:
                    if yend <= s[0]:   continue
                    if s[0] <  ystart: start = ystart
                    else:              start = s[0]
                    if yend <  s[1]:   end   = yend
                    else:              end   = s[1]
                    scores_y.extend([s[2]]*(end-start))
                scores += (_fn(scores_y,1.0/(yend-ystart)),)
            yield tuple(y)+scores

    if not(isinstance(trackScores,(list,tuple))): trackScores = [trackScores]
    if isinstance(trackFeatures,(list,tuple)):
        trackFeatures = concatenate(trackFeatures)
    if len(trackScores)>1:
        _fields = ["score"+str(i) for i in range(len(trackScores))]
    else:
        _fields = ["score"]
    _ts = [common.reorder(t,['start','end','score']) for t in trackScores]
    return FeatureStream(_stream(_ts,trackFeatures), trackFeatures.fields+_fields)

###############################################################################
def window_smoothing( trackList, window_size, step_size=1, stop_val=sys.maxint,
                      featurewise=False ):
    """
    Given a (list of) signal track(s) *trackList*, a *window_size* L (in base pairs by default,
    or in number of features if *featurewise* is True),  and a *step_size*,
    return as many signal tracks with, at each position p (multiple of *step_size*),
    the average score in the window [p-L/2, p+L/2]::

        X: __________666666666666____________
        R: ______12345666666666654321________ (not exact scores here)

    :param trackList: FeatureStream, or list of FeatureStream objects.
    :param window_size: (int) window size in bp.
    :param step_size: (int) step length. [1]
    :param stop_val: (int) sequence length. [sys.maxint]
    :param featurewise: (bool) bp (False), or number of features (True). [False]
    :rtype: FeatureStream

    Example of windows, window_size=9, step_size=3:

    [0,1,2,3,4,5,6,7,8,9), [3,4,5,6,7,8,9,10,11,12), ...
    """
    def _stepping_mean(track,score,denom):
        score = 0.0
        F = []
        score = 0.0
        nmid = window_size/2
        for x in track:
            F.append(x)
            score += x[2]
            if len(F)<window_size: continue
            score -= F.pop(0)[2]
            yield (F[nmid][0],F[nmid][1],score*denom)+F[nmid][3:]
            for shift in xrange(2,step_size):
                score -= F.pop(0)[2]

    def _running_mean(track,win_start,denom):
        score = 0.0
        F = []
        for x in track:
            F.append(x)
            fstart = F[0][0]
            fend = F[0][1]
            win_start = max(win_start,fstart-window_size)
            win_end = win_start+window_size
            lstart = F[-1][0]
            lend = F[-1][1]
            while win_end < lend:
                delta = 0
                steps = [fend-win_start,lend-win_end]
                if fstart>win_start: steps.append(fstart-win_start)
                else: delta -= F[0][2]
                if lstart>win_end:   steps.append(lstart-win_end)
                else: delta += F[-1][2]
                nsteps = min(steps)
                sst = -(win_start % step_size)%step_size
                sen = -(win_start+nsteps % step_size)%step_size
                win_center = (win_start+win_end)/2
                if abs(delta) > 1e-11:
                    delta *= denom
                    score += delta*sst
                    for step in xrange(sst,nsteps,step_size):
                        if score>1e-11 and win_center+step>=0 and win_center+step+step_size<=stop_val:
                            yield (win_center+step,win_center+step+step_size,score)
                        score += delta*step_size
                    score -= delta*sen
                else:
                    if score>1e-11 and win_center+sst>=0 and win_center+sen+nsteps<=stop_val:
                        yield (win_center+sst,win_center+sen+nsteps,score)
                win_start += nsteps
                win_end += nsteps
                if fend <= win_start:
                    F.pop(0)
                    if F:
                        fstart = F[0][0]
                        fend = F[0][1]
                        win_start = max(win_start,fstart-window_size)
                        win_end = win_start+window_size
                        if win_end > stop_val: break
        while F:
            delta = 0
            steps = [fend-win_start]
            if fstart>win_start: steps.append(fstart-win_start)
            else: delta -= F[0][2]
            nsteps = min(steps)
            sst = -(win_start % step_size)%step_size
            sen = -(win_start+nsteps % step_size)%step_size
            win_center = (win_start+win_end)/2
            if abs(delta) > 1e-11:
                delta *= denom
                score += delta*sst
                for step in xrange(sst,nsteps,step_size):
                    if score>1e-11 and win_center+step>=0 and win_center+step+step_size<=stop_val:
                        yield (win_center+step,win_center+step+step_size,score)
                    score += delta*step_size
                score -= delta*sen
            else:
                if score>1e-11 and win_center+sst>=0 and win_center+sen+nsteps<=stop_val:
                    yield (win_center+sst,win_center+sen+nsteps,score)
            win_start += nsteps
            win_end += nsteps
            if fend <= win_start:
                F.pop(0)
                if F:
                    fstart = F[0][0]
                    fend = F[0][1]
                    win_start = max(win_start,fstart-window_size)
                    win_end = win_start+window_size
                    if win_end > stop_val: break

    denom = 1.0/window_size
    win_start = -window_size
    _f = ['start','end','score']
    if featurewise:
        call = _stepping_mean
    else:
        call = _running_mean
    if isinstance(trackList,(list,tuple)):
        return [FeatureStream(call(common.reorder(t,_f),win_start,denom),fields=_f)
                for n,t in enumerate(trackList)]
    else:
        return FeatureStream(call(common.reorder(trackList,_f),win_start,denom),fields=_f)

###############################################################################
def normalize(trackList,method='total',field='score'):
    """Normalizes the scores in every stream from *trackList* using the given *method*.
    It assumes that each of the streams represents the same features, i.e. the n-th element
    of one stream corresponds to the n-th element of another.

    [!] This function will temporarily store everything in memory.

    :param trackList: FeatureStream, or list of FeatureStream objects.
    :param method: normalization method:
        * ``'total'`` divides every score vector by its sum (total number of reads) x 10^7 .
        * ``'deseq'`` applies DESeq's normalization ("size factors") - considering every track
            as belonging to a different group.
        * ``'quantile'`` applies quantile normalization.
    :param field: (str) name of the field containing the scores (must be the same for all streams).
    """
    if not isinstance(trackList,(list,tuple)):
        trackList = [trackList]
    allcontents = [list(t) for t in trackList]
    ncols = len(trackList)
    nlines = len(allcontents[0])
    assert all(len(t)==nlines for t in allcontents), "All streams must have the same number of elements."
    # Build the matrix
    allscores = zeros((ncols,nlines))
    for n,content in enumerate(allcontents):
        idx = trackList[n].fields.index(field)
        allscores[n] = asarray([x[idx] for x in content])
    # Normalize
    allscores = common.normalize(asarray(allscores),method)
    # Reinsert the new scores in the respective tracks
    for n,content in enumerate(allcontents):
        idx = trackList[n].fields.index(field)
        for k,x in enumerate(content):
            content[k] = x[:idx] + (allscores[n][k],) + x[idx+1:]
    res = [FeatureStream(t,fields=trackList[n].fields) for n,t in enumerate(allcontents)]
    if len(trackList) == 1:
        return res[0]
    else:
        return res

