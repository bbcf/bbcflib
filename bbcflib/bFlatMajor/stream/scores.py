# coding: utf-8

import sys
from math import floor, ceil
from bbcflib.bFlatMajor import common
from bbcflib.bFlatMajor.stream import concatenate
from bbcflib import btrack as track

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
    tracks = [track.FeatureStream(common.sentinelize(x,[sys.maxint]*len(x.fields)), x.fields)
              for x in trackList]
    tracks = [common.reorder(t,['start','end','score']) for t in tracks]
    fields = [f for f in tracks[0].fields if all([f in t.fields for t in tracks])] # common fields
    elements = [list(x.next()) for x in tracks]
    track_denom = 1.0/len(trackList)

    def _sum(scores,denom):
        return sum(scores)
    def _arithmetic_mean(scores,denom):
        return sum(scores)*denom
    def _geometric_mean(scores,denom):
        return (reduce(lambda x, y: x*y, scores))**denom
        ## more precise/less efficient: exp(sum([log(x) for x in scores])*denom)
    methods = {'arithmetic':_arithmetic_mean, 'geometric':_geometric_mean, 'sum':_sum}
    mean_fn = methods[method]
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
    return track.FeatureStream(_stream(tracks),fields)

###############################################################################
def filter_scores(trackScores,trackFeatures,method='sum'):
    """
    Extract from *trackScores* only the regions present in *trackFeatures*.
    Example::

        X: _____#########__________#############_______
        Y: __________666666666___2222776_444___________
        R: __________6666__________22776_444___________

    :param trackScores: (FeatureStream) one score track.
        If a list fo streams is provided, they will be merged (averaged scores).
    :param trackFeatures: (FeatureStream) one feature track.
        If a list fo streams is provided, they will be merged.
    :param method: (str) `merge_scores` *method* argument. ['sum']
    :rtype: FeatureStream
    """
    def _stream(ts,tf):
        X = common.sentinelize(ts, [sys.maxint]*len(ts.fields))
        S = [(-sys.maxint,-sys.maxint,0.0)]
        for y in tf:
            ystart = y[tf.fields.index('start')]
            yend = y[tf.fields.index('end')]
            xnext = S[-1]
            # Load into S all score items which intersect feature y
            while xnext[0] < yend:
                xnext = X.next()
                if xnext[1] > ystart: S.append(xnext)
            n = 0
            while S[n][1] <= ystart: n+=1
            S = S[n:]
            for s in S:
                if yend <= s[0]: continue
                start = ystart if s[0] < ystart else s[0]
                end   = yend   if yend < s[1]   else s[1]
                yield (start,end)+tuple(s[2:])

    if isinstance(trackFeatures,(list,tuple)): trackFeatures = concatenate(trackFeatures)
    if isinstance(trackScores,(list,tuple)): trackScores = merge_scores(trackScores,method)
    _ts = common.reorder(trackScores,['start','end'])
    return track.FeatureStream(_stream(_ts,trackFeatures), _ts.fields)

###############################################################################
def score_by_feature(trackScores,trackFeatures,fn='mean'):
    """
    For every feature from *trackFeatures*, get the list of all scores it contains
    and apply an operation *fn* on this list (by default, scores are averaged).
    The output is a stream similar to *trackFeatures* but with an additional `score` field
    for each stream in *trackScores*::

        X: ------##########--------------##########------
        Y: ___________666666666__________6666666666______
        R: ______[   3.   ]______________[   6.   ]______


        normalize = False:

        X : ------##########--------------##########------
        Y1: ___________666666666__________6666666666______
        Y2: ___222222_____________________333_____________
        R : ______[  30,6  ]______________[  60,9  ]______

    :param trackScores: (list of) one or several score track(s) (FeatureStream).
    :param trackFeatures: (FeatureStream) one feature track.
    :param fn: (str of function): operation applied to the list of scores from one feature.
        Can be one of 'sum','mean','median','min','max', or a custom function
    :rtype: FeatureStream
    """
    def _stream(ts,tf):
        X = [common.sentinelize(x, [sys.maxint]*len(x.fields)) for x in ts]
        S = [[(-sys.maxint,-sys.maxint,0.0)] for t in ts]
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
                score = 0.0
                scores_y = []
                for s in S[i]:
                    if yend <= s[0]:   continue
                    if s[0] <  ystart: start = ystart
                    else:              start = s[0]
                    if yend <  s[1]:   end   = yend
                    else:              end   = s[1]
                    scores_y.extend([s[2]]*(end-start))
                if fn == 'mean': score = sum(scores_y)/(yend-ystart)
                elif fn == 'sum': score = sum(scores_y)
                elif fn == 'min': score = min(scores_y)
                elif fn == 'max': score = max(scores_y)
                elif fn == 'median':
                    scores_y.sort()
                    score = (scores_y[(len(scores_y)-1)/2] + scores_y[len(scores_y)/2]) / 2.
                elif hasattr(fn,'__call__'): score = fn(scores_y)
                scores += (score,)
            yield tuple(y)+scores

    if not(isinstance(trackScores,(list,tuple))): trackScores = [trackScores]
    if isinstance(trackFeatures,(list,tuple)):
        trackFeatures = concatenate(trackFeatures)
    if len(trackScores)>1:
        _fields = ["score"+str(i) for i in range(len(trackScores))]
    else:
        _fields = ["score"]
    _ts = [common.reorder(t,['start','end','score']) for t in trackScores]
    return track.FeatureStream(_stream(_ts,trackFeatures), trackFeatures.fields+_fields)

###############################################################################
def window_smoothing( trackList, window_size, step_size=1, stop_val=sys.maxint,
                      featurewise=False ):
    """
    Given a (list of) signal track(s) *trackList*, a *window_size* (in base pairs by default,
    or in number of feature if *featurewise* is True),  and a *window_step*,
    returns new signal tracks with, at each position p (multiples of *step_size*),
    the average score in the window [p-L, p+L]::

        X: __________666666666666____________
        R: ______12345666666666654321________ (not exact scores here)

    :param trackList: FeatureStream, or list of FeatureStream objects.
    :param window_size: (int) window size in bp.
    :param step_size: (int) step length. [1]
    :param stop_val: (int) ? . [sys.maxint]
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
                if fstart>win_start: steps += [fstart-win_start]
                else: delta -= F[0][2]
                if lstart>win_end:   steps += [lstart-win_end]
                else: delta += F[-1][2]
                nsteps = min(steps)
                sst = -(win_start % step_size)%step_size
                sen = -(win_start+nsteps % step_size)%step_size
                win_center = (win_start+win_end)/2
                if abs(delta) > 1e-11:
                    delta *= denom
                    score += delta*sst
                    for step in xrange(sst,nsteps,step_size):
                        if score>0 and win_center+step>=0 and win_center+step+step_size<=stop_val:
                            yield (win_center+step,win_center+step+step_size,score)
                        score += delta*step_size
                    score -= delta*sen
                else:
                    if score>0 and win_center+sst>=0 and win_center+sen+nsteps<=stop_val:
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
                if fstart>win_start: steps += [fstart-win_start]
                else: delta -= F[0][2]
                nsteps = min(steps)
                sst = -(win_start % step_size)%step_size
                sen = -(win_start+nsteps % step_size)%step_size
                win_center = (win_start+win_end)/2
                if abs(delta) > 1e-11:
                    delta *= denom
                    score += delta*sst
                    for step in xrange(sst,nsteps,step_size):
                        if score>0 and win_center+step>=0 and win_center+step+step_size<=stop_val:
                            yield (win_center+step,win_center+step+step_size,score)
                        score += delta*step_size
                    score -= delta*sen
                else:
                    if score>0 and win_center+sst>=0 and win_center+sen+nsteps<=stop_val:
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
        return [track.FeatureStream(call(common.reorder(t,_f),win_start,denom),fields=_f)
                for n,t in enumerate(trackList)]
    else:
        return track.FeatureStream(call(common.reorder(trackList,_f),win_start,denom),fields=_f)
