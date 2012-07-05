import sys
from bbcflib.bFlatMajor import common
from bbcflib import btrack as track

def merge_scores(trackList, geometric=False):
    """
    Average several score tracks.

    X1: ▁▁▁▁▁▁▁▁▁▁█████████▁▁▁▁▁▁
    X2: ▁▁▁▁▁▅▅▅▅▅▅▅▅▅▅▁▁▁▁▁▁▁▁▁▁
    R:  ▁▁▁▁▁▂▂▂▂▂▇▇▇▇▇▅▅▅▅▁▁▁▁▁▁

    :param trackList: list of FeatureStream objects.
    :param geometric: (bool) set True to use the geometric mean instead of arithmetic.
    """
    tracks = [track.FeatureStream(common.sentinelize(x,[sys.maxint]*len(x.fields)), x.fields)
              for x in trackList]
    _fields = ['start','end','score']
    if all(['name' in t.fields for t in trackList]): _fields.append('name')
    tracks = [common.reorder(t,_fields) for t in tracks]
    elements = [list(x.next()) for x in tracks]
    track_denom = 1.0/len(trackList)

    def arithmetic_mean(scores,denom):
        return sum(scores)*denom
    def geometric_mean(scores,denom):
        return (reduce(lambda x, y: x*y, scores))**denom
## more precise/less efficient: exp(sum([log(x) for x in scores])*denom)
    mean_fn = (geometric and geometric_mean) or arithmetic_mean
    for i in xrange(len(tracks)-1, -1, -1):
        if elements[i][0] == sys.maxint:
            tracks.pop(i)
            elements.pop(i)
    def _stream(tracks):
        while tracks:
            start = min([x[0] for x in elements])
            end = min([x[0] for x in elements if x[0]>start]+[x[1] for x in elements])
            scores = [x[2] for x in elements if x[1]>start and x[0]<end]
            if 'name' in _fields:
                name = "|".join([x[3] for x in elements if not(x[3] is None) and x[1]>start and x[0]<end])
                yield (start, end, mean_fn(scores,track_denom), name)
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

    return track.FeatureStream(_stream(tracks),_fields)

###############################################################################
def mean_score_by_feature(trackScores,trackFeatures):
    """
    Given a score track X and a feature track Y, compute the mean of scores of
    every of Y's features in X. The output consists of a feature track
    similar to ``Y`` but with a new score value property for every feature.

    X: ──────▤▤▤▤▤▤▤▤▤▤──────────────▤▤▤▤▤▤▤▤▤▤──────
    Y: ▁▁▁▁▁▁▁▁▁▁▁█████████▁▁▁▁▁▁▁▁▁▁██████████▁▁▁▁▁▁
    R: ▁▁▁▁▁▁▅▅▅▅▅▅▅▅▅▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁██████████▁▁▁▁▁▁

    :param trackScores: (FeatureStream) score track.
    :param trackFeatures: (FeatureStream) feature track.
    """
    def _stream(ts,tf):
        X = [common.sentinelize(x, [sys.maxint]*len(x.fields)) for x in ts]
        F = [[(-sys.maxint,-sys.maxint,0.0)] for t in ts]
        for y in tf:
            ystart = y[tf.fields.index('start')]
            yend = y[tf.fields.index('end')]
            scores = ()
            for i in range(len(ts)):
                xnext = F[i][-1]
                while xnext[0] < yend:
                    xnext = X[i].next()
                    if xnext[1] > ystart: F[i].append(xnext)
                n = 0
                while F[i][n][1] <= ystart: n+=1
                F[i] = F[i][n:]
                score = 0.0
                for f in F[i]:
                    if yend <= f[0]:   continue
                    if f[0] <  ystart: start = ystart
                    else:              start = f[0]
                    if yend <  f[1]:   end   = yend
                    else:              end   = f[1]
                    score += (end-start)*f[2]
                scores += (score/(yend-ystart),)
            yield tuple(y)+scores
    if not(isinstance(trackScores,(list,tuple))): trackScores = [trackScores]

    if isinstance(trackFeatures,(list,tuple)): trackFeatures = trackFeatures[0]
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
    Given a signal track X and a window size in base pairs, return a new signal
    track with, at each position p, the mean of the scores in the window [p-L, p+L].
    Border cases are handled by zero padding and the signal's support is invariant.

    X: ▁▁▁▁▁▁▁▁▁▁████████████▁▁▁▁▁▁▁▁▁▁▁▁
    R: ▁▁▁▁▁▁▂▄▅▇████████████▇▅▄▂▁▁▁▁▁▁▁▁

    :param trackList: FeatureStream, or list of FeatureStream objects.
    :param window_size: (int) window size in bp.
    :param step_size: (int)  [1]
    :param stop_val: (int)  [sys.maxint]
    :param featurewise: (bool)  [False]

    Example of windows, window_size=9, step_size=3
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
