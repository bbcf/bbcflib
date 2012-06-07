import sys
from bbcflib.bFlatMajor import common
from bbcflib import btrack as track

def concatenate(trackList, fields=None):
###################################
    if len(trackList) == 1: return trackList[0]
    tracks = [track.FeatureStream(common.sentinelize(x,[sys.maxint]*len(x.fields)),x.fields)
              for x in trackList]

    def _find_min(ftple):
        nmin = 0
        xmin = ftple[0]
        for n,x in enumerate(ftple[1:]):
            for k in range(len(x)):
                if cmp(x[k],xmin[k])<0:
                    xmin = x
                    nmin = n+1
                    break
                if cmp(x[k],xmin[k])>0:
                    break
        return nmin

    def _knead(_t,N):
        current = [x.next()[:N] for x in _t]
        while (1):
            n = _find_min(current)
            if current[n][0] == sys.maxint: return
            yield current[n]
            current[n] = _t[n].next()[:N]
    
    if fields is None:
        fields = [f for f in tracks[0].fields if all(f in t.fields for t in tracks[1:])]
    _of = ['start','end']
    if 'chr' in fields: _of = ['chr']+_of
    if 'name' in fields: _of += ['name']
    _of += [f for f in fields if not(f in _of)]
    tl = [common.reorder(t,_of) for t in trackList]
    return track.FeatureStream(_knead(tl,len(fields)),fields=fields)

###############################################################################
def neighborhood(trackList, before_start=None, after_end=None, 
                 after_start=None, before_end=None, on_strand=False):
###################################
    def _generate_single(track,a,b,c,d):
        if a:
            for x in track:
                if on_strand and x[2]<0:
                    yield (x[0]-after_end,     x[0]+before_end+1) + x[2:]
                    yield (x[1]-after_start-1, x[1]+before_start) + x[2:]
                else:
                    yield (x[0]-before_start, x[0]+after_start+1) + x[2:]
                    yield (x[1]-before_end-1, x[1]+after_end)     + x[2:]
        elif b:
            for x in track:
                if on_strand and x[2]<0:
                    yield (x[1]-after_start-1, x[1]+before_start)  + x[2:]
                else:
                    yield (x[0]-before_start,  x[0]+after_start+1) + x[2:]
        elif c:
            for x in track:
                if on_strand and x[2]<0:
                    yield (x[0]-after_end,    x[0]+before_end+1) + x[2:]
                else:
                    yield (x[1]-before_end-1, x[1]+after_end)    + x[2:]
        elif d:
            for x in track:
                if on_strand and x[2]<0:
                    yield (x[0]-after_end,    x[1]+before_start) + x[2:]
                else:
                    yield (x[0]-before_start, x[1]+after_end)    + x[2:]

###################################
    _fields = ['start','end']
    if on_strand: _fields += ['strand']
    case1 = True
    case2 = True
    case3 = True
    case4 = True
    if before_start is None:
        case1 = case2 = case4 = False
    if after_end is None:
        case1 = case3 = case4 = False
    if after_start is None:
        case1 = case2 = False
    if before_end is None:
        case1 = case3 = False
    if isinstance(trackList,(list,tuple)):
        tl = [common.reorder(t,_fields) for t in trackList]
        return [track.FeatureStream(_generate_single(t,case1,case2,case3,case4),
                                    fields=t.fields) for t in tl]
    else: 
        tl = common.reorder(trackList,_fields)
        return track.FeatureStream(_generate_single(tl,case1,case2,case3,case4),
                                   fields=tl.fields)
###############################################################################
def _combine(tracks,N,fn,win_size,aggregate):
    fields = tracks[0].fields
    tracks = [common.sentinelize(t, [sys.maxint]*len(fields)) for t in tracks]
    init = [tracks[i].next() for i in range(N)]
    for i in xrange(N-1,-1,-1):
        if init[i][0] == sys.maxint: 
            N-=1
            tracks.pop(i)
            init.pop(i)
    if N == 0: return 
    available_tracks = range(N-1,-1,-1)
    current = [(init[i][0],i)+init[i][2:] for i in range(N)]+[(init[i][1],i) for i in range(N)]
    current.sort()

    start = current[0][0]
    activity = [False]*N
    z = [None]*N
    while current[0][0] == start:
        i = current[0][1]
        activity[i] = True
        z[i] = current.pop(0)[2:]

    k=1
    while available_tracks:
        # load *win_size* bp in memory
        to_remove = []
        for i in available_tracks:
            a = [0,0]
            limit = k*win_size
            while a[1] < limit:
                a = tracks[i].next()
                if a[0] == sys.maxint:
                    to_remove.append(i)
                else:
                    current.append((a[0],i)+a[2:])
                    current.append((a[1],i))
        for i in to_remove:
            available_tracks.remove(i)
        current.sort()

        # calculate boolean values for start-next interval
        while current and current[0][0] < limit:
            next = current[0][0]
            if fn(activity):
                feat_aggreg = tuple(aggregate.get(f,common.generic_merge)(tuple(zz[n] for zz in z if zz)) 
                                    for n,f in enumerate(fields[2:]))
                yield (start,next)+feat_aggreg
            while current and current[0][0] == next:
                i = current[0][1]
                activity[i] = not(activity[i])
                zi = current.pop(0)[2:]
                z[i] = activity[i] and zi or None
            start = next
        k+=1


def combine(trackList, fn, win_size=1000, 
            aggregate={'strand': common.strand_merge, 'chr': common.no_merge}):
    N = len(trackList)
    fields = ['start','end']
    if N<2: return trackList
    tracks = [common.fusion(common.reorder(t,fields=fields)) for t in trackList]
    return common.fusion(track.FeatureStream(_combine(tracks,N,fn,win_size,aggregate),
                                             fields=tracks[0].fields))

def exclude(x,indexList):
    return any([y for n,y in enumerate(x) if not(n in indexList)]) and all([not(y) for n,y in enumerate(x) if n in indexList])

def require(x,indexList):
    return any([y for n,y in enumerate(x) if not(n in indexList)]) and all([y for n,y in enumerate(x) if n in indexList])

def disjunction(x,indexList):
    complementList = [n for n in range(len(x)) if not(n in indexList)]
    return exclude(x,indexList) or exclude(x,complementList)

def intersection(x):
    return all(x)

def union(x):
    return any(x)

###############################################################################
def segment_features(trackList,nbins=10,upstream=None,downstream=None):
    """
    The features plus their upstream and downstream flanks must not overlap.
    """
    def _split_feat(_t):
        starti = _t.fields.index('start')
        endi   = _t.fields.index('end')
        if 'strand' in _t.fields:
            strandi = _t.fields.index('strand')
        else:
            strandi = -1
        uprel = False
        downrel = False
        if upstream is None:
            updist = 0
            upstep = 1
        elif upstream[0] > 1:
            updist = upstream[0]
            upstep = updist/upstream[1]
        else:
            uprel = True
        if downstream is None:
            downdist = 0
            downstep = 1
        elif downstream[0] > 1:
            downdist = downstream[0]
            downstep = downdist/downstream[1]
        else:
            downrel = True
        for x in _t:
            xs = list(x)
            xlen = (xs[endi]-xs[starti])
            xstep = xlen/nbins
            if uprel:
                updist = int(.5+upstream[0]*xlen)
                upstep = updist/upstream[1]
            if downrel:
                downdist = int(.5+downstream[0]*xlen)
                downstep = downdist/downstream[1]
            if strandi >= 0 and xs[strandi] > 0:
                allsteps = range(x[starti]-updist,x[starti],upstep)
                allsteps += range(x[starti],x[endi],xstep)
                allsteps += range(x[endi],x[endi]+downdist,downstep)
                start = allsteps[0]
                for s,n in enumerate(allsteps[1:]):
                    xs[starti] = start
                    xs[endi] = s
                    start = s
                    yield tuple(xs)+(n,)
            else:
                allsteps = range(x[starti]-downdist,x[starti],downstep)
                allsteps += range(x[endi],x[starti],xstep)
                allsteps += range(x[endi],x[endi]+updist,upstep)
                start = allsteps[0]
                ntot = len(allstart)-1
                for s,n in enumerate(allsteps[1:]):
                    xs[starti] = start
                    xs[endi] = s
                    start = s
                    yield tuple(xs)+(ntot-n,)
    if isinstance(trackList, (list,tuple)):
        return [track.FeatureStream(_split_feat(t), fields=t.fields+['bin'])
                for t in trackList]
    else:
        return track.FeatureStream(_split_feat(trackList),
                                   fields=trackList.fields+['bin'])
