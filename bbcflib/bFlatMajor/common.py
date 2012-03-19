from bbcflib import btrack as track
####################################################################
def sentinelize(iterable, sentinel):
    '''
    Add an item to the end of an iterable

    >>> list(sentinelize(range(4), 99))
    [0, 1, 2, 3, 99]
    '''
    for item in iterable: yield item
    yield sentinel

####################################################################
def reorder(stream,fields):
    if not(all([f in stream.fields for f in fields])):
        raise ValueError("Need %s fields in stream."%(", ".join(fields)))
    if all(stream.fields[n] == f for n,f in enumerate(fields)):
        return stream
    _inds = [stream.fields.index(f) for f in fields]+[n for n,f in enumerate(stream.fields) if f not in fields]
    _flds = [stream.fields[n] for n in _inds]
    return track.FeatureStream((tuple(x[n] for n in _inds) for x in stream), fields=_flds)

####################################################################
def unroll( stream, start, end, fields=['score'] ):
    if not(isinstance(fields,(list,tuple))): fields = [fields]
    s = reorder(stream,['start','end']+fields)
    def _unr(s):
        pos = start
        for x in s:
            if x[1]<=pos: next
            while pos<min(x[0],end):
                yield (0,)+x[3:]
                pos+=1
            while pos<min(x[1],end):
                yield x[2:]
                pos+=1
            if pos>=end: break
    return track.FeatureStream(_unr(s),fields=s.fields[2:])
####################################################################
def strand_merge(x): 
    return all(x[0]==y for y in x[1:]) and x[0] or 0

def no_merge(x): 
    return x[0]

def generic_merge(x):
    if isinstance(x[0],(int, long, float, complex)):
        return sum(x)
    if isinstance(x[0],basestring):
        return "|".join(x)
    if isinstance(x[0],tuple):
        for y in x[1:]:
            x[0] += tuple(y)
        return x[0]

aggreg_functions = {'strand': strand_merge, 'chr': no_merge}

def fusion(stream,aggregate=aggreg_functions):
    def _fuse(s):
        try:
            x = list(s.next())
        except StopIteration:
            return
        has_chr = 'chr' in s.fields
        if has_chr: chridx = s.fields.index('chr')
        for y in s:
            new_chr = has_chr and (x[chridx] != y[chridx])
            if y[0] < x[1] and not(new_chr):
                x[1] = max(x[1], y[1])
                x[2:] = [aggregate.get(f,generic_merge)((x[n+2],y[n+2]))
                         for n,f in enumerate(s.fields[2:])]
            else:
                yield tuple(x)
                x = list(y)
        yield tuple(x)
    _s = reorder(stream,['start','end'])
    return track.FeatureStream( _fuse(_s), _s.fields )

####################################################################
