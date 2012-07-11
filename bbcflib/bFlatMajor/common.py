from bbcflib.btrack import FeatureStream
import sys

####################################################################
def sentinelize(iterable, sentinel):
    """Append *sentinel* at the end of *iterable* (avoid StopIteration error)."""
    for item in iterable: yield item
    yield sentinel

####################################################################
def reorder(stream,fields):
    """Reorders *stream.fields* so that *fields* come first.

    :param stream: FeatureStream object.
    :param fields: list of field names.
    :rtype: FeatureStream
    """
    if not(hasattr(stream, 'fields')) or stream.fields is None:
        return stream
    if not(all([f in stream.fields for f in fields])):
        raise ValueError("Need %s fields in stream."%(", ".join(fields)))
    if all(stream.fields[n] == f for n,f in enumerate(fields)):
        return stream
    _inds = [stream.fields.index(f) for f in fields]+[n for n,f in enumerate(stream.fields) if f not in fields]
    _flds = [stream.fields[n] for n in _inds]
    return FeatureStream((tuple(x[n] for n in _inds) for x in stream), fields=_flds)

####################################################################
def unroll( stream, start, end, fields=['score'] ):
    """Creates a stream of *end*-*start* items with appropriate *fields* values at every base position.
    For example, ``unroll([(10,12,0.5), (14,15,1.2)], start=9, end=16)`` returns::

        FeatureStream([(0,),(0.5,),(0.5,),(0,),(0,),(1.2,),(0,)])
                        9     10     11    12   13    14    15

    :param stream: FeatureStream object.
    :param start, end: (int) bounds of the region to return.
    :param fields: list of field names **in addition to 'start','end'**. [['score']]
    :rtype: FeatureStream
    """
    if not(isinstance(fields,(list,tuple))): fields = [fields]
    s = reorder(stream,['start','end']+fields)
    def _unr(s):
        pos = start
        for x in s:
            if x[1] <= pos: next
            while pos < min(x[0],end):
                yield (0,)+x[3:]
                pos+=1
            while pos < min(x[1],end):
                yield x[2:]
                pos+=1
            if pos >= end: break
        while pos < end:
            yield (0,)+x[3:]
            pos+=1
    return FeatureStream(_unr(s),fields=s.fields[2:])

####################################################################
def sorted_stream(stream,chrnames=[],fields=['chr','start','end']):
    """Sorts a stream according to *fields* values. Will load the entire stream in memory.
    The order of names in *chrnames* is used to to sort the 'chr' field if available.

    :param stream: FeatureStream object.
    :param chrnames: list of chrmosome names.
    :param fields: list of field names. [['chr','start','end']]
    :rtype: FeatureStream
    """
    s = reorder(stream,fields)
    sort_list = []
    feature_list = []
    for n,f in enumerate(s):
        if f[0] in chrnames: fi1 = chrnames.index(f[0])
        else: fi1 = f[0]
        sort_list.append((fi1,f[1],f[2],n))
        feature_list.append(f)
    sort_list.sort()
    def _sorted_stream(l1,l2):
        for t in l1:
            yield l2[t[-1]]
    return FeatureStream(_sorted_stream(sort_list,feature_list), stream.fields)

####################################################################
def shuffled(stream, chrlen=sys.maxint, repeat_number=1, sorted=True):
    """Yields randomly located features of the same length as the original stream.

    :param stream: FeatureStream object.
    :param chrlen: (int) chromosome length.
    :param repeat_number: (int) ??? [1]
    :rtype: FeatureStream
    """
    import random
    _f = ['start','end']
    features = reorder(stream,_f)
    def _shuffled(_s):
        randpos = []
        for feat in _s:
            feat_len = feat[1]-feat[0]
            for s in xrange(repeat_number):
                if len(randpos) == 0:
                    randpos = [random.randint(0, chrlen-feat_len) for i in xrange(10000)]
                start = randpos.pop()
                yield (start,start+feat_len)+feat[2:]
    if sorted:
        return sorted_stream(FeatureStream(_shuffled(features),features.fields),fields=_f)
    else:
        return FeatureStream(_shuffled(features),features.fields)

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
        x0 = x[0]
        for y in x[1:]:
            x0 += tuple(y)
        return x0

aggreg_functions = {'strand': strand_merge, 'chr': no_merge}

def fusion(stream,aggregate=aggreg_functions):
    """Fuses overlapping features in *stream* and applies *aggregate[f]* function to each field $f$.

    Example::

        [('chr1',10,15,'A',1),('chr1',13,18,'B',-1),('chr1',18,25,'C',-1)]

        yields

        (10, 18, 'chr1', 'A|B', 0)
        (18, 25, 'chr1', 'C', -1)

    :param stream: FeatureStream object.
    :rtype: FeatureStream
    """
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
    return FeatureStream( _fuse(_s), _s.fields )

def cobble(stream,aggregate=aggreg_functions):
    """Fuses overlapping features in *stream* and applies `aggregate[f]` function to each field `f`.

    Example::

        [('chr1',10,15,'A',1),('chr1',13,18,'B',-1),('chr1',18,25,'C',-1)]

        yields

        (10, 13, 'chr1', 'A', 1)
        (13, 15, 'chr1', 'AB', 0)
        (15, 18, 'chr1', 'B', -1)
        (18, 25, 'chr1', 'C', -1)

    This is to avoid having overlapping coordinates of features from both DNA strands,
    which some genome browsers cannot handle for quantitative tracks.

    :param stream: FeatureStream object.
    :rtype: FeatureStream
    """
    def _intersect(A,B):
        """Return *z*, the part that must replace A in *toyield*, and
        *rest*, that must reenter the loop instead of B."""
        rest = None
        if B[1] < A[2]:           # has an intersection
            if B[2] < A[2]:
                if B[1] == A[1]:  # same left border, A is bigger
                    z = [(A[0],B[1],B[2],A[3]+B[3]), (A[0],B[2],A[2],A[3])]
                elif B[1] < A[1]: # B shifted to the left wrt A
                    z = [(A[0],B[1],A[1],B[3]), (A[0],A[1],B[2],A[3]+B[3]), (A[0],B[2],A[2],A[3])]
                else:             # B embedded in A
                    z = [(A[0],A[1],B[1],A[3]), (A[0],B[1],B[2],A[3]+B[3]), (A[0],B[2],A[2],A[3])]
            elif B[2] == A[2]:
                if B[1] == A[1]:  # identical
                    z = [(A[0],A[1],A[2],A[3]+B[3])]
                elif B[1] < A[1]: # same right border, B is bigger
                    z = [(A[0],B[1],A[1],B[3]), (A[0],A[1],B[2],A[3]+B[3])]
                else:             # same right border, A is bigger
                    z = [(A[0],A[1],B[1],A[3]), (A[0],B[1],B[2],A[3]+B[3])]
            else:
                if B[1] == A[1]:  # same left border, B is bigger
                    z = [(A[0],A[1],A[2],A[3]+B[3])]
                    rest = (A[0],A[2],B[2],B[3])
                elif B[1] < A[1]: # B contains A
                    z = [(A[0],B[1],A[1],B[3]), (A[0],A[1],A[2],A[3]+B[3])]
                    rest = (A[0],A[2],B[2],B[3])
                else:             # B shifted to the right wrt A
                    z = [(A[0],A[1],B[1],A[3]), (A[0],B[1],A[2],A[3]+B[3])]
                    rest = (A[0],A[2],B[2],B[3])
        else: z = None            # no intersection
        return z, rest

    def _fuse(stream):
        x = stream.next()
        toyield = [x]
        while 1:
            try:
                x = stream.next()
                intersected = False
                for y in toyield:
                    replace, rest = _intersect(y,x)
                    if replace:
                        intersected = True
                        iy = toyield.index(y)
                        y = toyield.pop(iy)
                        for k in range(len(replace)):
                            toyield.insert(iy+k, replace[k])
                        x = rest
                        if not x: break
                if not intersected:
                    while toyield:
                        y = toyield.pop(0)
                        yield tuple(y)
                    toyield = [x]
                elif x:
                    toyield.append(x)
            except StopIteration:
                while toyield:
                    y = toyield.pop(0)
                    yield tuple(y)
                break

    _s = reorder(stream,['start','end'])
    return FeatureStream( _fuse(_s), _s.fields )

####################################################################
