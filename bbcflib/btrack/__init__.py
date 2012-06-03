"""
Examples::

    import track

    track.convert("data/test.bed","test0.sql",assembly='mm9')

    chrmeta = {'chr2':{'length':4000000}}
    info = {'datatype':'features'}
    track_in = track.track("data/test.bed",chrmeta=chrmeta)
    track_out = track.track("test1.sql",fields=['start','end','name'],chrmeta=chrmeta,info=info)
    track_out.write(track_in.read())
    track_out.close()
    track_in.close()

    track_in = track.track("data/HoxD13_4C_FB.sql")
    track_out = track.track("test2.wig")
    selection = [{'chr':'chr1','start':(7568000,9607000)},{'chr':'chr2','end':(3907400,4302000)}]
    track_out.write(track_in.read(selection=selection),mode='overwrite')
    track_out.close()
    track_in.close()
    
    track_in = track.track("data/Gene_TxS_chr2.bed.gz",chrmeta='mm9',format='bed')
    for x in track_in.read():
        print x #('chr2', 3030497, 3032496, 'ENSMUST00000072955_txS')
        break

    for x in track.split_field(track_in.read(),['name','extension'],'name','_'):
        print x #['chr2', 3030497, 3032496, 'ENSMUST00000072955', 'txS']
        break

    for n,x in enumerate(bFlatMajor.common.sorted_stream(track_in.read_shuffled()),track_in.chrmeta.keys()):
        print x
        if n>10: break

    selection = {'chr':'chr2','start':(7540000,75650000)}
    track_features = track.track("data/Bricks_HoxD4_FB_05_chr2.bed")
    track_scores = track.track("data/HoxD13_4C_FB.sql",readonly=True)
    score = float()
    length = float()
    for x in track_scores.read(selection=track_features.read(selection=selection),
                               fields=['start','end','score']):
        score += x[2]*(x[1]-x[0])
        length += (x[1]-x[0])
    track_features.close()
    track_scores.close()
    print score/length

"""

__all__ = ['Track','track','FeatureStream','strand_to_int','int_to_strand','format_float','format_int',
           'ucsc_to_ensembl','ensembl_to_ucsc',]

import sys, os, re

_track_map = {
    'sql': ('bbcflib.btrack.sql','SqlTrack'),
    'db':  ('bbcflib.btrack.sql','SqlTrack'),
    'text':('bbcflib.btrack.text','TextTrack'),
    'txt': ('bbcflib.btrack.text','TextTrack'),
    'bed': ('bbcflib.btrack.text','BedTrack'),
    'bedGraph': ('bbcflib.btrack.text','BedGraphTrack'),
    'bedgraph': ('bbcflib.btrack.text','BedGraphTrack'),
    'wig': ('bbcflib.btrack.text','WigTrack'),
    'gff': ('bbcflib.btrack.text','GffTrack'),
    'gtf': ('bbcflib.btrack.text','GffTrack'),
    'sga': ('bbcflib.btrack.text','SgaTrack'),
    'bigWig': ('bbcflib.btrack.bin','BigWigTrack'),
    'bigwig': ('bbcflib.btrack.bin','BigWigTrack'), 
    'bw':  ('bbcflib.btrack.bin','BigWigTrack'),
    'bam': ('bbcflib.btrack.bin','BamTrack'),
}

def track(path, format=None, **kwargs):
    if format is None:
        format = os.path.splitext(path)[1][1:]
        if format in ['gz','gzip']:
            format = os.path.splitext(path.strip("."+format))[1][1:]
        if format == '': 
            with open(path, 'r') as file:
                rstart = file.read(15)
                if rstart == "SQLite format 3": format='sql'
                else:
                    while rstart.startswith("#"):
                        rstart = file.readline()
                        rstart = file.read(1)
                    rstart += file.readline()
                    head = re.search(r'track\s+type=(\S+)',rstart)
                    if head: format = head.groups()[0]
    if not(format in _track_map): 
        raise Exception("The format '%s' is not supported."%format)
    __import__(_track_map[format][0])
    return getattr(sys.modules[_track_map[format][0]],
                   _track_map[format][1])(path,**kwargs)

def convert( source, target, assembly=None, chrmeta=None ):
    if isinstance(source, tuple):
        tsrc = track(source[0], format=source[1], assembly=assembly, chrmeta=chrmeta)
    else:
        tsrc = track(source, assembly=assembly, chrmeta=chrmeta)
    if isinstance(target, tuple):
        ttrg = track(target[0], format=target[1], assembly=assembly, chrmeta=tsrc.chrmeta, fields=tsrc.fields)
    else:
        ttrg = track(target, assembly=assembly, chrmeta=tsrc.chrmeta, fields=tsrc.fields)
    ttrg.write( tsrc.read() )
    ttrg.close()
    tsrc.close()

def concat_fields( stream, infields, outfield='name', separator='|', as_tuple=False ):
    _infields = [f for f in stream.fields if not(f in infields)]
    in_out_indx = [stream.fields.index(f) for f in _infields]
    to_extend = []
    if not(outfield in _infields): 
        _infields += [outfield]
        to_extend = [None]
    out_indx = _infields.index(outfield)
    in_indx = [stream.fields.index(f) for f in infields]
    def _concat(stream):
        for x in stream:
            y = [x[i] for i in in_out_indx]+to_extend
            if as_tuple:
                y[out_indx] = tuple((x[i] for i in in_indx))
            else:
                y[out_indx] = separator.join([str(x[i]) for i in in_indx])
            yield tuple(y)
    return FeatureStream(_concat(stream),_infields)

def map_chromosomes( stream, assembly, keep=False ):
    if not('chr' in stream.fields): return stream 
    ic = stream.fields.index('chr')
    chrom_map = {}
    for k,c in assembly.chromosomes.iteritems():
        cname = c['name']
        chrom_map[cname] = cname
        if cname.startswith('chr') and len(cname)>3: chrom_map[cname[3:]] = cname
        chrom_map[k[0]] = cname
        chrom_map[str(k[1])+"."+str(k[2])] = cname
        if c['synonyms']: 
            for s in c['synonyms'].split(','): chrom_map[s] = cname
    if keep:
        return FeatureStream((x[:ic]+(chrom_map.get(x[ic],x[0]),)+x[ic+1:] 
                              for x in stream),stream.fields)
    else:
        return FeatureStream((x[:ic]+(chrom_map[x[ic]],)+x[ic+1:] 
                              for x in stream if x[ic] in chrom_map),stream.fields)

def split_field( stream, outfields, infield='name', separator=';',
                 header_split=None, strip_input=False ):
    _outfields = stream.fields+[f for f in outfields if not(f in stream.fields)]
    in_indx = stream.fields.index(infield)
    in_out = [_outfields.index(f) for f in outfields]
    more_len = len(_outfields)-len(stream.fields)
    def _split(stream):
        for x in stream:
            y = list(x)+[None for f in range(more_len)]
            xsplit = x[in_indx].split(separator)
            if header_split:
                if not(header_split in v): continue
                xmore = dict([re.search(r'\s*(\S+)'+header_split+'(\S+)',v).groups()
                              for v in xsplit])
                for n,f in enumerate(outfields):
                    y[in_out[n]] = xmore.get(f,'').strip('"')
                if strip_input:
                    y[in_indx] = separator.join([str(k)+header_split+str(v)
                                                 for k,v in xmore.iteritems()
                                                 if not(k in outfields)])
            else:
                for n,v in enumerate(xsplit):
                    if n >= len(in_out): break
                    y[in_out[n]] = v
                if strip_input:
                    y[in_indx] = separator.join(xsplit[n:])
            yield tuple(y)
    return FeatureStream(_split(stream),_outfields)

def score_threshold( source, threshold=0.0, lower=False, fields='score' ):
    if not(isinstance(fields,(list,tuple))):
           fields = [fields]
    if lower:
        selection = dict((f,(sys.float_info.min,threshold)) for f in fields)
    else:
        selection = dict((f,(threshold,sys.float_info.max)) for f in fields)
    tsrc = track(source,fields=['chr','start','end','score'])
    if isinstance(source,(list,tuple)):
        return [t.read(selection=selection) for t in source]
    else:
        return source.read(selection=selection)

def strand_to_int(strand=None):
    if strand == '+': return 1
    if strand == '-': return -1
    return 0

def int_to_strand(num=0):
    num = int(num)
    if num > 0: return '+'
    if num < 0: return '-'
    return '.'

def format_float(f=float()):
    return '%.4g' % float(f)

def format_int(i=int()):
    return '%i' % int(i)

def ucsc_to_ensembl(start):
    return format_int(int(start)+1)

def ensembl_to_ucsc(start):
    return format_int(int(start)-1)

################################################################################

class Track(object):
    def __init__(self, path,**kwargs):
        self.path = path
        self.filehandle = None
        self.format = None
        self.fields = []
        self.types = {}
        self.assembly = kwargs.get('assembly')
        self.chrmeta = self._get_chrmeta(kwargs.get('chrmeta'))
        self.info = self._get_info(info=kwargs.get('info'))

    def _get_chrmeta(self,chrmeta=None):
        from bbcflib import genrep
        if isinstance(chrmeta,dict):
            return chrmeta
        if isinstance(chrmeta,basestring):
            self.assembly = chrmeta
        if self.assembly is None:
            return {}
        return genrep.Assembly(self.assembly).chrmeta

    def _get_info(self,info=None):
        pass

    def __enter__(self):
        """Called when evaluating the 'with' statement."""
        return self

    def __exit__(self, errtype, value, traceback):
        self.close()

    def __iter__(self):
        """Iterates on the list of features."""
        return iter(self.read())

    def open(self):
        pass

    def close(self):
        pass

    def read(self, **kw):
        pass

    def write(self, **kw):
        pass

    def read_shuffled(self, repeat_number=1, **kwargs):
        ''' Yields randomly located features of the same length as the original track.'''
        import random
        feat_iter = self.read(**kwargs)
        try:
            chri = feat_iter.fields.index('chr')
            starti = feat_iter.fields.index('start')
            endi = feat_iter.fields.index('end')
        except:
            raise Error('Need chr, start and end fields to shuffle.')
        def _shuffled(iter):
            chr_cur = None
            for feat in iter:
                chr_len = self.chrmeta[feat[chri]]['length']
                if chr_cur != feat[chri]:
                    randpos = []
                    chr_cur = feat[chri]
                feat_len = feat[endi]-feat[starti]
                featnew = list(feat)
                for s in range(repeat_number):
                    if len(randpos) == 0:
                        randpos = [random.randint(0, chr_len-feat_len)
                                   for i in range(10000)]
                    featnew[starti] = randpos.pop()
                    featnew[endi] = featnew[starti]+feat_len
                    yield tuple(featnew)
        return FeatureStream(_shuffled(feat_iter),feat_iter.fields)

    

################################################################################
class FeatureStream(object):
    """Contains an iterator yielding features and an extra
       fields attribute.

       @param data: the iterator (or cursor) itself.
       @param fields: the list of fields
    """

    def __init__(self, data, fields=None):
        self.data = data
        if not fields and hasattr(data, 'description'): 
            fields = [x[0] for x in data.description]
        self.fields = fields

    def __iter__(self): return self.data

    def next(self): return self.data.next()

################################################################################
