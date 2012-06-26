"""
Examples::

    import btrack as track

    track.convert("data/test.bed","test0.sql",chrmeta='mm9')

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

    from bFlatMajor.common import shuffled
    for n,x in enumerate(shuffled(track_in.read('chr2'),chrlen=chrmeta['chr2']['length'])):
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

def track( path, format=None, **kwargs):
    """
    Guess file format.
    Return a Track object of the corresponding subclass (e.g. BedTrack).

    :param path: (str) path to a track-like file.
    :param format: (str) format of the file.
    :param **kwargs: (dict) parameters of the Track subclass' constructor.
    """
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

def convert( source, target, chrmeta=None, info=None ):
    """
    Converts a file from one format to another. Format can be explicitly specified::

        convert(('file1','bed'), ('file2','sql')) ,

    otherwise it is guessed first from file extension::

        convert('file1.bed', 'file2.sql')

    or in the worst case, by reading the first lines of the file.

    :param source: (str or tuple) path to the source file, or tuple of the form (path, format).
    :param target: (str or tuple) path to the target file, or tuple of the form (path, format).
    """
    if isinstance(source, tuple):
        tsrc = track(source[0], format=source[1], chrmeta=chrmeta)
    else:
        tsrc = track(source, chrmeta=chrmeta)
    if isinstance(target, tuple):
        ttrg = track(target[0], format=target[1], chrmeta=tsrc.chrmeta, fields=tsrc.fields, info=info)
    else:
        ttrg = track(target, chrmeta=tsrc.chrmeta, fields=tsrc.fields, info=info)
    ttrg.write( tsrc.read() )
    ttrg.close()
    tsrc.close()

def concat_fields( stream, infields, outfield='name', separator='|', as_tuple=False ):
    """
    Concatenate fields of a stream. Ex.:

    ('chr1', 12, 'aa', 'bb') -> ('chr1', 12, 'aa|bb')     # as_tuple=False
    ('chr1', 12, 'aa', 'bb') -> ('chr1', 12, ('aa','bb')) # as_tuple=True

    :param stream: FeatureStream object.
    :param infields: (list of str) list of fields to concatenate.
    :param outfield: (str) name of the new field created by concatenation of *infields*
        (can be an already existing one).
    :param separator: (str) char to add between entries from concatenated fields.
    :param as_tuple: (bool) join concatenated field entries in a tuple instead of a
        separator in a single string.
    :rtype: FeatureStream object.
    """
    _infields = [f for f in stream.fields if not(f in infields)] # untouched fields
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
    """
    Translate the chromosome identifiers in *stream* into chromosome names of the type 'chr5'.

    :param stream: FeatureStream object.
    :param assembly: genrep.Assembly object.
    :param keep: (bool) keep all features (True) or only those which chromosome identifier
        is recognized (False) [False].
    """
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
    """
    Split one field of a stream containing multiple information, into multiple fields. Ex.:

    ('chr1', 12, 'aa;bb;cc') -> ('chr1', 12, 'aa', 'bb', 'cc')

    :param stream: FeatureStream object.
    :param outfields: (list of str) list of new fields to be created.
    :param infield: (str) name of the field to be splitted.
    :param separator: (str) char separating the information in *infield*'s entries [';'].
    :param header_split: ?
    :param strip_input: (bool) ?
    """
    _outfields = stream.fields+[f for f in outfields if not(f in stream.fields)]
    in_indx = stream.fields.index(infield)
    in_out = [_outfields.index(f) for f in outfields]
    more_len = len(_outfields)-len(stream.fields)
    def _split(stream):
        for x in stream:
            y = list(x)+[None for f in range(more_len)]
            xsplit = x[in_indx].split(separator)
            if header_split:
                xmore = dict([re.search(r'\s*(\S+)'+header_split+'(\S*)',v+header_split).groups()
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
    """
    Filter the features of a track which score is above or below a certain threshold.

    :param source: Track instance (or a subclass), or a list/tuple of them.
    :param threshold: (float) threshold above which features are not retained (?)
    :param lower: (bool) higher (False) or lower (True) bound.
    :param fields: (str or list of str) names of the fields to apply the filter to.
    """
    if not(isinstance(fields,(list,tuple))):
           fields = [fields]
    if lower:
        selection = dict((f,(sys.float_info.min,threshold)) for f in fields)
    else:
        selection = dict((f,(threshold,sys.float_info.max)) for f in fields)
    #tsrc = track(source,fields=['chr','start','end','score'])
    if isinstance(source,(list,tuple)):
        return [t.read(selection=selection) for t in source]
    else:
        return source.read(selection=selection)

def strand_to_int(strand=None):
    """Convert +/- into 1/-1 notation for DNA strands."""
    if strand == '+': return 1
    if strand == '-': return -1
    return 0

def int_to_strand(num=0):
    """Convert 1/-1 into +/- notation for DNA strands."""
    num = int(num)
    if num > 0: return '+'
    if num < 0: return '-'
    return '.'

def format_float(f=float()):
    """Limit printing of a float to 4 decimals after the comma. :rtype: (str)"""
    return '%.4g' % float(f)

def format_int(i=int()):
    """:rtype: (str)"""
    return '%i' % int(i)

def ucsc_to_ensembl(start):
    """Add +1 to start coordinates going from UCSC to Ensembl annotation."""
    return format_int(int(start)+1)

def ensembl_to_ucsc(start):
    """Substract 1 to start coordinates going from Ensembl to UCSC annotation."""
    return format_int(int(start)-1)

################################################################################

class Track(object):
    """
    Metaclass regrouping the track properties.
    """
    def __init__(self, path, **kwargs):
        self.path = path
        self.filehandle = None
        self.format = kwargs.get("format")
        self.fields = kwargs.get("fields",[])
        self.types = {}
        self.assembly = kwargs.get('assembly')
        self.chrmeta = self._get_chrmeta(kwargs.get('chrmeta'))
        self.info = self._get_info(info=kwargs.get('info'))

    def _get_chrmeta(self,chrmeta=None):
        if isinstance(chrmeta,dict):
            return chrmeta
        if isinstance(chrmeta,basestring):
            self.assembly = chrmeta
        if self.assembly is None:
            return {}
        from bbcflib import genrep
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

    def save(self):
        self.open()
        self.close()

    def read(self, **kw):
        pass

    def write(self, **kw):
        pass

################################################################################
class FeatureStream(object):
    """
    Contains an iterator yielding features, and an extra fields attribute.

    :param data: an iterator, or cursor.
    :param fields: (list of str) the list of field names.
    """

    def __init__(self, data, fields=None):
        self.data = data
        if not fields and hasattr(data, 'description'):
            fields = [x[0] for x in data.description]
        self.fields = fields

    def __iter__(self): return self.data

    def next(self): return self.data.next()

################################################################################
