"""
Documentation `here <http://bbcf.epfl.ch/bbcflib/tutorial_track.html>`_.
"""

__all__ = ['Track','track','FeatureStream','convert',
           'strand_to_int','int_to_strand','format_float','format_int',
           'ucsc_to_ensembl','ensembl_to_ucsc']

import sys, os, re

_track_map = {
    'sql': ('bbcflib.track.sql','SqlTrack'),
    'db':  ('bbcflib.track.sql','SqlTrack'),
    'text':('bbcflib.track.text','TextTrack'),
    'txt': ('bbcflib.track.text','TextTrack'),
    'bed': ('bbcflib.track.text','BedTrack'),
    'bedGraph': ('bbcflib.track.text','BedGraphTrack'),
    'bedgraph': ('bbcflib.track.text','BedGraphTrack'),
    'sga': ('bbcflib.track.text','SgaTrack'),
    'wig': ('bbcflib.track.text','WigTrack'),
    'gff': ('bbcflib.track.text','GffTrack'),
    'gtf': ('bbcflib.track.text','GffTrack'),
    'bigWig': ('bbcflib.track.bin','BigWigTrack'),
    'bigwig': ('bbcflib.track.bin','BigWigTrack'),
    'bw':  ('bbcflib.track.bin','BigWigTrack'),
    'bam': ('bbcflib.track.bin','BamTrack'),
    'sam': ('bbcflib.track.text','SamTrack'),
    'fps': ('bbcflib.track.text','FpsTrack'),
}

def track( path, format=None, **kwargs):
    """
    Guess file format and return a Track object of the corresponding subclass (e.g. BedTrack).

    :param path: (str) name of/path to a track-like file. If the file does not exist yet,
        a new track-like file of the requested *format* will be created at this location
        on closure, if data is added to the track (using `write()`).
    :param format: (str) format of the file. If not provided, the format is set according
        to the file's extension.
    :param **kwargs: (dict) parameters of the Track subclass' constructor.
        Typically `assembly` or `chrmeta`.
    """
    assert isinstance(path, basestring), "Expected string or unicode, found %s." % type(path)
    if format is None:
        path2, format = os.path.splitext(path)
        format = format.lstrip('.')
        if format in ['gz','gzip']:
            format = os.path.splitext(path2)[1].lstrip('.')
        if format == '':
            if os.path.exists(path):
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
            else:
                raise ValueError("Format of file %s not known." %path)
    if not(format in _track_map):
        raise Exception("Format '%s' is not supported."%format)
    __import__(_track_map[format][0])
    return getattr(sys.modules[_track_map[format][0]],
                   _track_map[format][1])(path,**kwargs)

def convert( source, target, chrmeta=None, info=None, mode='write', clip=False ):
    """
    Converts a file from one format to another. Format can be explicitly specified::

        convert(('file1','bed'), ('file2','sql')) ,

    otherwise it is guessed first from file extension::

        convert('file1.bed', 'file2.sql')

    or in the worst case, by reading the first lines of the file.

    :param source: (str or tuple) path to the source file, or tuple of the form (path, format).
    :param target: (str or tuple) path to the target file, or tuple of the form (path, format).
    :param chrmeta: (dict) to specify manually 'chrmeta' for both input and output tracks. [None]
    :param info: (dict) info that will be available as an attribute of the output track. [None]
    :param mode: (str) writing mode: either 'write', 'append' or 'overwrite'. ['write']
    """
    if isinstance(source, tuple):
        tsrc = track(source[0], format=source[1], chrmeta=chrmeta)
    else:
        tsrc = track(source, chrmeta=chrmeta)
    _f = tsrc.fields
    if not('chr' in _f): _f = ['chr']+_f
    if isinstance(target, tuple):
        ttrg = track(target[0], format=target[1], chrmeta=tsrc.chrmeta, fields=_f, info=info)
    else:
        ttrg = track(target, chrmeta=tsrc.chrmeta, fields=_f, info=info)
    try:
        ttrg.write( tsrc.read(), mode=mode, clip=clip )
    finally:
        ttrg.close()
        tsrc.close()
    return ttrg

def strand_to_int(strand=''):
    """Convert +/- into 1/-1 notation for DNA strands."""
    return {'1':1, '-1':-1, '+':1, '-':-1, 'fwd':1, 'rev':-1}.get(str(strand),0)

def int_to_strand(num=0):
    """Convert 1/-1 into +/- notation for DNA strands."""
    return {'1':'+', '-1':'-', '+':'+', '-':'-'}.get(str(num),'.')

def format_float(f=float()):
    """Return a formatted string from a float or a string representing a float.
    Limit to 4 decimals after the comma."""
    return '%.4g' % float(f)

def format_int(i=int()):
    """Return a formatted string from an integer or a string representing an integer."""
    return '%i' % int(i)

def ucsc_to_ensembl(stream):
    """Shifts start coordinates 1 base to the right, to map UCSC to Ensembl annotation."""
    istart = stream.fields.index('start')
    return FeatureStream((item[:istart]+(item[istart]+1,)+item[istart+1:]
                          for item in stream),fields=stream.fields)

def ensembl_to_ucsc(stream):
    """Shifts start coordinates 1 base to the left, to map Ensembl to UCSC annotation."""
    istart = stream.fields.index('start')
    return FeatureStream((item[:istart]+(item[istart]-1,)+item[istart+1:]
                          for item in stream),fields=stream.fields)

def check(source, out=sys.stdout,
          check_sorted=True, check_duplicates=True, check_zerosize=True, **kwargs):
    """Read a track-like file *source* once to check that it is correctly formatted
    (according to ``*source*.format``), plus other optional filters.

    :param source: (str) name of the file.
    :param check_sorted: verify that the file is sorted: chromosomes are grouped,
        and regions of each chromosome are sorted w.r.t. 'start' and 'end' in ascending order.
    :param check_duplicates: verify that no line is repeated exactly (except empty lines).
    :param check_zerosize: verify that no region has length zero (start==end).
    :param **kwargs: ``track`` keyword arguments.
    """
    if isinstance(source, basestring):
        t = track(source, chrmeta=kwargs.get('chrmeta'), **kwargs)
        filename = source
    else:
        t = source
        filename = os.path.basename(t.path)
    is_chr = 'chr' in t.fields
    is_start = 'start' in t.fields
    is_end = 'end' in t.fields
    chr_idx = t.fields.index('chr') if is_chr else None
    start_idx = t.fields.index('start') if is_start else None
    end_idx = t.fields.index('end') if is_end else None
    visited_chr = []
    last_chr = None
    last_start = last_end = 0
    lastrow = ''
    n = 0
    s = t.read()
    while 1:
        n += 1
        try:
            row = s.next()
        except StopIteration:
            return True
        except Exception, e:
            out.write("Check format: line %s of %s is not compatible with format %s. \
                       \nException raised: %s" % (n,source,t.format,e))
            return False
        if check_duplicates and row == lastrow and len(row)!=0:
            out.write("Check duplicates: %s: duplicate at line %d. \n\n" % (filename,n))
            return False
        chr   = row[chr_idx] if is_chr else None
        start = row[start_idx] if is_start else 0
        end   = row[end_idx] if is_end else 0
        if check_zerosize and start == end:
            out.write("Check zero size: %s: empty region at line %d.\n\n" % (filename,n))
            return False
        if check_sorted:
            if chr != last_chr:
                if chr in visited_chr:
                    out.write("Check order: %s: error at line %d: \
                               \n\tChromosome %s appears twice in the file.\n\n" % (filename,n,chr))
                    return False
                visited_chr.append(chr)
                last_start = last_end = 0
            elif start < last_start:
                out.write("Check order: %s: error at line %d: \
                           \n\tStart position %d < %d.\n\n" % (filename,n,start,last_start))
                return False
            elif start == last_start and end < last_end:
                out.write("Check order: %s: error at line %d: \
                           \n\tEnd position %d < %d.\n\n" % (filename,n,end,last_end))
                return False
        last_chr = chr
        last_start = start
        last_end = end
        lastrow = row

def stats(source, out=sys.stdout, plot=True, wlimit=80, **kwargs):
    """Prints stats about the track. Draws a plot of the scores distribution (if any)
    directly to the console.

    It does not load the whole file in memory, but the distribution of its scores,
    which are expected to be limited in variety for count data (normalized or not).

    :param source: (str) name of the file. Can also be a Track instance.
    :param out: writable/file object (default: stdout), or a dict (will be updated).
    :param wlimit: max width of the distribution plot - console screen -, in number of chars. [80]
    :param **kwargs: ``track`` keyword arguments.
    """
    def median(vals,distr,nfeat):
        smedian = cumul = 0
        lastv = vals[0]
        for v in vals:
            cumul += distr[v]
            if cumul > 0.5*nfeat:
                smedian = v if len(vals)%2==1 else (v+lastv)/2.
                break
            lastv = v
        return smedian

    def stats_from_distr(distr,nfeat):
        if nfeat==0: return (None,)*6
        vals = sorted(distr.keys())
        total = float(sum(k*distr[k] for k in vals))
        smin = vals[0]; smax = vals[-1]
        smean = total/nfeat
        stdev = (sum((x-smean)**2 for x in vals)/nfeat)**(0.5)
        smedian = median(vals,distr,nfeat)
        return total,smin,smax,smean,stdev,smedian

    def console_distr_plot(distr,out,hlimit,wlimit,binw):
        assert isinstance(wlimit,int) and wlimit > 1, "wlimit must be an integer."
        vals = sorted(distr.keys())
        maxval = vals[-1]
        nbins = int(maxval//binw + 1)
        bvals = [0]*nbins
        bscores = [0]*nbins
        for b in range(min(hlimit,nbins)):
            bvals[b] = []
            for v in vals:
                if v < (b+1)*binw:
                    bvals[b].append(v)
                else:
                    vals = vals[len(bvals[b]):]
                    break
            bscores[b] = sum(distr[v] for v in bvals[b])
        wlimit = wlimit-8 # to leave place for a legend
        max_bscore = max(bscores)
        if max_bscore == 0:
            out.write("| (0)\n")
            return
        legends = ["%.1f"%(b*binw) for b in range(min(hlimit,nbins))]
        lmargin = max(len(L) for L in legends)
        legends = [L+" "*(lmargin-len(L)) for L in legends]
        for b in range(min(hlimit,nbins)):
            nblocks = int(bscores[b] * wlimit/max_bscore +0.5)
            out.write(legends[b] + "|" + "#"*nblocks + " (%d)\n"%bscores[b])

    def score_stats(s):
        distr = {} # distribution of scores
        ldistr = {} # distribution of feat lengths
        score_idx = s.fields.index('score')
        st_idx= s.fields.index('start')
        en_idx = s.fields.index('end')
        nfeat = 0
        for x in s:
            nfeat += 1
            v = x[score_idx]
            distr[v] = distr.get(v,0.0) + 1
            w = x[en_idx]-x[st_idx]
            ldistr[w] = ldistr.get(w,0.0) + 1
        return (nfeat,distr,ldistr,
                stats_from_distr(distr,nfeat), stats_from_distr(ldistr,nfeat), )

    def feat_stats(s):
        ldistr = {} # distribution of feat lengths
        st_idx= s.fields.index('start')
        en_idx = s.fields.index('end')
        for nfeat,x in enumerate(s):
            w = x[en_idx]-x[st_idx]
            ldistr[w] = ldistr.get(w,0.0) + 1
        return nfeat,ldistr,stats_from_distr(ldistr,nfeat)

    def total_coverage(t):
        from bbcflib.gfminer.common import fusion
        total_cov = 0
        kw = dict(kwargs)
        s = t.read(**kw)
        fs = fusion(s)
        st_idx = fs.fields.index('start')
        en_idx = fs.fields.index('end')
        total_cov += sum(x[en_idx]-x[st_idx] for x in fs)
        return total_cov

    if isinstance(source, basestring):
        t = track(source, **kwargs)
    else:
        t = source
    total_cov = total_coverage(t)
    s = t.read(**kwargs)
    is_score = 'score' in s.fields
    if is_score:
        nfeat,distr,ldistr,stat,lstat = score_stats(s)
        if isinstance(out,dict):
            out['feat_stats'] = (nfeat,ldistr,lstat,total_cov)
            out['score_stats'] = (distr,stat)
    else:
        nfeat,ldistr,lstat = feat_stats(s)
        if isinstance(out,dict):
            out['feat_stats'] = (nfeat,ldistr,lstat,total_cov)
    if isinstance(out,dict):
        return out
    if nfeat == 0:
        out.write("Empty content\n\n")
        return
    out.write("  Number of features: %d\n" % nfeat)
    out.write("  Sum of feature lengths: %d\n" % lstat[0])
    out.write("  Total coverage (# of non-null positions): %s\n" % total_cov)
    out.write("  Min length: %.2f\n" % lstat[1])
    out.write("  Max length: %.2f\n" % lstat[2])
    out.write("  Mean length: %.2f\n" % lstat[3])
    out.write("  Lengths standard deviation: %.2f\n" % lstat[4])
    out.write("  Median length: %.2f\n" % lstat[5])
    out.write("Distribution of lengths:\n\n")
    if plot: console_distr_plot(ldistr,out,hlimit=20,wlimit=wlimit,binw=200)
    out.write('\n')
    if is_score:
        out.write("-----------\n")
        out.write("  Total score: %.2f\n" % stat[0])
        out.write("  Min score: %.2f\n" % stat[1])
        out.write("  Max score: %.2f\n" % stat[2])
        out.write("  Mean score: %.2f\n" % stat[3])
        out.write("  Scores standard deviation: %.2f\n" % stat[4])
        out.write("  Median score: %.2f\n" % stat[5])
        out.write("Distribution of scores:\n\n")
        if plot: console_distr_plot(distr,out,hlimit=10,wlimit=wlimit,binw=1)
    out.write('\n')
    return out

################################################################################

class Track(object):
    """
    Metaclass regrouping the track properties. Subclasses for each specific format
    are respectively in `track/text.py`, `track/bin.py`, `track/sql.py`,
    and are instanciated when ``track.track()`` is called on a file.

    .. attribute:: path

        Path to the file the Track was generated from.

    .. attribute:: filehandle

        The Python opened file object from the file found in *self.path*. Can read() and write() it.

    .. attribute:: format

        Format of the file the track was generated from.

    .. attribute:: fields

        Fields defining the info contained in the track items.

    .. attribute:: assembly

        GenRep assembly ID.

    .. attribute:: chrmeta

        A dictionary with information about the species' chromosomes, or a genrep assembly name.

    .. attribute:: info

        A dictionary with meta-data about the track, e.g. data type, such as::

            {'datatype': 'signal'}

    """
    def __init__(self, path, **kwargs):
        self.path = path
        self.filehandle = None
        self.format = kwargs.get("format")
        self.fields = kwargs.get("fields",[])
#        self.types = {}
        self.assembly = kwargs.get('assembly')
        self.chrmeta = self._get_chrmeta(kwargs.get('chrmeta'))
        self.info = self._get_info(info=kwargs.get('info'))
        self.index = {}

    def _get_chrmeta(self,chrmeta=None):
        """:param chrmeta: (str or dict) assembly name, or dict of the type {chr: {'length': 1234}}."""
        if isinstance(chrmeta,dict):
            return chrmeta
        if isinstance(chrmeta,basestring) and not(str(chrmeta) == "guess"):
            self.assembly = chrmeta
        if self.assembly is None:
            return {}
        from bbcflib import genrep
        if genrep.GenRep().assemblies_available(self.assembly):
            self.assembly = genrep.Assembly(self.assembly)
            return self.assembly.chrmeta
        else:
            self.assembly = None
            return {}

    def _get_info(self,info=None):
        return info or {}

    def __enter__(self):
        """Called when evaluating the 'with' statement."""
        return self

    def __exit__(self, errtype, value, traceback):
        self.close()

    def __iter__(self):
        """Iterates over the list of features."""
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

    def readline(self, **kw):
        return self.read(**kw).next()

    def write(self, **kw):
        pass

    def column_by_name(self, fields=[], num=True):
        """
        Finds a column with name in `fields`.
        Returns its index (if `num`is True) or its name.
        """
        if isinstance(fields,basestring): fields=[fields]
        _f = [f for f in fields if f in self.fields]
        if len(_f) == 0: return None
        if num: return self.fields.index(_f[0])
        else: return _f[0]

    @property
    def name(self):
        "Returns an appropriate name for the track"
        name = os.path.basename(self.path)
        if name.endswith(".gzip"): name = name[:-5]
        if name.endswith(".gz"): name = name[:-3]
        return self.info.get('name',os.path.splitext(name)[0])


################################################################################
class FeatureStream(object):
    """
    Contains an iterator yielding features, and an extra fields attribute.
    It can be constructed from either an iterator, a cursor, a list or a tuple.

    Example::

        stream = FeatureStream([('chr',1,2),('chr',3,4)])
        stream = FeatureStream((('chr',1,2),('chr',3,4)))
        stream = FeatureStream(iter([('chr',1,2),('chr',3,4)]))

        def gen():
            for k in range(2):
                yield ('chr',2*k+1,2*k+2)

        stream = FeatureStream(gen())

    Example of usage::

        >>> stream = FeatureStream([('chr',1,2),('chr',3,4)], fields=['chromosome','start','end'])
        >>> stream.next()
        ('chr', 1, 2)
        >>> stream.next()
        ('chr', 3, 4)
        >>> stream.data
        <listiterator object at 0x10183b650>

        >>> stream = FeatureStream([('chr',1,2),('chr',3,4)], fields=['chromosome','start','end'])
        >>> for s in stream: print s
        ('chr', 1, 2)
        ('chr', 3, 4)

    .. attribute:: data

        An iterator, cursor, list or tuple. Each item is a tuple with as many members as the number of *fields*.

    .. attribute:: fields

        The list of field names.

    .. method:: __iter__()

        ``iter(self)`` returns self.data, which is an iterator itself.

    .. method:: next()

        Iterating over the stream is iterating over its data.

    """

    def __init__(self, data, fields=None):
        if isinstance(data,(list,tuple)):
            data = iter(data)
        if not fields:
            if hasattr(data, 'description'):
                fields = [x[0] for x in data.description]
            else: raise ValueError("Must specify a 'fields' attribute for %s." % self.__str__())
        self.data = data
        self.fields = fields

    def __iter__(self):
        return self.data

    def next(self):
        return self.data.next()

################################################################################
