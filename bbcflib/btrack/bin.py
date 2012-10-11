from bbcflib.btrack import *
from bbcflib.common import program_exists
import subprocess, tempfile, os


class BinTrack(Track):
    """
    Generic Track class for binary files.
    """
    def __init__(self,path,**kwargs):
        Track.__init__(self,path,**kwargs)
        self.format = kwargs.get("format",'bin')
        self.fields = kwargs.get("fields",['chr','start','end'])

    def _run_tool(self, tool_name, args):
        if not program_exists(tool_name):
            raise OSError("Program not found in $PATH: %s" % tool_name)
        proc = subprocess.Popen([tool_name]+args, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr: raise OSError("%s exited with message: %s" % (tool_name,stderr))

    def _make_selection(self, selection):
        reg = [None,None,None]
        if isinstance(selection, dict):
            if 'chr' in selection:
                reg[0] = selection['chr']
            if 'start' in selection:
                if isinstance(selection['start'],tuple):
                    reg[1:] = [str(selection['start'][0]),
                               str(selection['start'][1])]
                elif 'end' in selection and not(isinstance(selection['end'],tuple)):
                    reg[1:] = [str(selection['start']),
                               str(selection['end'])]
            elif 'end' in selection and isinstance(selection['end'],tuple):
                reg[1:] = [str(selection['end'][0]),
                           str(selection['end'][1])]
        elif isinstance(selection, basestring):
            reg[0] = selection
        return reg


############################# BigWig via UCSC tools ##############################

class BigWigTrack(BinTrack):
    """
    BinTrack class for BigWig files (extension ".bigWig", ".bigwig" or ".bw").

    Fields are::

        ['chr','start','end','score']

    will use *bedGraphToBigWig* (write) and *bigWigToBedGraph* (read) and use
    the BedGraphTrack class.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bigWig'
        kwargs['fields'] = ['chr','start','end','score']
        BinTrack.__init__(self,path,**kwargs)
        self.bedgraph = None
        self.chrfile = None

    def open(self):
        if self.bedgraph is None:
            tmp = tempfile.NamedTemporaryFile(dir='./')
            self.bedgraph = tmp.name
            tmp.close()

    def close(self):
        if self.chrfile and self.bedgraph:
            self._run_tool('bedGraphToBigWig', [self.bedgraph, self.chrfile.name, self.path])
        if not(self.chrfile is None):
            os.remove(self.chrfile.name)
            self.chrfile = None
        if not(self.bedgraph is None):
            os.remove(self.bedgraph)
            self.bedgraph = None

    def read(self, selection=None, fields=None, **kw):
        """
        :param selection: list of dict of the type
            `[{'chr':'chr1','start':(12,24)},{'chr':'chr3','end':(25,45)},...]`,
            where tuples represent ranges.
        :param fields: (list of str) list of field names.
        """
        self.open()
        if not(fields): fields = self.fields
        fields = [f for f in self.fields if f in fields]
        reg = self._make_selection(selection)
        options = []
        if reg[0]: options += ["-chrom="+reg[0]]
        if reg[1]: options += ["-start="+reg[1]]
        if reg[2]: options += ["-end="+reg[2]]
        self._run_tool('bigWigToBedGraph', options+[self.path, self.bedgraph])
        t = track(self.bedgraph,format='bedGraph',chrmeta=self.chrmeta,info=self.info)
        return t.read(selection=selection,fields=fields,**kw)

    def write(self, source, **kw):
        if self.chrfile is None:
            self.chrfile = tempfile.NamedTemporaryFile(dir='./',delete=False)
            for c,v in self.chrmeta.iteritems():
                self.chrfile.write("%s %i\n"%(c,v['length']))
            self.chrfile.close()
        self.open()
        kw['mode'] = 'append'
        with track(self.bedgraph,format='bedgraph',chrmeta=self.chrmeta) as f:
            f.write(source,**kw)

################################ Bam via pysam ################################

try:
    import pysam
    class BamTrack(BinTrack):
        """
        BinTrack class for Bam files (extension ".bam").

        Fields are::

            ['chr','start','end','score','name','strand','flag','qual','tags']

        uses *pysam* to read the binary bam file and extract the relevant fields.
        Write is not implemented in this class.
        """
        def __init__(self,path,**kwargs):
            kwargs['format'] = 'bam'
            kwargs['fields'] = ['chr','start','end','score','name','strand',
                                'flag','qual','tags']
            BinTrack.__init__(self,path,**kwargs)
            self.filehandle = None
            self.open()
            for h in self.filehandle.header["SQ"]:
                self.chrmeta[h["SN"]] = {'length':h["LN"]}
            self.close()
            self.open()

        def open(self):
            self.filehandle = pysam.Samfile(self.path, "rb")
            self.references = self.filehandle.references # reference sequences (@SN)
            self.lengths = self.filehandle.lengths # reference sequences' length (@LN)
            self.fetch = self.filehandle.fetch
            self.pileup = self.filehandle.pileup

        def close(self):
            self.filehandle.close()

        def read(self, selection=None, fields=None, **kw):
            """
            :param selection: list of dict of the type
                `[{'chr':'chr1','start':(12,24)},{'chr':'chr3','end':(25,45)},...]`,
                where tuples represent ranges.
            :param fields: (list of str) list of field names.
            """
            self.open()
            if not(os.path.exists(self.path+".bai")):
                self._run_tool('samtools', ["index",self.path])
            if not(isinstance(selection,(list,tuple))):
                selection = [selection]
            if fields is None: fields = self.fields
            else: fields = [f for f in fields if f in self.fields]
            srcl = [self.fields.index(f) for f in fields]

            def _bamrecord(stream, srcl):
                for sel in selection:
                    reg = self._make_selection(sel)
                    for read in stream.fetch(reference=reg[0],start=reg[1],end=reg[2]):
                        row = [self.filehandle.getrname(read.rname),
                               read.pos, read.pos+read.rlen,
                               read.mapq, read.qname, (read.is_reverse and -1 or 1),
                               read.flag, read.qual, read.tags]
                        yield tuple([row[n] for n in srcl])
                self.close()
            return FeatureStream(_bamrecord(self.filehandle,srcl),fields)

        def write(self, source, fields, **kw):
            raise NotImplementedError("Writing to bam is not implemented.")

        def count(self, regions, on_strand = False):
            """
            Count the number of reads falling in a given set of *regions*.
            Return a FeatureStream with one element per region, its score being the number of reads
            overlapping (even partially) this region.

            :param regions: any iterable over of tuples of the type `(chr,start,end)`. `chr` has to be
                present in the BAM file's header (see `self.references`). `start` and `end` are 0-based
                coordinates, counting from the beginning of feature `chr` (see `self.lengths`).
            :param on_strand: (bool) restrict to reads on same strand as region.
            :rtype: FeatureStream with fields ['chr','start','end','score'].
            """
            class Counter(object):
                def __init__(self,on_str,reg_str):
                    self.n = 0
                    self.on_str = on_str
                    self.reg_str = reg_str
                def __call__(self, alignment):
                    if self.on_str:
                        if (alignment.is_reverse and self.reg_str>0) or \
                           (not alignment.is_reverse and self.reg_str<0): return
                    self.n += 1

            if isinstance(regions,FeatureStream):
                _f = regions.fields + ['score']
                _si = regions.fields.index('strand') if 'strand' in regions.fields else -1
            else:
                if on_strand:
                    _f = ['chr','start','end','strand','score']
                    _si = 3
                else:
                    _f = ['chr','start','end','score']
                    _si = -1
            def _count(regions):
                for x in regions:
                    if _si > 0: _st = x[_si]
                    else: _st = None
                    c = Counter(on_strand,_st)
                    self.filehandle.fetch(*x[:3], callback=c)
                    #The callback (c.n += 1) is executed for each alignment in a region
                    yield x + (c.n,)

            return FeatureStream(_count(regions),fields=_f)

        def coverage(self, region, strand=None):
            """
            Calculates the number of reads covering each base position within a given *region*.
            Return a FeatureStream where the score is the number of reads overlapping this position.

            :param region: tuple `(chr,start,end)`. `chr` has to be
                present in the BAM file's header (see `self.references`). `start` and `end` are 0-based
                coordinates, counting from the beginning of feature `chr` (see `self.lengths`).
            :strand: if not None, computes a strand-specific coverage ('+' or 1 for forward strand,
                '-' or -1 for reverse strand).
            :rtype: FeatureStream with fields ['chr','start','end','score'].
            """
            if strand is not None:
                pplus = self.filehandle.pileup(*region,mask=16)
                if str(strand) in ['-','-1']: pplus = iter([(x.pos,x.n) for x in pplus])
            pboth = self.filehandle.pileup(*region)
            chr,start,end = region
            _f = ['chr','start','end','score']

            def _coverage(pboth,pplus=None):
                s = start
                e = start
                score = 0
                score1 = 0
                if pplus is None:
                    p1 = (end,0)
                else:
                    p1 = pplus.next()
                for p0 in pboth:
                    if p0.pos < s: continue
                    if p0.pos >= end: break
                    while p0.pos >= p1[0]:
                        try:
                            score1 = p1[1]
                            p1 = pplus.next()
                        except StopIteration:
                            p1 = (end,0)
                    if p0.pos == e+1 and p0.n-score1 == score:
                        e += 1
                    else:
                        if score > 0: yield (chr,s,e+1,score)
                        s = p0.pos
                        e = s
                        score = p0.n-score1
                if score > 0:
                    yield (chr,s,e+1,score)
            if str(strand) in ['+','1']:
                return FeatureStream(_coverage(pplus), fields=_f)
            elif str(strand) in ['-','-1']:
                return FeatureStream(_coverage(pboth,pplus), fields=_f)
            else:
                return FeatureStream(_coverage(pboth), fields=_f)


except ImportError: print "Warning: 'pysam' not installed, 'bam' format unavailable."

