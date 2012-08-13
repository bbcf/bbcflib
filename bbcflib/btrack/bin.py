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
        proc = subprocess.Popen([tool_name]+args, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr: raise Exception("%s exited with message: %s" % (tool_name,stderr))

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
            if not program_exists('bedGraphToBigWig'):
                raise OSError("Program not found in $PATH: %s" % 'bedGraphToBigWig')
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
        if not program_exists('bigWigToBedGraph'):
            raise OSError("Program not found in $PATH: %s" % 'bigWigToBedGraph')
        self.open()
        if not(fields): fields = self.fields
        fields = [f for f in self.fields if f in fields]
        reg = self._make_selection(selection)
        options = []
        if reg[0]: options += ["-chrom="+reg[0]]
        if reg[1]: options += ["-start="+reg[1]]
        if reg[2]: options += ["-end="+reg[2]]
        self._run_tool('bigWigToBedGraph', options+[self.path, self.bedgraph])
        return track(self.bedgraph,format='bedGraph',
                     chrmeta=self.chrmeta,info=self.info).read(fields=fields,**kw)

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

        def count(self, regions):
            """
            Count the number of reads falling in a given set of *regions*.
            Return a dictionary of the type `{name: count}`.

            :param regions: any iterable over of tuples of the type `(name,start,end)`. `name` has to be
                present in the BAM file's header (see `self.references`). `start` and `end` are 0-based
                coordinates, counting from the beginning of feature `name` (see `self.lengths`).
            :rtype: dict
            """
            class Counter(object):
                def __init__(self):
                    self.n = 0
                def __call__(self, alignment):
                    self.n += 1

            counts = {}
            c = Counter()
            for x in regions:
                self.filehandle.fetch(x[0],x[1],x[2], callback=c)
                #The callback (c.n += 1) is executed for each alignment in a region
                counts[x[0]] = c.n
                c.n = 0
            return counts

        def coverage(self, region):
            """
            Calculates the number of reads covering each base position within a given *region*.
            Return a dict of the form {pos: coverage}

            :param region: tuple `(name,start,end)`. `name` has to be
                present in the BAM file's header (see `self.references`). `start` and `end` are 0-based
                coordinates, counting from the beginning of feature `name` (see `self.lengths`).
            """
            coverage = {}
            pile = self.filehandle.pileup(region[0],region[1],region[2])
            for p in pile:
                pos = p.pos
                if pos >= region[1] and pos < region[2]:
                    coverage[pos] = p.n
            for pos in range(region[1],region[2]):
                coverage.get(pos,0)
            return coverage


except ImportError: print "Warning: 'pysam' not installed, 'bam' format unavailable."

