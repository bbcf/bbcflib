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

        def open(self):
            self.filehandle = pysam.Samfile(self.path, "rb")

        def close(self):
            self.filehandle.close()

        def read(self, selection=None, fields=None, **kw):
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

except ImportError: print "Warning: 'pysam' not installed, 'bam' format unavailable."
