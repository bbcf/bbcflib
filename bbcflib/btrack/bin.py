from bbcflib.btrack import *
import subprocess, tempfile, os


class BinTrack(Track):
    def __init__(self,path,**kwargs):
        Track.__init__(self,path,**kwargs)
        self.format = kwargs.get("format",'bin')
        self.fields = kwargs.get("fields",['chr','start','end'])

    def _run_tool(self, tool_name, args):
        proc = subprocess.Popen([tool_name]+args, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr: raise Exception("%s exited with message: %s"%(tool_name,stderr))

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
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bigWig'
        kwargs['fields'] = ['chr','start','end','score']
        BinTrack.__init__(self,path,**kwargs)
        self.bedgraph = None

    def open(self):
        tmp = tempfile.NamedTemporaryFile(dir='./')
        self.bedgraph = tmp.name
        tmp.close()

    def close(self):
        if os.path.exists(self.bedgraph): os.remove(self.bedgraph)

    def read(self, selection=None, fields=None, **kw):
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
                     chrmeta=self.chrmeta,info=self.info).read(fields=fields)

    def write(self, source, fields=None, **kw):
        chrfile = tempfile.NamedTemporaryFile(dir='./',delete=False)
        for c,v in self.chrmeta.iteritems():
            chrfile.write("%s %i\n"%(c,v['length']))
        chrfile.close()
        if fields is None: fields = self.fields
        if hasattr(source, 'fields'):
            srcfields = source.fields
        elif hasattr(source, 'description'):
            srcfields = [x[0] for x in source.description]
        elif not(fields is None):
            srcfields = fields
        else:
            srcfields = self.fields
        srcl = [srcfields.index(f) for f in fields]
        self.open()
        with open(self.bedgraph,'w') as f:
            for x in source:
                f.write("\t".join([str(x[i]) for i in srcl])+"\n")
        self._run_tool('bedGraphToBigWig', [self.bedgraph, chrfile.name, self.path])
        self.close()
        os.remove(chrfile.name)

################################ Bam via pysam ################################

try:
    import pysam
    class BamTrack(BinTrack):
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
