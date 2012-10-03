from bbcflib.btrack import *
import re, urllib2, gzip, os, sys

_in_types = {'start':        int,
             'end':          int,
             'score':        float,
             'strand':       strand_to_int,
             'thick_start':  int,
             'thick_end':    int,
             'block_count':  int,
             'frame':        int}

_out_types = {'start':  format_int,
              'end':    format_int,
              'score':  format_float,
              'strand': int_to_strand}

################################ GENERIC TEXT ####################################

class TextTrack(Track):
    """
    Generic Track class for text files (extension ".txt" or ".text").
    Additional attributes:

    .. attribute:: separator

       Character separating fields in the file (default "\t").

    .. attribute:: intypes

       Dictionary with keys field names and values functions that will be called on each item when reading the file.

    .. attribute:: outtypes

       Dictionary with keys field names and values functions that will be called on each item when writing the file.

    When reading a file, lines beginning with "browser", "track" or "#" are skipped.
    The *info* attribute will be filled with "key=value" pairs found on a "track" line at the top of the file.
    The *open" method takes the argument *mode" which can be *read* (default), *write*, *append* or *overwrite*.
    Path can be a url, or a gzipped file.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = kwargs.get("format",'txt')
        kwargs['fields'] = kwargs.get("fields",['chr','start','end'])
        self.separator = kwargs.get('separator',"\t")
        Track.__init__(self,path,**kwargs) # super(TextTrack, self).__init__(self,path,**kwargs)
        if kwargs.get('ucsc_to_ensembl',False):
            kwargs.setdefault('outtypes', {})
            kwargs['outtypes']['start'] = ucsc_to_ensembl
        if kwargs.get('ensembl_to_ucsc',False):
            kwargs.setdefault('outtypes', {})
            kwargs['outtypes']['start'] = ensembl_to_ucsc
        self.intypes = dict((k,v) for k,v in _in_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('intypes'),dict): self.intypes.update(kwargs["intypes"])
        self.outtypes = dict((k,v) for k,v in _out_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('outtypes'),dict): self.outtypes.update(kwargs["outtypes"])

    def _get_chrmeta(self,chrmeta=None):
        """
        :param chrmeta: (str or dict) assembly name, or dict of the type {chr: {'length': 1234}}.
            If set to `"guess"`, the whole file will be read in order to find the chromosome names
            and lengths.
        """
        _chrmeta = Track._get_chrmeta(self,chrmeta)
        if _chrmeta or not(os.path.exists(self.path)): return _chrmeta
        elif chrmeta == "guess" and 'chr' in self.fields and 'end' in self.fields:
            self.intypes = {'end': int}
            for row in self.read(fields=['chr','end']):
                if not(row[0] in _chrmeta):
                    _chrmeta[row[0]] = {'length': row[1]}
                elif row[1] > _chrmeta[row[0]]['length']:
                    _chrmeta[row[0]]['length'] = row[1]
        return _chrmeta

    def _get_info(self,info=None):
        """
        Read the header of *self.filehandle* to gather some information about the track.
        Update an existing *info* dict, if provided.
        """
        if info: return info
        if not(os.path.exists(self.path)): return {}
        self.open()
        _info = {}
        for row in self.filehandle:
            if row.startswith("browser") or \
                    row.startswith("#"): continue
            if row.startswith("track"):
                for x in row.strip().split():
                    key_val = re.search(r'(\S+)=(\S+)',x)
                    if key_val: _info[key_val.groups()[0]] = key_val.groups()[1]
            break
        self.close()
        return _info

    def _check_type(self,val,field):
        if field in self.intypes:
            return self.intypes[field](val)
        else:
            return val

    def _select_values(self,row,selection):
        """
        Check whether all elements in a *row* pass through the *selection* filter.

        :row: (list) splitted row from file - elements correspond to fields items.
        :param selection: dict of the form {field_name: value} or {field_name: (start,end)}.
        :rtype: boolean
        """
        tests = []
        for k,v in selection.iteritems():
            fi = self.fields.index(k)
            if isinstance(v,(list,tuple)):
                if k == 'chr':
                    tests.append(str(row[fi]) in v)
                else:
                    tests.append(float(row[fi]) >= float(v[0]))
                    tests.append(float(row[fi]) <= float(v[1]))
            else:
                tests.append(str(row[fi]) == str(v))
        return all(tests)

    def open(self,mode='read'):
        """Unzip (if necessary) and open the file with the given *mode*.

        :param mode: (str) one of 'read', 'write', 'append' or 'overwrite'
        """
        isgzip = False
        if self.path.endswith(".gz") or self.path.endswith(".gzip"):
            isgzip = True
        if mode == 'read':
            if os.path.exists(self.path):
                if isgzip:
                    self.filehandle = gzip.open(self.path, 'rb')
                else:
                    self.filehandle = open(self.path, 'r')
            elif self.path.startswith("http://") or \
                    self.path.startswith("https://") or \
                    self.path.startswith("ftp://"):
                if isgzip:
                    raise TypeError("No support for gzipped files through URL yet.")
                else:
                    self.filehandle = urllib2.urlopen(self.path)
            else:
                raise ValueError("Could not find file %s."%self.path)
        elif mode in ['write','overwrite','append']:
            if mode == 'write' and os.path.exists(self.path):
                raise ValueError("File %s already exists, use 'overwrite' or 'append' modes."%self.path)
            else:
                if isgzip:
                    if mode == 'append':
                        self.filehandle = gzip.open(self.path, 'ab')
                    else:
                        self.filehandle = gzip.open(self.path, 'wb')
                elif mode == 'append':
                    self.filehandle = open(self.path,'a')
                else:
                    self.filehandle = open(self.path,'w')
        else:
            raise ValueError("Possible modes are 'read', 'write', 'append' and 'overwrite'.")


    def close(self):
        if self.filehandle: self.filehandle.close()

    def _read(self, fields, index_list, selection, skip):
        self.open('read')

        if skip and selection:
            chr_selected = [s.get('chr')[0] for s in selection if s.get('chr')]
            chr_toskip = [self.index.get(c) for c in self.index if c not in chr_selected]
            skip = iter(chr_toskip+[[sys.maxint,sys.maxint]])
            toskip = skip.next()
        start = end = 0

        try:
            while 1:
                if skip and selection:
                    start = self.filehandle.tell()
                    if start >= toskip[0] and start <= toskip[1]:
                        self.filehandle.seek(toskip[1])
                        start = self.filehandle.tell()
                    elif start > toskip[1]:
                        toskip = skip.next()
                    end = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row.startswith("browser") or \
                   row.startswith("track") or \
                   row.startswith("#"): continue
                splitrow = row.strip(' \r\n').split(self.separator)
                if selection:
                    if not any(self._select_values(splitrow,s) for s in selection):
                        continue
                    else:
                        chr_idx = self.fields.index('chr')
                        chr = splitrow[chr_idx]
                        if self.index.get(chr):
                            self.index[chr][1] = end
                        else:
                            self.index[chr] = [start,end]
                start = end
                yield tuple(self._check_type(splitrow[index_list[n]],f)
                            for n,f in enumerate(fields))
        except (ValueError,IndexError) as ve:
            raise ValueError("Bad line in file %s:\n %s%s\n" % (self.path,row,ve))
        self.close()

    def read(self, selection=None, fields=None, skip=False, **kw):
        """
        :param selection: list of dict of the type
            `[{'chr':'chr1','start':(12,24)},{'chr':'chr3','end':(25,45)},...]`,
            where tuples represent ranges, or a FeatureStream.
        :param fields: (list of str) list of field names (columns) to read.
        :param skip: (bool) assuming that lines are grouped by chromosome name,
            chromosomes that have already been read (by selection) will be skipped
            (unless it is present in *selection*) to increase reading speed. [False]
        """
        if fields is None:
            fields = self.fields
            ilist = range(len(self.fields))
        else:
            try:
                ilist = [self.fields.index(f) for f in fields]
            except ValueError:
                raise ValueError("No such field '%s' in track %s."%(f,self.path))
        if isinstance(selection,basestring):
            selection = [selection]
        if isinstance(selection,(list,tuple)) and isinstance(selection[0],basestring):
            selection = {'chr': [str(x) for x in selection]}
        if isinstance(selection,dict):
            selection = [selection]
        if isinstance(selection,FeatureStream):
            chr_idx = selection.fields.index('chr')
            start_idx = selection.fields.index('start')
            end_idx = selection.fields.index('end')
            sel2 = []
            for feat in selection:
                sel2.append({'chr': feat[chr_idx],
                             'start': (-1,feat[end_idx]),
                             'end': (feat[start_idx],sys.maxint)})
            selection = sel2
        return FeatureStream(self._read(fields,ilist,selection,skip),fields)

    def _format_fields(self,vec,row,source_list,target_list):
        """
        :param vec: (list of types) types of the elements according to self.fields.
        :param row: (list of str) splitted row.
        :param source_list: (list of int) list of indices of the given fields in the source stream fields.
        :param target_list: (list of int) list of indices of the given fields in self.fields.
        """
        for i,j in enumerate(target_list):
            vec[j] = self.outtypes.get(self.fields[j],str)(row[source_list[i]])
        return self.separator.join(vec)

    def write(self, source, fields=None, mode='write', chrom=None, **kw):
        """
        Add data to the track. Effectively writes in the related file.

        :param source: (FeatureStream) data to be added to the track.
        :param fields: list of field names.
        :param mode: (str) file opening mode - one of 'write','overwrite','append'. ['write']
        :param chrom: (str) a chromosome name.
        """
        if self.written and mode=='write':
            mode='append'
        if self.separator is None:
            self.separator = "\t"
        if hasattr(source, 'fields'):
            srcfields = source.fields
        elif hasattr(source, 'description'):
            srcfields = [x[0] for x in source.description]
        elif fields is None:
            srcfields = self.fields
        else:
            srcfields = fields
        if fields is None:
            fields = self.fields
        fields = [f for f in fields if f in self.fields and f in srcfields]

        voidvec = [self.outtypes.get(f,str)() for f in self.fields]
        if not(chrom is None):
            if self.chrmeta and not(chrom in self.chrmeta):
                raise ValueError("Chromosome %s not found in %s." %(chrom,self.chrmeta))
            if 'chr' in self.fields:
                voidvec[self.fields.index('chr')] = chrom
        srcl = [srcfields.index(f) for f in fields]
        trgl = [self.fields.index(f) for f in fields]
        self.open(mode)
        if 'chr' in srcfields: chridx = srcfields.index('chr')
        else: chridx = 0
        if kw.get('clip'):
            sidx = srcfields.index('start')
            eidx = srcfields.index('end')
        for row in source:
            if kw.get('clip'):
                chrsize = self.chrmeta.get(chrom, self.chrmeta.get(row[chridx],{})).get('length',sys.maxint)
                start = max(0,row[sidx])
                end = min(row[eidx],chrsize)
                if end <= start: continue
                row = row[:sidx]+(start,)+row[(sidx+1):]
                row = row[:eidx]+(end,)+row[(eidx+1):]
            self.filehandle.write(self._format_fields(voidvec,row,srcl,trgl)+"\n")
        self.written = True
        self.close()

    def make_header(self, info=None, mode='write'):
        """
        If *self* is an empty track, this function can be used to write a header in place
        of the first line of its related file. By default, only `track type=<track type>` is added.

        :param info: (dict) information to be written. Keys can be:
            'name','description','visibility','color','itemRgb'.
        :param mode: (str) writing mode - one of 'write','overwrite','append'.
        """
        if isinstance(info,dict): self.info.update(info)
        header = "track type=%s " % self.format
        _keys = ["name","description","visibility","color","itemRgb"]
        header += " ".join(["%s=%s" %(k,self.info[k])
                            for k in _keys if k in self.info])
        self.open(mode)
        self.filehandle.write(header+"\n")
        self.close()
        self.written = True

################################ Bed ##########################################

class BedTrack(TextTrack):
    """
    TextTrack class for Bed files (extension ".bed").

    Default fields are::

        ['chr','start','end','name','score','strand',
        'thick_start','thick_end','item_rgb',
        'block_count','block_sizes','block_starts']

    This list will be shortened depending on the number of items found in the first line of the file.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bed'
        _allf = ['chr','start','end','name','score','strand',
                 'thick_start','thick_end','item_rgb',
                 'block_count','block_sizes','block_starts']
        _parf = ['chr','start','end']+kwargs.get('fields',_allf)
        _nf = max([n for n,f in enumerate(_allf) if f in _parf])+1
        kwargs['fields'] = _allf[:_nf]
        TextTrack.__init__(self,path,**kwargs)
        if not(os.path.exists(self.path)): return
        self.open()
        rowlen = None
        for row in self.filehandle:
            if row.startswith("browser") or \
                    row.startswith("track") or \
                    row.startswith("#"): continue
            splitrow = row.strip(' \r\n').split(self.separator)
            rowlen = len(splitrow)
            break
        self.close()
        if rowlen is None: return
        [self.intypes.pop(f) for f in self.fields[rowlen:] if f in self.intypes]
        [self.outtypes.pop(f) for f in self.fields[rowlen:] if f in self.outtypes]
        self.fields = self.fields[:rowlen]

################################ BedGraph ##########################################

class BedGraphTrack(TextTrack):
    """
    TextTrack class for BedGraph files (extension ".bedGraph" or ".bedgraph").

    Fields are::

        ['chr','start','end','score']

    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bedGraph'
        kwargs['fields'] = ['chr','start','end','score']
        TextTrack.__init__(self,path,**kwargs)

################################### Sga ############################################

class SgaTrack(TextTrack):
    """
    TextTrack class for `SGA <http://ccg.vital-it.ch/chipseq/sga_specs.html>`_ files (extension ".sga").

    Fields are::

        ['chr','start','end','name','strand','counts']

    Chromosome names are automatically converted to refseq ids if the assembly attribute is specified.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'sga'
        kwargs['fields'] = ['chr','start','end','name','strand','counts']
        def _sga_strand(num=0):
            if num in ['+','-']: return num
            num = int(num)
            if num > 0: return '+'
            if num < 0: return '-'
            return '0'
        def _score_to_counts(x=0.0): return "%i"%float(x)
        kwargs['intypes'] = {'counts': int}
        kwargs['outtypes'] = {'strand': _sga_strand, 'counts': _score_to_counts}
        TextTrack.__init__(self,path,**kwargs)

    def _read(self, fields, index_list, selection, skip):
        self.open('read')

        if skip and selection:
            chr_selected = [s.get('chr')[0] for s in selection if s.get('chr')]
            chr_toskip = [self.index.get(c) for c in self.index if c not in chr_selected]
            skip = iter(chr_toskip+[[sys.maxint,sys.maxint]])
            toskip = skip.next()
        start = end = 0

        rowdata = {'+': ['',-1,-1,'','+',''], # chr,start,end,name,strand,counts
                   '-': ['',-1,-1,'','-',''],
                   '0': ['',-1,-1,'','.','']}
        while 1:
            if skip and selection:
                start = self.filehandle.tell()
                if start >= toskip[0] and start <= toskip[1]:
                    self.filehandle.seek(toskip[1])
                    start = self.filehandle.tell()
                elif start > toskip[1]:
                    toskip = skip.next()
                end = self.filehandle.tell()
            row = self.filehandle.readline()
            if not row: break
            if row.startswith("#"): continue
            splitrow = row.strip(' \r\n').split(self.separator)
            yieldit = True
            chr,name,pos,strand,counts = splitrow
            if float(counts) == 0: continue
            pos = int(pos)
            rowdata[strand][0] = chr
            if pos-1 == rowdata[strand][2] and \
                    counts == rowdata[strand][5] and \
                    name == rowdata[strand][3]:
                rowdata[strand][2] = pos
                yieldit = False
            if selection:
                if not (self._select_values(rowdata[strand],selection)):
                    continue
                else:
                    chr_idx = self.fields.index('chr')
                    chr = splitrow[chr_idx]
                    if self.index.get(chr):
                        self.index[chr][1] = end
                    else:
                        self.index[chr] = [start,end]
            start = end
            if not(yieldit): continue
            if rowdata[strand][1]>=0:
                yield tuple(self._check_type(rowdata[strand][index_list[n]],f)
                            for n,f in enumerate(fields))
            rowdata[strand][1] = pos-1
            rowdata[strand][2] = pos
            rowdata[strand][3] = name
            rowdata[strand][5] = counts
        for rd in rowdata.values():
            if rd[1]>=0:
                yield tuple(self._check_type(rd[index_list[n]],f)
                            for n,f in enumerate(fields))

    def write(self, source, **kw):
        """Rename field 'score' to 'counts' just for writing, then restablish."""
        sidx = -1
        if 'score' in source.fields:
            sidx = source.fields.index('score')
            source.fields[sidx] = 'counts'
        TextTrack.write(self,source,**kw)
        if sidx > -1:
            source.fields[sidx] = 'score'

    def _format_fields(self,vec,row,source_list,target_list):
        """See `TextTrack._format_fields` ."""
        rowres = ['',0,0,'',0,0]
        for k,n in enumerate(source_list):
            rowres[target_list[k]] = row[n]
        rowres[0] = self.chrmeta.get(rowres[0],{}).get('ac',rowres[0])
        feat = []
        for pos in range(int(rowres[1]),int(rowres[2])):
            x = [rowres[0],rowres[3],
                 self.outtypes.get("start",str)(pos+1),
                 self.outtypes.get("strand",str)(rowres[4]),
                 self.outtypes.get("counts",str)(rowres[5])]
            feat.append(self.separator.join(x))
        return "\n".join(feat)

################################### Wig ############################################

class WigTrack(TextTrack):
    """
    TextTrack class for Wig files (extension ".wig").

    Fields are::

        ['chr','start','end','score']

    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'wig'
        kwargs['fields'] = ['chr','start','end','score']
        TextTrack.__init__(self,path,**kwargs)

    def _read(self, fields, index_list, selection, skip):
        self.open('read')

        if skip and selection:
            chr_selected = [s.get('chr')[0] for s in selection if s.get('chr')]
            chr_toskip = [self.index.get(c) for c in self.index if c not in chr_selected]
            skip = iter(chr_toskip+[[sys.maxint,sys.maxint]])
            toskip = skip.next()
        start = end = 0

        fixedStep = None
        chrom = None
        span = 1
        start = None
        end = None
        step = None
        score = None
        try:
            rowdata = ['',-1,-1,'']
            while 1:
                if skip and selection:
                    start = self.filehandle.tell()
                    if start >= toskip[0] and start <= toskip[1]:
                        self.filehandle.seek(toskip[1])
                        start = self.filehandle.tell()
                    elif start > toskip[1]:
                        toskip = skip.next()
                    end = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row.startswith("#"): continue
                row = row.strip(' \r\n')
                if row.startswith("browser") or row.startswith("track"):
                    fixedStep = None
                    chrom = None
                    span = 1
                    start = None
                    end = None
                    score = None
                    step = None
                    continue
                if row.startswith("fixedStep"):
                    fixedStep = True
                    chrom,start,step = re.search(r'chrom=(\S+)\s+start=(\d+)\s+step=(\d+)',row).groups()
                    start = int(start)
                    step = int(step)
                    rowdata[0] = chrom
                    rowdata[1] = -1
                    rowdata[2] = -1
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt:
                        span = int(s_patt.groups()[0])
                    end = -1
                    start -= step
                    continue
                if row.startswith("variableStep"):
                    fixedStep = False
                    chrom = re.search(r'chrom=(\S+)',row).groups()[0]
                    rowdata[0] = chrom
                    rowdata[1] = -1
                    rowdata[2] = -1
                    start = -1
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt:
                        span = int(s_patt.groups()[0])
                    end = start+span
                    continue
                if fixedStep is None: continue
                splitrow = row.split(self.separator)
                yieldit = True
                if fixedStep:
                    score = splitrow[0]
                    start += step
                    end = start+span
                    if start == rowdata[2] and score == rowdata[3]:
                        yieldit = False
                        rowdata[2] = end
                else:
                    score = splitrow[1]
                    start = int(splitrow[0])
                    end = start+span
                    if start == rowdata[2] and score == rowdata[3]:
                        yieldit = False
                        rowdata[2] = end
                if not(yieldit): continue
                if selection:
                    if not (self._select_values(rowdata,selection)):
                        continue
                    else:
                        chr_idx = self.fields.index('chr')
                        chr = splitrow[chr_idx]
                        if self.index.get(chr):
                            self.index[chr][1] = end
                        else:
                            self.index[chr] = [start,end]
                start = end
                if rowdata[1] >= 0:
                    yield tuple(self._check_type(rowdata[index_list[n]],f)
                                for n,f in enumerate(fields))
                rowdata[1] = start
                rowdata[2] = end
                rowdata[3] = score
                if rowdata[1] >= 0:
                    yield tuple(self._check_type(rowdata[index_list[n]],f)
                                for n,f in enumerate(fields))
        except ValueError:
            raise ValueError("Bad line in file %s:\n %s\n" % (self.path,row))
        self.close()
        if fixedStep is None:
            raise IOError("Please specify 'fixedStep' or 'variableStep'.")

    def _format_fields(self,vec,row,source_list,target_list):
        """See `TextTrack._format_fields` ."""
        chrom = row[source_list[0]]
        start = self.outtypes.get('start',str)(row[source_list[1]])
        span = int(row[source_list[2]])-int(row[source_list[1]])
        score = self.outtypes.get('score',str)(row[source_list[3]])
        head = ''
        if span != vec[1]:
            head = self.separator.join(["variableStep","chrom=%s"%chrom,"span=%s"%span])+"\n"
            vec[1] = span
        feat = self.separator.join([start,score])
        return head+feat

################################ GFF ##########################################

class GffTrack(TextTrack):
    """
    TextTrack class for GFF files (extension ".gff" or ".gtf").

    Fields are::

        ['chr','source','name','start','end','score','strand','frame','attributes']

    with 9th field optional.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'gff'
        kwargs['fields'] = ['chr','source','name','start','end','score','strand','frame','attributes']
        TextTrack.__init__(self,path,**kwargs)
        self.intypes.pop('score') # can be '.' (UCSC)
        self.intypes.pop('frame') # can be '.' (UCSC)
        if not(os.path.exists(self.path)): return
        rowlen = 9
        self.open()
        for row in self.filehandle:
            if row.startswith("browser") or \
                    row.startswith("track") or \
                    row.startswith("#"): continue
            splitrow = row.strip(' \r\n').split(self.separator)
            rowlen = len(splitrow)
            try: float(splitrow[5])
            except ValueError: print "Warning: invalid scores in file %s ." % self.path
            try: int(splitrow[7])
            except ValueError: print "Warning: invalid frame in file %s ." % self.path
            break
        self.close()
        if rowlen > 9 or rowlen < 8:
            raise ValueError("Gff should have 8 or 9 fields.")
        if rowlen == 8:
            self.intypes.pop('attributes')
            self.outtypes.pop('attributes')
            self.fields = self.fields[:8]

