from bbcflib.track import *
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

       Dictionary with keys field names and values functions that will be called on each item when reading the file
       (e.g. check the entry type).

    .. attribute:: outtypes

       Dictionary with keys field names and values functions that will be called on each item when writing the file
       (e.g. format numerics to strings).

    .. attribute:: header

       Indicates the presence of a header.
       * `None` to skip all consecutive lines starting with "browser", "track" or "#" (default).
       * `False` if there is no header.
       * `True` if it is made of a standard unique line with the same number of fields as
         the rest of the file (R-like). Then the header will be used to guess the track fields.
       * An int N to indicate that the first N lines of the file should be skipped.
       * A string to indicate that all header lines start with this string.

    .. attribute:: written

       Boolean indicating whether the self.filehandle has already been written.
       If it has, the default writing mode chages from 'write' to 'append'
       (used after writing a header, for instance).

    When reading a file, all lines beginning with "browser", "track" or "#" are skipped.
    The *info* attribute will be filled with "key=value" pairs found on a "track" line at the top of the file.
    The *open* method takes the argument *mode* which can be 'read' (default), 'write', 'append' or 'overwrite'.
    Path can also be a url, or a gzipped file.
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = kwargs.get("format",'txt')
        self.separator = kwargs.get('separator',"\t")
        self.header = kwargs.get('header',None)
        Track.__init__(self,path,**kwargs)
        if os.path.exists(self.path) and os.path.getsize(self.path): 
            self.separator = self._check_sep()
        self.fields = self._get_fields(kwargs.get('fields'))
        self.intypes = dict((k,v) for k,v in _in_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('intypes'),dict): self.intypes.update(kwargs["intypes"])
        self.outtypes = dict((k,v) for k,v in _out_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('outtypes'),dict): self.outtypes.update(kwargs["outtypes"])
        self.written = False

    def _check_sep(self):
        """Checks if default separator works, otherwise tries ' ' and '\t'."""
        if self.separator is None: return None
        if type(self.filehandle) == file and not self.filehandle.closed:
            return None
        self.open('read')
        self._skip_header()
        while 1:
            row = self.filehandle.readline()
            if not row: break
            if row.count(self.separator) > 1: break
            if row.count(" ") > 1:
                self.separator = " "
                break
            if row.count("\t") > 1:
                self.separator = "\t"
                break
        self.close()
        return self.separator

    def _get_fields(self,fields):
        """R-like fields guessing according if `header=True` is specified.
        Else just return keyword arguments or default value."""
        if fields:
            pass
        elif self.header is True:
            if type(self.filehandle) == file and not self.filehandle.closed:
                return ['chr','start','end']
            self.open()
            header = self.filehandle.readline().split(self.separator)
            self._skip_header()
            first = self.filehandle.readline().split(self.separator)
            if len(header) == len(first):
                fields = [header[0].strip('# ')]+header[1:-1]+[header[-1][:-1]]
            self.close()
        else:
            fields = ['chr','start','end']
        return fields

    def _get_chrmeta(self,chrmeta=None):
        """
        :param chrmeta: (str or dict) assembly name, or dict of the type {chr: {'length': 1234}}.
            If set to `"guess"`, the whole file will be read in order to find the chromosome names
            and lengths.
        """
        _chrmeta = Track._get_chrmeta(self,chrmeta)
        if _chrmeta or not(os.path.exists(self.path) and os.path.getsize(self.path)):
            return _chrmeta
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
        if not(os.path.exists(self.path) and os.path.getsize(self.path)):
            return {}
        self.open()
        _info = {}
        if type(self.filehandle) == file and 'r' not in self.filehandle.mode:
            return _info
        while 1:
            row = self.filehandle.readline()
            if not row: break
            if row[:7]=="browser" or row[0] in ['#','@']: continue
            if row[:5]=="track":
                tok = row[5:].strip().split("=")
                key = tok.pop(0).strip("\'\" ")
                val = tok.pop().strip("\'\" ")
                for t in tok:
                    if t[0] in ["'",'"']:
                        p = re.search(t[0]+r'(.*)'+t[0]+r'\s+(\S+)',t)
                    else:
                        p = re.search(r'(\S+)\s+(\S+)',t)
                    if p:
                        _info[key] = p.groups()[0]
                        key = p.groups()[1].strip("\'\" ")
                    else:
                        key = ''
                        break
                if key and val: _info[key] = val
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
        :param selection: dict of the form {field_name: value} or {field_name: (min,max)}.
        :rtype: boolean
        """
        tests = []
        for k,v in selection.iteritems():
            if k == 'length':
                fi1 = self.fields.index('start')
                fi2 = self.fields.index('end')
            else:
                fi = self.fields.index(k)
            if isinstance(v,(list,tuple)):
                if k == 'chr':
                    tests.append(str(row[fi]) in v)
                elif k == 'length':
                    tests.append(int(row[fi2])-int(row[fi1]) >= int(v[0]))
                    tests.append(int(row[fi2])-int(row[fi1]) <= int(v[1]))
                else:
                    tests.append(float(row[fi]) >= float(v[0]))
                    tests.append(float(row[fi]) <= float(v[1]))
            else:
                if k == 'length':
                    tests.append(int(row[fi2])-int(row[fi1]) == int(v))
                else:
                    tests.append(str(row[fi]) == str(v))
        return all(tests)

    def _index_chr(self,start,end,splitrow):
        chr = splitrow[self.fields.index('chr')]
        if self.index.get(chr): self.index[chr][1] = end
        else:                   self.index[chr] = [start,end]

    def _init_skip(self,selection):
        chr_selected = [s.get('chr')[0] for s in selection if s.get('chr')]
        chr_toskip = sorted([self.index.get(c) for c in self.index if c not in chr_selected])
        chr_toskip = iter(chr_toskip+[[sys.maxint,sys.maxint]])
        return chr_toskip

    def _skip(self,start,next_toskip,chr_toskip):
        if start >= next_toskip[0] and start <= next_toskip[1]:
            self.filehandle.seek(next_toskip[1])
            start = self.filehandle.tell()
            next_toskip = chr_toskip.next()
        end = self.filehandle.tell()
        return start,end,next_toskip

    def open(self,mode='read'):
        """Unzip (if necessary) and open the file with the given *mode*.

        :param mode: (str) one of 'read', 'write', 'append' or 'overwrite'
        """
        if self.filehandle and not self.filehandle.closed: return
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

    def _skip_header(self):
        """If *header* is an int N, skip the first N lines. If it is a string or a list of strings,
        skip the first consecutive lines starting with strings in *header*. If `True`, skips one.
        Else, does nothing."""
        p = 0
        header = self.header
        if isinstance(header,str):
            header = [header]
        if header is True: # skip 1 line
            self.filehandle.readline()
            p = self.filehandle.tell()
        elif isinstance(header, int): # skip N lines
            for k in range(header):
                row = self.filehandle.readline()
            p = self.filehandle.tell()
        elif isinstance(header,(list,tuple)): # skip all lines starting with *str*
            row = header[0]
            while any(row.startswith(h) for h in header):
                p = self.filehandle.tell()
                row = self.filehandle.readline()
        self.filehandle.seek(p)
        row = "#"
        while row and (row[0] in ['#','@'] or row[:5]=='track' or row[:7]=='browser'):
            p = self.filehandle.tell()
            row = self.filehandle.readline()
        self.filehandle.seek(p)

    def _read(self, fields, index_list, selection, skip):
        self.open('read')
        if skip and selection:
            chr_toskip = self._init_skip(selection)
            next_toskip = chr_toskip.next()
        fstart = fend = 0
        self._skip_header()
        try:
            while 1:
                fstart = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row[0] in ['#','@']: continue
                if row[:5]=='track' or row[:7]=='browser': break
                splitrow = [s.strip() for s in row.split(self.separator)]
                if not any(splitrow): continue
                if selection:
                    if skip:
                        fstart,fend,next_toskip = self._skip(fstart,next_toskip,chr_toskip)
                        self._index_chr(fstart,fend,splitrow)
                    if not any(self._select_values(splitrow,s) for s in selection):
                        continue
                fstart = fend
                yield tuple(self._check_type(splitrow[index_list[n]],f)
                            for n,f in enumerate(fields))
        except (ValueError,IndexError) as ve:
            raise ValueError("Bad line in file %s:\n %s%s\n" % (self.path,row,ve))

    def read(self, selection=None, fields=None, skip=False, **kw):
        """
        :param selection: list of dict of the type
            `[{'chr':'chr1','start':(12,24)},{'chr':'chr3','end':(25,45)},...]`,
            where tuples represent ranges, or a FeatureStream.
        :param fields: (list of str) list of field names (columns) to read.
        :param skip: (bool) assuming that lines are grouped by chromosome name,
            increases reading speed when looping over selections of several/all chromosomes.
            The first time lines corresponding to a chromosome are read, their position in the
            file is recorded (self.index). In the next iterations, only lines corresponding to
            chromosomes either yet unread or present in *selection* will be read. [False]
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
        Prepares for writing:
        formats a row (tuple) into a well-formatted string according to the output track format.

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
        initial_separator = self.separator
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
        self.separator = initial_separator
        self.close()

    def make_header(self, *args, **kw):
        """
        If *self* is an empty track, this function can be used to write a header in place
        of the first line of its related file.
        Info can be given as a dictionary *info*, as keyword arguments to the function,
        or as a string. If *info* is a string, it will be written as it is.
        In other cases, the header line will start with 'track' and
        each pair of key/value in *info* is added as 'key=value'. Example::

            make_header(type='bedGraph',name='aaa')
            # or
            make_header(info={'type':'bedGraph','name':'aaa'})

        writes on top of the file::

            track type=bedgraph name=aaa

        :param info: (dict) information to be written. Keys can be:
            'name','description','visibility','color','itemRgb'.
        :param mode: (str) writing mode - one of 'write','overwrite','append'.
        """
        self.open(kw.get('mode','write'))
        if len(args)>0:
            info = args[0]
        else:
            info = kw.get('info')
        if isinstance(info,basestring):
            self.filehandle.write(info+'\n')
            self.header = len(info.split('\r\n'))
        else:
            if isinstance(info,dict): self.info.update(info)
            else: self.info.update(kw)
            header = "track "
            _keys = ["name","type","description","visibility","color","itemRgb"]
            header += " ".join(["%s='%s'"%(k,self.info[k]) for k in _keys if k in self.info])
            self.filehandle.write(header+"\n")
            self.header = 1
        self.close()
        self.written = True

    def get_range(self, selection=None, fields=None):
        """
        Returns the range of values for the given selection.
        If `fields` is None, returns min and max positions, otherwise min and max
        field values.
        """
        if fields is None:
            _f = ["start","end"]
            rback = [None,None]
            for row in self.read(selection=selection,fields=_f):
                if rback[0] is None: rback = list(row)
                if rback[0] > row[0]: rback[0] = row[0]
                if rback[1] < row[1]: rback[1] = row[1]
        else:
            if isinstance(fields,basestring): fields = [fields]
            _f = [f for f in fields if f in self.fields]
            if len(_f) == 0:
                raise ValueError("Fields %s not in track: %s" % (fields,self.fields))
            rback = [None,None]*len(_f)
            for row in self.read(selection=selection,fields=_f):
                if rback[0] is None:
                    for n,x in enumerate(row):
                        rback[2*n] = x
                        rback[2*n+1] = x
                for n,x in enumerate(row):
                    if rback[2*n] > x: rback[2*n] = x
                    if rback[2*n+1] < x: rback[2*n+1] = x
        return rback

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
        if not(os.path.exists(self.path) and os.path.getsize(self.path)): return
        if type(self.filehandle) == file and not 'r' in self.filehandle.mode:
            return
        self.open()
        rowlen = None
        while 1:
            row = self.filehandle.readline()
            if not(row.strip()) or row[0]=='#': continue
            if row[:5]=='track' or row[:7]=='browser': continue
            rowlen = len(row.split(self.separator))
            break
        self.filehandle.seek(0)
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

        ['chr','start','end','name','strand','score'] (when read)

        ['chr','name','end','strand','score']         (when written)

    Scores are rounded to the upper integer when written (but are supposed to be integer originally).
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'sga'
        kwargs['fields'] = ['chr','start','end','name','strand','score']
        def _sga_strand(num=0):
            if num in ['+','-']: return num
            num = int(num)
            if num > 0: return '+'
            if num < 0: return '-'
            return '0'
        def _format_score(x=0.0):
            return "%i" % (x+0.5,) # round to upper int
        kwargs['outtypes'] = {'strand': _sga_strand, 'score': _format_score}
        TextTrack.__init__(self,path,**kwargs)

    def _read(self, fields, index_list, selection, skip):
        self.open('read')
        if skip and selection:
            chr_toskip = self._init_skip(selection)
            next_toskip = chr_toskip.next()
        fstart = fend = 0
        translate_strand = {'+':1, '-':-1, '0':0, 1:1, -1:-1, 0:0}
        sga_fields = ['chr','name','end','strand','score']
        while True:
            fstart = self.filehandle.tell()
            row = self.filehandle.readline()
            if not row: break
            if row[0]=="#": continue
            splitrow = [self._check_type(s.strip(),sga_fields[n])
                        for n,s in enumerate(row.split(self.separator))]
            if not any(splitrow): continue
            chrom,name,pos,strand,score = splitrow
            strand = translate_strand[strand]
            rowdata = (chrom,pos-1,pos,name,strand,score)
            if selection:
                if skip:
                    fstart,fend,next_toskip = self._skip(fstart,next_toskip,chr_toskip)
                    self._index_chr(fstart,fend,rowdata)
                if not any(self._select_values(rowdata,s) for s in selection):
                    continue
            fstart = fend
            yield tuple(rowdata[ind] for ind in index_list)

    def _format_fields(self,vec,row,source_list,target_list):
        """'Bucher' conversion expecting the source to be a result of `bam2wig -q 1`.
        Each entry represents a read start."""
        rowres = ['',0,0,'',0,0]
        for k,n in enumerate(source_list):
            rowres[target_list[k]] = row[n]
        feat = []
        for pos in range(int(rowres[1]),int(rowres[2])):
            x = [rowres[0],rowres[3] or '--',
                 self.outtypes.get("start",str)(pos+1),
                 self.outtypes.get("strand",str)(rowres[4]),
                 self.outtypes.get("score",str)(rowres[5])]
            if x[4] == '0': continue
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
            chr_toskip = self._init_skip(selection)
            next_toskip = chr_toskip.next()
        fstart = fend = 0
        fixedStep = None
        chrom = start = end = step = score = None
        span = 1
        row = ""
        try:
            rowdata = ['',-1,-1,'']
            while 1:
                fstart = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row[0]=="#": continue
                if row[:7]=="browser" or row[:5]=="track":
                    if rowdata[1] >= 0:
                        if (not selection) or any(self._select_values(rowdata,s) for s in selection):
                            yield tuple(self._check_type(rowdata[index_list[n]],f)
                                        for n,f in enumerate(fields))
                    fixedStep = None
                    chrom = start = end = score = step = None
                    span = 1
                    continue
                if row[:9]=="fixedStep":
                    if rowdata[1] >= 0:
                        if (not selection) or any(self._select_values(rowdata,s) for s in selection):
                            yield tuple(self._check_type(rowdata[index_list[n]],f)
                                        for n,f in enumerate(fields))
                    fixedStep = True
                    rowdata = ['',-1,-1,'']
                    chrom,start = re.search(r'chrom=(\S+)\s+start=(\d+)',row).groups()
                    start = int(start)-1
                    step = 1
                    s_patt = re.search(r'step=(\d+)',row)
                    if s_patt: step = max(1,int(s_patt.groups()[0]))
                    rowdata[0] = chrom
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt: span = max(1,int(s_patt.groups()[0]))
                    end = -1
                    start -= step
                    continue
                if row[:12]=="variableStep":
                    if rowdata[1] >= 0:
                        if (not selection) or any(self._select_values(rowdata,s) for s in selection):
                            yield tuple(self._check_type(rowdata[index_list[n]],f)
                                        for n,f in enumerate(fields))
                    fixedStep = False
                    rowdata = ['',-1,-1,'']
                    chrom = re.search(r'chrom=(\S+)',row).groups()[0]
                    rowdata[0] = chrom
                    start = -1
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt: span = int(s_patt.groups()[0])
                    end = start+span
                    continue
                if fixedStep is None: continue
                splitrow = [s.strip() for s in row.split()]
                if not any(splitrow): continue
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
                    start = int(splitrow[0])-1
                    end = start+span
                    if start == rowdata[2] and score == rowdata[3]:
                        yieldit = False
                        rowdata[2] = end
                if not(yieldit): continue
                fstart = fend
                if selection:
                    if skip:
                        fstart,fend,next_toskip = self._skip(fstart,next_toskip,chr_toskip)
                        self._index_chr(fstart,fend,rowdata)
                    if not any(self._select_values(rowdata,s) for s in selection):
                        rowdata[1] = start
                        rowdata[2] = end
                        rowdata[3] = score
                        continue
                if rowdata[1] >= 0:
                    yield tuple(self._check_type(rowdata[index_list[n]],f)
                                for n,f in enumerate(fields))
                rowdata[1] = start
                rowdata[2] = end
                rowdata[3] = score
            if rowdata[1] >= 0:
                if (not selection) or any(self._select_values(rowdata,s) for s in selection):
                    yield tuple(self._check_type(rowdata[index_list[n]],f)
                                for n,f in enumerate(fields))
        except ValueError as ve:
            raise ValueError("Bad line in file %s:\n %s%s\n" % (self.path,row,ve))
        self.close()
        if fixedStep is None:
            raise IOError("Please specify 'fixedStep' or 'variableStep'.")

    def _format_fields(self,vec,row,source_list,target_list):
        chrom = row[source_list[0]]
        start = self.outtypes.get('start',str)(row[source_list[1]]+1)
        span = int(row[source_list[2]])-int(row[source_list[1]])
        score = self.outtypes.get('score',str)(row[source_list[3]])
        head = ''
        if span != vec[1] or chrom != vec[0]:
            head = " ".join(["variableStep","chrom=%s"%chrom,"span=%i"%span])+"\n"
            vec[1] = span
            vec[0] = chrom
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
        def _gff_score(x):
            if x == '.': return 0.0
            else: return float(x)
        def _gff_frame(x):
            if x == '.': return '.'
            else: return int(x)
        self.intypes.update({'score': _gff_score, 'frame': _gff_frame})
        self.outtypes.pop('score')
        if not(os.path.exists(self.path) and os.path.getsize(self.path)): return
        rowlen = 9
        self.open()
        tostrip = ' \r\n'+(self.separator or '\t')
        while 1:
            row = self.filehandle.readline().strip(tostrip)
            if not row: continue
            if row[0] in ['#','@'] or row[:5]=='track' or row[:7]=='browser': continue
            rowlen = len(row.split(self.separator))
            break
        self.filehandle.seek(0)
        if rowlen > 9 or rowlen < 8:
            raise ValueError("Gff should have 8 or 9 fields.")
        if rowlen == 8:
            #self.intypes.pop('attributes')
            #self.outtypes.pop('attributes')
            self.fields = self.fields[:8]

################################ SAM ##########################################

class SamTrack(TextTrack):
    """
    TextTrack class for SAM files (extension ".sam").

    Fields are::

        ['name','flag','chr','start','end','mapq','cigar','rnext','pnext','tlen','seq','qual']

    according to the SAM `specification <http://samtools.sourceforge.net/SAM1.pdf>_`.
    Here, 'name' is the read name (QNAME); 'chr' holds for the reference sequence name (RNAME);
    'start' is the leftmost mapping position (POS); 'end' is 'start' plus the length of the read.
    These can be followed by a few optional tags that can be specified as follows::

        track("myfile.sam", tags=['XA','MD','NM'])
    """
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'sam'
        kwargs['fields'] = ['name','flag','chr','start','end','mapq','cigar','rnext','pnext','tlen','seq','qual'] \
                           + kwargs.get('tags',[])
        TextTrack.__init__(self,path,**kwargs)
        if not(os.path.exists(self.path) and os.path.getsize(self.path)): return
        self.intypes.update({'flag':int, 'mapq':int, 'pnext':int, 'tlen':int})

    def _read(self, fields, index_list, selection, skip):
        self.open('read')
        if skip and selection:
            chr_toskip = self._init_skip(selection)
            next_toskip = chr_toskip.next()
        fstart = fend = 0
        try:
            while 1:
                fstart = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row[0]=="@": continue
                splitrow = [s.strip() for s in row.split(self.separator)]
                if not any(splitrow): continue
                if selection:
                    if skip:
                        fstart,fend,next_toskip = self._skip(fstart,next_toskip,chr_toskip)
                        self._index_chr(fstart,fend,splitrow)
                    if not any(self._select_values(splitrow,s) for s in selection):
                        continue
                splitrow = splitrow[:4]+[int(splitrow[3])+len(splitrow[9])]+splitrow[4:] # end = start + read length
                fstart = fend
                yield tuple(self._check_type(splitrow[index_list[n]],f) for n,f in enumerate(fields))
        except (ValueError,IndexError) as ve:
            raise ValueError("Bad line in file %s:\n %s%s\n" % (self.path,row,ve))
        self.close()

    def _format_fields(self,vec,row,source_list,target_list):
        for i,j in enumerate(target_list):
            vec[j] = self.outtypes.get(self.fields[j],str)(row[source_list[i]])
        vec.pop(4) # remove 'end'
        return self.separator.join(vec)

################################ Fps ##########################################

class FpsTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'fps'
        kwargs['fields'] = ['chr','start','end','name','score','strand']
        TextTrack.__init__(self,path,**kwargs)

    def _read(self, fields, index_list, selection, skip):
        self.open('read')
        if skip and selection:
            chr_toskip = self._init_skip(selection)
            next_toskip = chr_toskip.next()
        fstart = fend = 0
        try:
            while 1:
                fstart = self.filehandle.tell()
                row = self.filehandle.readline()
                if not row: break
                if row[:2] != "FP": continue
                splitrow = [s.strip() for s in row.split()]
                splitrow = [splitrow[1],int(splitrow[5]),int(splitrow[5])+1,splitrow[0],splitrow[7],splitrow[4][1]]
                if selection:
                    if skip:
                        fstart,fend,next_toskip = self._skip(fstart,next_toskip,chr_toskip)
                        self._index_chr(fstart,fend,splitrow)
                    if not any(self._select_values(splitrow,s) for s in selection):
                        continue
                fstart = fend
                yield tuple(self._check_type(splitrow[index_list[n]],f) for n,f in enumerate(fields))
        except (ValueError,IndexError) as ve:
            raise ValueError("Bad line in file %s:\n %s%s\n" % (self.path,row,ve))
        self.close()

