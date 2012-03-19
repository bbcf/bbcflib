from bbcflib.btrack import *
import re, urllib2, gzip, os

#### remark: to do ucsc<->ensembl conversion, define _in_type/_out_type for 'start'

_in_types = {'start':        int,
             'end':          int,
             'score':        float,
             'strand':       strand_to_int,
             'thick_start':  int,
             'thick_end':    int,
             'block_count':  int,
             'frame':        int}

_out_types = {'start': format_int,
              'end': format_int,
              'score': format_float,
              'strand': int_to_strand}

################################ GENERIC TEXT ####################################

class TextTrack(Track):
    def __init__(self,path,**kwargs):
        Track.__init__(self,path,**kwargs)
        self.format = kwargs.get("format",'txt')
        self.fields = kwargs.get("fields",['chr','start','end'])
        self.intypes = dict((k,v) for k,v in _in_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('intypes'),dict): self.intypes.update(kwargs["intypes"])
        self.outtypes = dict((k,v) for k,v in _out_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('outtypes'),dict): self.outtypes.update(kwargs["outtypes"])
        self.filehandle = None
        self.separator = kwargs.get('separator',"\t")

    def _get_info(self,info=None):
        if info: return info
        if not(os.path.exists(self.path)): return
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

    def order_stream(self,stream):
        order = ['chr','start','end']
        indices = [stream.fields.index(o) for o in order if o in stream.fields]
        sort_list = []
        feature_list = []
        chrnames = self.chrmeta.keys()
        for n,f in enumerate(stream):
            fi1 = chrnames.index(f[indices[0]])
            sort_list.append((fi1,f[indices[1]],f[indices[1]],n))
            feature_list.append(f)
        sort_list.sort()
        def _sorted_stream(l1,l2,n):
            for t in l1:
                yield l2[t[n]]
        return FeatureStream(_sorted_stream(sort_list,feature_list,len(indices)), 
                             stream.fields)

    def open(self,mode='read'):
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
                self.filehandle = urllib2.urlropen(self.path)
            else: 
                raise ValueError("Couldn't find the file %s."%self.path)
        elif mode in ['write','overwrite','append']:
            if mode == 'write' and os.path.exists(self.path):
                raise ValueError("File %s exists, use 'overwrite' or 'append' modes."%self.path)
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

    def _read(self, fields, index_list, selection=None):
        self.open('read')
        try:
            for row in self.filehandle:
                if row.startswith("browser") or \
                        row.startswith("track") or \
                        row.startswith("#"): continue
                splitrow = row.strip().split(self.separator)
                if selection and not(self._select_values(splitrow,selection)): 
                    continue
                yield tuple(self._check_type(splitrow[index_list[n]],f) 
                            for n,f in enumerate(fields))
        except ValueError:
            raise ValueError("Bad line in file %s:\n %s\n"%(self.path,row))
        self.close()
                
    def read(self, selection=None, fields=None):
        if fields is None:
            fields = self.fields
            ilist = range(len(self.fields))
        else:
            try:
                ilist = [self.fields.index(f) for f in fields]
            except ValueError:
                raise ValueError("No such field %s in %s."%(f,self.path))
        if isinstance(selection,basestring):
            selection = [str(selection)]
        if isinstance(selection,(list,tuple)): 
            selection = {'chr':selection}
        return FeatureStream(self._read(fields,ilist,selection),fields)

    def _format_fields(self,vec,row,source_list,target_list):
        for i,j in enumerate(target_list):
            vec[j] = self.outtypes.get(self.fields[j],str)(row[source_list[i]])
        return self.separator.join(vec)

    def write(self, source, fields=None, mode='write', chrom=None):
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
        else:
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
        for row in source:
            self.filehandle.write(self._format_fields(voidvec,row,srcl,trgl)+"\n")
        self.close()
        
################################ Bed ##########################################

class BedTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bed'
        kwargs['fields'] = kwargs.get('fields',
                                      ['chr','start','end','name','score','strand',
                                       'thick_start','thick_end','item_rgb',
                                       'block_count','block_sizes','block_starts'])
        TextTrack.__init__(self,path,**kwargs)
        if not(os.path.exists(self.path)): return
        self.open()
        rowlen = None
        for row in self.filehandle:
            if row.startswith("browser") or \
                    row.startswith("track") or \
                    row.startswith("#"): continue
            splitrow = row.strip().split(self.separator)
            rowlen = len(splitrow)
            break
        self.close()
        if rowlen is None: return
        [self.intypes.pop(f) for f in self.fields[rowlen:] if f in self.intypes]
        [self.outtypes.pop(f) for f in self.fields[rowlen:] if f in self.outtypes]
        self.fields = self.fields[:rowlen]
        

################################ BedGraph ##########################################

class BedGraphTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'bedGraph'
        kwargs['fields'] = ['chr','start','end','score']
        TextTrack.__init__(self,path,**kwargs)

################################### Wig ############################################

class WigTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'wig'
        kwargs['fields'] = ['chr','start','end','score']
        TextTrack.__init__(self,path,**kwargs)

    def _read(self, fields, index_list, selection=None):
        self.open('read')
        fixedStep = None
        chrom = None
        span = 1
        start = None
        step = None
        try:
            for row in self.filehandle:
                if row.startswith("#"): continue
                if row.startswith("browser") or \
                        row.startswith("track"): 
                    fixedStep = None
                    chrom = None
                    span = 1
                    start = None
                    step = None
                    continue
                if row.startswith("variableStep"): 
                    fixedStep = False 
                    chrom = re.search(r'chrom=(\S+)',row).groups()[0]
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt:
                        span = s_patt.groups()[0]
                    continue
                if row.startswith("fixedStep"):
                    fixedStep = True
                    chrom,start,step = re.search(r'chrom=(\S+)\s+start=(\d+)\s+step=(\d+)',row).groups()
                    s_patt = re.search(r'span=(\d+)',row)
                    if s_patt:
                        span = s_patt.groups()[0]
                    continue
                if fixedStep is None: continue
                splitrow = row.strip().split(self.separator)
                if fixedStep:
                    data = (chrom,start,start+span,splitrow[0])
                    start += step
                else:
                    start = splitrow[0]
                    data = (chrom,start,start+span,splitrow[1])
                if selection and not(self._select_values(data,selection)): 
                    continue
                yield tuple(self._check_type(data[index_list[n]],f) 
                            for n,f in enumerate(fields))
        except ValueError:
            raise ValueError("Bad line in file %s:\n %s\n"%(self.path,row))
        self.close()
        if fixedStep is None:
            raise IOError("Please specify 'fixedStep' or 'variableStep'.")

    def _format_fields(self,vec,row,source_list,target_list):
        chrom = row[source_list[0]]
        start = self.outtypes.get('start',str)(row[source_list[1]])
        span = row[source_list[2]]-row[source_list[1]]
        score = self.outtypes.get('score',str)(row[source_list[3]])
        head = ''
        if span != vec[1]:
            head = self.separator.join(["variableStep","chrom=%s"%chrom,"span=%s"%span])+"\n"
            vec[1] = span
        feat = self.separator.join([start,score])
        return head+feat


################################ GFF ##########################################
class GffTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'gff'
        kwargs['fields'] = ['chr','source','name','start','end','score','strand','frame','attributes']
        TextTrack.__init__(self,path,**kwargs)
        if not(os.path.exists(self.path)): return
        rowlen = 9
        self.open()
        for row in self.filehandle:
            if row.startswith("browser") or \
                    row.startswith("track") or \
                    row.startswith("#"): continue
            splitrow = row.strip().split(self.separator)
            rowlen = len(splitrow)
            break
        self.close()
        if rowlen > 9 or rowlen < 8:
            raise Error("Gff should have 8 or 9 fields.")
        if rowlen == 8:
            self.intypes.pop('attributes')
            self.outtypes.pop('attributes')
            self.fields = self.fields[:8]

