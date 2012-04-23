from bbcflib.btrack import *
import re, urllib2, gzip, os, sys

#### remark: to do ucsc<->ensembl conversion, define _in_type/_out_type for 'start'

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
                if selection and not(any(self._select_values(splitrow,s) for s in selection)): 
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

################################### Sga ############################################

class SgaTrack(TextTrack):
    def __init__(self,path,**kwargs):
        kwargs['format'] = 'sga'
        kwargs['fields'] = ['chr','start','end','name','strand','counts']
        kwargs['intypes'] = {'counts': int}
        TextTrack.__init__(self,path,**kwargs)
        self.chromosomes = {}
        if self.assembly:
            chdict = genrep.Assembly(self.assembly).chromosomes
            self.chromosomes = dict((v['name'],str(key[1])+"."+str(key[2])) 
                                    for k,v in chdict.iteritems()) 

    def _read(self, fields, index_list, selection=None):
        self.open('read')
        rowdata = {'+': ['',-1,-1,'','+',''],
                   '-': ['',-1,-1,'','-',''],
                   '0': ['',-1,-1,'','.','']}
        for row in self.filehandle:
            row = row.strip()
            if not(row) or row.startswith("#"): continue
            row = row.split(self.separator)
            yieldit = True
            strand = row[3]
            rowdata[strand][0] = row[0]
            start = int(row[2])
            counts = row[4]
            name = row[1]
            if start-1 == rowdata[strand][2] and \
                    counts == rowdata[strand][5] and \
                    name == rowdata[strand][3]:
                rowdata[strand][2] = start
                yieldit = False
            if selection and not(self._select_values(rowdata[strand],selection)): 
                continue
            if not(yieldit): continue
            if rowdata[strand][1]>=0:
                yield tuple(self._check_type(rowdata[strand][index_list[n]],f) 
                            for n,f in enumerate(fields))
            rowdata[strand][1] = start-1
            rowdata[strand][2] = start
            rowdata[strand][3] = name
            rowdata[strand][5] = counts
        for rd in rowdata.values():
            if rd[1]>=0:
                yield tuple(self._check_type(rd[index_list[n]],f) 
                            for n,f in enumerate(fields))


    def _format_fields(self,vec,row,source_list,target_list):
        chrom = row[source_list[0]]
        chrom = self.chromosomes.get(chrom,chrom)
        start = row[source_list[1]]
        end = row[source_list[2]]
        name = row[source_list[3]]
        strand = row[source_list[4]]
        counts = row[source_list[5]]
        feat = []
        for pos in range(start,end):
            x = [chrom,name,self.outtypes.get("start",str)(pos+1),
                 self.outtypes.get("strand",str)(strand),
                 self.outtypes.get("counts",str)(counts)]
            feat.append(self.separator.join(x))
        return "\n".join(feat)

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
        end = None
        step = None
        score = None
        try:
            rowdata = ['',-1,-1,'']
            for row in self.filehandle:
                row = row.strip()
                if not(row) or row.startswith("#"): continue
                if row.startswith("browser") or \
                        row.startswith("track"): 
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
                if selection and not(self._select_values(rowdata,selection)): 
                    continue
                if rowdata[1]>=0:
                    yield tuple(self._check_type(rowdata[index_list[n]],f) 
                                for n,f in enumerate(fields))
                rowdata[1] = start
                rowdata[2] = end
                rowdata[3] = score
            if rowdata[1]>=0:
                yield tuple(self._check_type(rowdata[index_list[n]],f) 
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

