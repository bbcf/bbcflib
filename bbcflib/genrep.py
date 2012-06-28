"""
======================
Module: bbcflib.genrep
======================

This module provides an interface to GenRep repositories.
It provides two classes: the ``Assembly`` class provides a representation of a particular entry in GenRep,
the ``GenRep`` class allows to switch to any potential GenRep repository and
handles all queries. To retrieve an ``Assembly`` named ``ce6``, we write::

    from bbcflib import genrep
    a = genrep.Assembly( assembly='ce6' )

To switch to another instance of genrep::

    g = genrep.GenRep( url=my_url, root=my_path )
    if g.assemblies_available( 'ce6' ):
        a = genrep.Assembly( assembly='ce6', genrep=g )

Assemblies in GenRep are also assigned unique integer IDs.  The unique
integer ID for assembly ``ce6`` is 14.  We can use these IDs anywhere
we would use the name, so the third line in the prevous code could
equally well be written::

    a = genrep.Assembly(14)

"""

# Built-in modules #
import urllib2, json, os, re
from datetime import datetime
from operator import itemgetter

# Internal modules #
from bbcflib.common import normalize_url, unique_filename_in
from bbcflib import btrack as track
from bbcflib.bFlatMajor.common import shuffled as track_shuffle

# Other modules #
import sqlite3

# Constants #
default_url = 'http://bbcftools.vital-it.ch/genrep/'
default_root = '/db/genrep'

################################################################################
class Assembly(object):
    def __init__(self, assembly=None, genrep=None, intype=0):
        """A representation of a GenRep assembly.
        To get an assembly from the repository, call the Assembly
        constructor with either the integer assembly ID or the string assembly
        name.  This returns an Assembly object::

            a = Assembly(3)
            b = Assembly('mm9')

        An Assembly has the following fields:

        .. attribute:: id

        An integer giving the assembly ID in GenRep.

        .. attribute:: name

        A string giving the nameassembly of the assembly in GenRep.

        .. attribute:: index_path

        The absolute path to the bowtie index for this assembly.

        .. attribute:: chromosomes

        A dictionary of chromosomes in the assembly.  The dictionary
        values are tuples of the form (chromsome id, RefSeq locus,
        RefSeq version), and the values are dictionaries with the keys
        'name' and 'length'.

        .. attribute:: bbcf_valid

        Boolean.

        .. attribute:: updated_at

        .. attribute:: created_at

        ``datetime`` objects.

        .. attribute:: nr_assembly_id

        .. attribute:: genome_id

        .. attribute:: source_id

        .. attribute:: intype

        All integers. ``intype`` is '0' for genomic data, '1' for exons and '2' for transcripts.

        .. attribute:: source_name

        .. attribute:: md5

        """
        if genrep is None: genrep = GenRep()
        self.genrep = genrep
        self.intype = int(intype)
        if not(assembly is None):
            self.set_assembly(assembly)

    def set_assembly(self, assembly):
        """Reset the Assembly attributes to correspond to *assembly*.
        *assembly* may be an integer giving the assembly ID, or a string giving the assembly name.
        """
        try:
            assembly = int(assembly)
            assembly_info = json.load(urllib2.urlopen(urllib2.Request(
                            """%s/assemblies/%d.json""" % (self.genrep.url, assembly))))
            chromosomes = json.load(urllib2.urlopen(urllib2.Request(
                            """%s/chromosomes.json?assembly_id=%d""" %(self.genrep.url, assembly))))
        except:
            assembly_info = json.load(urllib2.urlopen(urllib2.Request(
                            """%s/assemblies.json?assembly_name=%s""" %(self.genrep.url, assembly))))[0]
            chromosomes = json.load(urllib2.urlopen(urllib2.Request(
                            """%s/chromosomes.json?assembly_name=%s""" %(self.genrep.url, assembly))))
        self._add_info(**dict((str(k),v) for k,v in assembly_info['assembly'].iteritems()))
        root = os.path.join(self.genrep.root,"nr_assemblies/bowtie")
        if self.intype == 1:
            root = os.path.join(self.genrep.root,"nr_assemblies/exons_bowtie")
        elif self.intype == 2:
            root = os.path.join(self.genrep.root,"nr_assemblies/cdna_bowtie")
        self.index_path = os.path.join(root,self.md5)
        for c in chromosomes:
            chrom = dict((str(k),v) for k,v in c['chromosome'].iteritems())
            cnames = chrom.pop('chr_names')
            chrom['name'] = (str(x['chr_name']['value']) for x in cnames \
                             if x['chr_name']['assembly_id'] == self.id).next()
            self._add_chromosome(**chrom)
        return None

    def _add_info(self, **kwargs):
        self.__dict__.update(kwargs)
        self.created_at = datetime.strptime(self.created_at,'%Y-%m-%dT%H:%M:%SZ')
        self.updated_at = datetime.strptime(self.updated_at,'%Y-%m-%dT%H:%M:%SZ')
        self.source_name = str(self.source_name)
        self.md5 = str(self.md5)
        self.chromosomes = {}

    def _add_chromosome(self, **kwargs):
        chr_keys = ['id','refseq_locus','refseq_version']
        key = tuple([kwargs[k] for k in chr_keys])
        self.chromosomes[key] = dict((k,kwargs[k]) for k in kwargs.keys() if not(k in chr_keys))

    def map_chromosome_names(self, names):
        """
        Finds keys in the `chromosomes` dictionary that corresponds to the names or ids given as `names`.
        Returns a dictionary, such as::

            assembly.map_chromosome_names([3,5,6,47])

            {'3': (2701, u'NC_001135', 4),
            '47': None,
            '5': (2508, u'NC_001137', 2),
            '6': (2580, u'NC_001138', 4)}

        """
        if not(isinstance(names,(list,tuple))):
            names = [names]
        url = "%s/chromosomes.json?assembly_name=%s&identifier=" %(self.genrep.url,self.name)
        mapped = {}
        chr_keys = ['id','refseq_locus','refseq_version']
        for n in names:
            chr_info = json.load(urllib2.urlopen(urllib2.Request(url+str(n))))
            if len(chr_info):
                mapped[str(n)] = tuple([chr_info[0]["chromosome"][k] for k in chr_keys])
            else:
                mapped[str(n)] = None
        return mapped

    def get_links(self,params):
        """
        Returns urls to features. Example::

            assembly.get_links({'name':'ENSMUSG00000085692', 'type':'gene'})

        returns the dictionary
        ``{"Ensembl":"http://ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000085692"}`
        If ``params`` is a string, then it is assumed to be the ``name`` parameter with ``type=gene``.
        """
        if isinstance(params,basestring):
            params = {'name':str(params),'type':'gene'}
        request = "&".join(["md5=%s"%self.md5]+[k+"="+v for k,v in params.iteritems()])
        url = urllib2.urlopen(urllib2.Request(
                """%s/nr_assemblies/get_links/%s.json?%s""" %(self.genrep.url,self.nr_assembly_id,request)))
        return url.read()

    def fasta_from_regions(self, regions, out=None, chunk=50000, shuffled=False):
        """
        Get a fasta file with sequences corresponding to the features in the
        bed or sqlite file.

        Returns the name of the file and the total sequence size.

        If 'out' is a (possibly empty) dictionary, will return the filled dictionary,
        if 'regions' is a dictionary {'chr': [[start1,end1],[start2,end2]]}
        or a list [['chr',start1,end1],['chr',start2,end2]],
        will simply iterate through its items instead of loading a track from file.
        """
        if out == None:
            out = unique_filename_in()
        def _push_slices(slices,start,end,name,cur_chunk):
            if end>start:
                slices['coord'].append([start,end])
                slices['names'].append(name)
                cur_chunk += end-start
            return slices,cur_chunk
        def _flush_slices(slices,chrid,chrn,out):
            names = slices['names']
            coord = slices['coord']
            if isinstance(out,str):
                with open(out,"a") as f:
                    for i,s in enumerate(self.genrep.get_sequence(chrid,coord)):
                        f.write(">"+names[i]+"|"+chrn+":"+str(coord[i][0])+"-"+str(coord[i][1])+"\n"+s+"\n")
            else:
                out[chrn].extend([s for s in self.genrep.get_sequence(chrid,coord)])
            return {'coord':[],'names':[]}
        slices = {'coord':[],'names':[]}
        size = 0
        if isinstance(regions,list):
            reg_dict = {}
            for reg in regions:
                chrom = reg[0]
                if not(chrom in reg_dict):
                    reg_dict[chrom] = []
                reg_dict[chrom].append(reg[1:])
            regions = reg_dict
        if isinstance(regions,dict):
            cur_chunk = 0
            for cid,chrom in self.chromosomes.iteritems():
                if not(chrom['name'] in regions): continue
                if isinstance(out,dict): out[chrom['name']] = []
                for row in regions[chrom['name']]:
                    s = max(row[0],0)
                    e = min(row[1],chrom['length'])
                    slices,cur_chunk = _push_slices(slices,s,e,'',cur_chunk)
                    if cur_chunk > chunk:
                        size += cur_chunk
                        slices = _flush_slices(slices,cid[0],chrom['name'],out)
                        cur_chunk = 0
                size += cur_chunk
                slices = _flush_slices(slices,cid[0],chrom['name'],out)
        else:
            with track.track(regions, chrmeta=self.chrmeta) as t:
                cur_chunk = 0
                for cid,chrom in self.chromosomes.iteritems():
                    features = t.read(selection=chrom['name'],
                                      fields=["start","end","name"])
                    if shuffled:
                        features = track_shuffle( features, chrlen=chrom['length'],
                                                  repeat_number=1, sorted=False )
                    for row in features:
                        s = max(row[0],0)
                        e = min(row[1],chrom['length'])
                        name = re.sub('\s+','_',row[2])
                        slices,cur_chunk = _push_slices(slices,s,e,name,cur_chunk)
                        if cur_chunk > chunk:
                            size += cur_chunk
                            slices = _flush_slices(slices,cid[0],chrom['name'],out)
                            cur_chunk = 0
                    size += cur_chunk
                    slices = _flush_slices(slices,cid[0],chrom['name'],out)
        return (out,size)

    def statistics(self, output=None, frequency=False, matrix_format=False):
        """
        Return (di-)nucleotide counts or frequencies for an assembly, writes in file ``output`` if provided.
        Example of result:
        {
            "TT": 13574667
            "GG": 3344762
            "CC": 3365555
            "AA": 13571722
            "A": 32370285
            "TA": 6362526
            "GT": 4841536
            "AC": 4846697
            "N": 0
            "C": 17781115
            "TC": 6228639
            "GA": 6231575
            "CG": 3131283
            "GC: 3340219
            "CT": 5079814
            "AG": 5075950
            "G": 17758095
            "TG": 6206098
            "CA": 6204462
            "AT": 8875914
            "T": 32371931
        }
        Total = A + T + G + C

        if matrix_format is True output is like:
         >Assembly: sacCer2
        1   0.309798640038793   0.308714120881750   0.190593944221299   0.190893294858157
        """
        request = urllib2.Request("%s/nr_assemblies/%d.json?data_type=counts" % (self.genrep.url, self.nr_assembly_id))
        stat    = json.load(urllib2.urlopen(request))
        total   = float(stat["A"]+stat["T"]+stat["G"]+stat["C"])
        if frequency:
            stat = dict((k,x/total) for k,x in stat.iteritems())
        else:
            stat.update
        if output == None:
            return stat
        else:
            with open(output, "w") as f:
                if matrix_format:
                    f.write(">Assembly: %s\n" % self.name)
                    f.write("%s\t%s\t%s\t%s" %(stat["A"],stat["C"],stat["G"],stat["T"]))
                    f.write("\n")
                else:
                    f.write("#Assembly: %s\n" % self.name)
                    [f.write("%s\t%s\n" % (x,stat[x])) for x in ["A","C","G","T"]]
                    f.write("#\n")
                    [[f.write("%s\t%s\n" % (x+y,stat[x+y])) for y in ["A","C","G","T"]] for x in ["A","C","G","T"]]
            return output

    def fasta_path(self, chromosome=None):
        """
        Returns the path to the compressed fasta file, for the whole assembly or for a single chromosome.
        """
        root = os.path.join(self.genrep.root,"nr_assemblies/fasta")
        path = os.path.join(root,self.md5+".tar.gz")
        if chromosome != None:
            chr_id = str(chromosome[0])+"_"+str(chromosome[1])+"."+str(chromosome[2])
            root = os.path.join(self.genrep.root,"chromosomes/fasta")
            path = os.path.join(root,chr_id+".fa.gz")
        elif self.intype == 1:
            root = os.path.join(self.genrep.root,"nr_assemblies/exons_fasta")
            path = os.path.join(root,self.md5+".fa.gz")
        elif self.intype == 2:
            root = os.path.join(self.genrep.root,"nr_assemblies/cdna")
            path = os.path.join(root,self.md5+".fa.gz")
        return path

    def get_sqlite_url(self):
        '''
        Returns the url of the sqlite file containing gene annotations.
        '''
        return '%s/data/nr_assemblies/annot_tracks/%s.sql' %(self.genrep.url, self.md5)

    def sqlite_path(self):
        '''
        Returns the path to the sqlite file containing genes annotations.
        '''
        root = os.path.join(self.genrep.root,"nr_assemblies/annot_tracks")
        return os.path.join(root,self.md5+".sql")

    def get_features_from_gtf(self,h,chr=None,method="dico"):
        '''
        Return a dictionary *data* of the form
        {key:[[values],[values],...]} containing the result of an SQL request which
        parameters are given as a dictionary *h*. All [values] correspond to a line in the SQL.

        :param chr: (str, or list of str) chromosomes on which to perform the request. By default,
        every chromosome is searched.

        Available keys for h, and possible values:
        "keys":       "*,*,..."            (fields to SELECT and pass as a key of *data*)
        "values":     "*,*,..."            (fields to SELECT and pass as respective values of *data*)
        "conditions": "*:#,*:#,..."        (filter (SQL `WHERE`))
        "uniq":       "whatever"           (SQL `DISTINCT` if specified, no matter what the -string- value is)
        "at_pos":     "12,36,45,1124,..."  (to select only features overlapping this list of positions)

        where
        * holds for any column name in the database
        # holds for any value in the database

        Note: giving several field names to "keys" permits to select unique combinations of these fields.
        The corresponding keys of *data* are a concatenation (by ';') of these fields.
        '''
        data = {}
        if isinstance(chr,list): chromosomes = chr
        elif isinstance(chr,str): chromosomes = chr.split(',')
        else: chromosomes = [None]
        if not(method in ["dico","boundaries"]): return data
        for chr_name in chromosomes:
            request = self.genrep.url+"/nr_assemblies/get_%s?md5=%s" %(method,self.md5)
            request += "&".join(['']+["%s=%s" %(k,v) for k,v in h.iteritems()])
            if chr_name: request += "&chr_name="+chr_name
            data.update(json.load(urllib2.urlopen(request)))
        return data

    def get_gene_mapping(self):
        """Return a dictionary {geneID: (geneName, start, end, length, strand, chromosome)}
        Note that the gene's length is not the sum of the lengths of its exons."""
        gene_mapping = {}
        dbpath = self.sqlite_path()
        if os.path.exists(dbpath):
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chr in self.chrnames:
                sql = """SELECT DISTINCT gene_id,gene_name,MIN(start),MAX(end),strand
                       FROM '%s' WHERE (type='exon') GROUP BY gene_id""" %chr
                cursor.execute(sql)
                for g,name,start,end,strand in cursor:
                    gene_mapping[str(g)] = (str(name),start,end,-1,strand,chr)
                # Find gene lengths
                sql = """SELECT DISTINCT gene_id,start,end
                       FROM '%s' WHERE (type='exon') ORDER BY strand,start,end""" %chr
                cursor.execute(sql)
                try: fg,start,fend = cursor.fetchone()  # initialize
                except TypeError: continue
                length = fend-start
                for g,start,end in cursor:
                    if g == fg:                    # if still in the same gene
                        if start >= fend:
                            length += end-start    # new exon
                            fend = end
                        elif end <= fend: pass     # embedded exon
                        else:
                            length += end-fend     # overlapping exon
                            fend = end
                    else:                          # if new gene
                        gmap = list(gene_mapping[str(fg)])
                        gmap[3] = length
                        gene_mapping.update({str(fg):tuple(gmap)})
                        length = end-start
                        fend = end
                        fg = g
                gmap = list(gene_mapping[str(g)])
                gmap[3] = length
                gene_mapping.update({str(g):tuple(gmap)})
        else:
            h = {"keys":"gene_id", "values":"gene_name,start,end,strand", "conditions":"type:exon", "uniq":"1"}
            for chr in self.chrnames:
                resp = self.get_features_from_gtf(h,chr)
                for k,v in resp.iteritems():
                    start = min([x[1] for x in v])
                    end = max([x[2] for x in v])
                    name = str(v[0][0])
                    strand = int(v[0][3])
                    gene_mapping[str(k)] = (name,start,end,-1,strand,chr)
                # Find gene lengths
                resp_iter = resp.iteritems()
                try: fg,init = resp_iter.next() # initialize
                except StopIteration: continue
                name,start,fend,strand,chr = init[0]
                #start+=1
                length = fend-start
                for g,v in resp_iter:
                    v.sort(key=itemgetter(1,2)) # sort w.r.t. start, then end
                    for x in v:
                        start = x[1]; end = x[2]
                        #start+=1
                        if start >= fend:
                            length += end-start # new exon
                            fend = end
                        elif end <= fend: pass  # embedded exon
                        else:
                            length += end-fend  # overlapping exon
                            fend = end
                    gmap = list(gene_mapping[str(g)])
                    gmap[3] = length
                    length = 0
                    gene_mapping.update({str(g):tuple(gmap)})
                    fend = gmap[2]
        return gene_mapping

    def get_transcript_mapping(self):
        """Return a dictionary ``{transcript ID: (gene ID,start,end,length,strand,chromosome)}``"""
        transcript_mapping = {}
        dbpath = self.sqlite_path()
        if os.path.exists(dbpath):
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chr in self.chrnames:
                sql = """SELECT DISTINCT transcript_id,gene_id,MIN(start),MAX(end),SUM(end-start),strand
                       FROM '%s' WHERE (type='exon') GROUP BY transcript_id""" %chr
                cursor.execute(sql)
                for t,g,start,end,length,strand in cursor:
                    transcript_mapping[str(t)] = (str(g),start,end,length,strand,chr)
        else:
            h = {"keys":"transcript_id", "values":"gene_id,start,end,strand", "conditions":"type:exon", "uniq":"1"}
            for chr in self.chrnames:
                resp = self.get_features_from_gtf(h,chr)
                for k,v in resp.iteritems():
                    start = min([x[1] for x in v])
                    end = max([x[2] for x in v])
                    length = sum([x[2]-x[1] for x in v])
                    gid = str(v[0][0])
                    strand = int(v[0][3])
                    transcript_mapping[str(k)] = (gid,start,end,length,strand,chr)
        return transcript_mapping

    def get_exon_mapping(self):
        """Return a dictionary ``{exon ID: ([transcript IDs],gene ID,start,end,strand,chromosome)}``"""
        exon_mapping = {}
        dbpath = self.sqlite_path()
        if os.path.exists(dbpath):
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chr in self.chrnames:
                sql = """SELECT DISTINCT exon_id,transcript_id FROM '%s' WHERE (type='exon') AND exon_id IS NOT NULL""" %chr
                cursor.execute(sql)
                T={}
                for e,t in cursor:
                    T.setdefault(str(e),[]).append(str(t))
                sql = """SELECT DISTINCT exon_id,gene_id,start,end,strand FROM '%s' WHERE (type='exon') AND exon_id IS NOT NULL""" %chr
                cursor.execute(sql)
                for e,g,start,end,strand in cursor:
                    exon_mapping[str(e)] = (T[str(e)],str(g),start,end,strand,chr)
        else:
            h = {"keys":"exon_id", "values":"gene_id,transcript_id,start,end,strand", "conditions":"type:exon", "uniq":"1"}
            for chr in self.chrnames:
                resp = self.get_features_from_gtf(h,chr)
                for k,v in resp.iteritems():
                    start = int(v[0][2])
                    end = int(v[0][3])
                    tid = [str(x[1]) for x in v]
                    gid = str(v[0][0])
                    strand = int(v[0][4])
                    exon_mapping[str(k)] = (tid,gid,start,end,strand,chr)
        return exon_mapping

    def get_exons_in_trans(self):
        """Return a dictionary ``{transcript ID: list of exon IDs it contains}``"""
        exons_in_trans = {}
        dbpath = self.sqlite_path()
        if os.path.exists(dbpath):
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chr in self.chrnames:
                sql = """SELECT DISTINCT transcript_id,exon_id from '%s' WHERE (type='exon') AND exon_id IS NOT NULL""" %chr
                cursor.execute(sql)
                for t,e in cursor:
                    exons_in_trans.setdefault(str(t),[]).append(str(e))
        else:
            h = {"keys":"transcript_id", "values":"exon_id", "conditions":"type:exon", "uniq":"1"}
            data = self.get_features_from_gtf(h)
            for k,v in data.iteritems():
                exons_in_trans[str(k)] = [str(x[0]) for x in v]
        return exons_in_trans

    def get_trans_in_gene(self):
        """Return a dictionary ``{gene ID: list of transcript IDs it contains}``"""
        trans_in_gene = {}
        dbpath = self.sqlite_path()
        if os.path.exists(dbpath):
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chr in self.chrnames:
                sql = """SELECT DISTINCT gene_id,transcript_id FROM '%s' WHERE (type='exon')""" %chr
                cursor.execute(sql)
                for g,t in cursor:
                    trans_in_gene.setdefault(str(g),[]).append(str(t))
        else:
            h = {"keys":"gene_id", "values":"transcript_id", "conditions":"type:exon", "uniq":"1"}
            data = self.get_features_from_gtf(h)
            for k,v in data.iteritems():
                trans_in_gene[str(k)] = [str(x[0]) for x in v]

        return trans_in_gene

    def gene_coordinates(self,id_list):
        """
        Creates a BED-style stream from a list of gene ids.
        """
        dbpath = self.sqlite_path()
        chromlist = self.chrnames
        _fields = ['chr','start','end','name','strand']
        def _sql_query():
            _ids = "','".join(id_list)
            sql1 = "SELECT DISTINCT MIN(start) AS gstart,MAX(end) AS gend,gene_id,gene_name,strand FROM '"
            sql2 = "' WHERE gene_id IN ('%s') GROUP BY gene_id ORDER BY gstart,gend,gene_id" %_ids
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chrom in chromlist:
                cursor.execute(sql1+chrom+sql2)
                for x in cursor:
                    name = "%s|%s" %x[2:4]
                    yield (chrom,int(x[0]),int(x[1]),str(name))+tuple(x[4:])

        def _web_query():
            _ids = "|".join(id_list)
            sort_list = []
            webh = { "names": "gene_id,gene_name,strand,chr_name",
                     "conditions": "gene_id:%s" %_ids,
                     "uniq":"1" }
            resp = self.get_features_from_gtf(webh,method='boundaries')
            for k,v in resp.iteritems():
                if not(v): continue
                start,end = v
                gene_id,gene_name,strand,chr_name = [str(y).strip() for y in k.split('; ')]
                name = "%s|%s" %(gene_id,gene_name)
                sort_list.append((chr_name,start,end,name,int(strand)))
            sort_list.sort()
            for k in sort_list: yield k
        if os.path.exists(dbpath):
            _db_call = _sql_query
        else:
            _db_call = _web_query
        return track.FeatureStream(_db_call(),fields=_fields)


    def annot_track(self,annot_type='gene',chromlist=None,biotype=["protein_coding"]):
        if chromlist is None: chromlist = self.chrnames
        elif isinstance(chromlist,basestring): chromlist = [chromlist]
        dbpath = self.sqlite_path()
        _fields = ['chr','start','end','name','strand']
        nmax = 4
        biosel = ''
        if not(biotype is None):
            biosel = "AND biotype IN ('"+"','".join(biotype)+"')"
        if annot_type == 'gene':
            flist = "gene_id,gene_name,strand"
            sql1 = "SELECT DISTINCT MIN(start) AS gstart,MAX(end) AS gend,"+flist+" FROM '"
            sql2 = "' WHERE type='exon' %s GROUP BY gene_id ORDER BY gstart,gend,gene_id" %biosel
            webh = { "keys": "gene_id",
                     "values": "start,end,"+flist,
                     "conditions": "type:exon",
                     "uniq":"1" }
        elif annot_type in ['CDS','exon']:
            flist = "exon_id,gene_id,gene_name,strand,frame"
            sql1 = "SELECT DISTINCT start,end,"+flist+" FROM '"
            sql2 = "' WHERE type='%s' %s ORDER BY start,end,exon_id" %(annot_type,biosel)
            webh = { "keys": "exon_id",
                     "values": "start,end,"+flist,
                     "conditions": "type:"+annot_type,
                     "uniq":"1" }
            nmax = 5
            _fields += ['frame']
        elif annot_type == 'transcript':
            if biosel:
                biosel = "WHERE "+biosel[4:]
            flist = "transcript_id,gene_name,strand"
            sql1 = "SELECT DISTINCT MIN(start) AS tstart,MAX(end) AS tend,"+flist+" FROM '"
            sql2 = "' %s GROUP BY transcript_id ORDER BY tstart,tend,transcript_id" %biosel
            webh = { "keys": "transcript_id",
                     "values": "start,end,"+flist,
                     "conditions": "type:exon",
                     "uniq":"1" }
        else:
            raise TypeError("Annotation track type %s not implemented." %annot_type)

        def _sql_query():
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            for chrom in chromlist:
                cursor.execute(sql1+chrom+sql2)
                for x in cursor:
                    name = "|".join([str(y) for y in x[2:nmax]])
                    yield (chrom,int(x[0]),int(x[1]),name)+tuple(x[nmax:])
        def _web_query():
            for chrom in chromlist:
                sort_list = []
                for bt in biotype:
                    wh = webh.copy()
                    wh["conditions"]+=",biotype:"+bt
                    resp = self.get_features_from_gtf(wh,chrom)
                    for k,v in resp.iteritems():
                        start = min([x[0] for x in v])
                        end = max([x[1] for x in v])
                        name = "|".join([str(y) for y in v[0][2:nmax]])
                        sort_list.append((start,end,name)+tuple(v[0][nmax:]))
                sort_list.sort()
                for k in sort_list: yield (chrom,)+k
        if os.path.exists(dbpath):
            _db_call = _sql_query
        else:
            _db_call = _web_query
        return track.FeatureStream(_db_call(),fields=_fields)

    def gene_track(self,chromlist=None,biotype=["protein_coding"]):
        return self.annot_track(annot_type='gene',chromlist=chromlist,biotype=biotype)

    def exon_track(self,chromlist=None,biotype=["protein_coding"]):
        return self.annot_track(annot_type='exon',chromlist=chromlist,biotype=biotype)

    def transcript_track(self,chromlist=None,biotype=["protein_coding"]):
        return self.annot_track(annot_type='transcript',chromlist=chromlist,biotype=biotype)

    @property
    def chrmeta(self):
        """
        Returns a dictionary of chromosome meta data looking something like:
        {'chr1': {'length': 249250621},'chr2': {'length': 135534747},'chr3': {'length': 135006516}}
        """
        return dict([(v['name'], dict([('length', v['length'])])) for v in self.chromosomes.values()])

    @property
    def chrnames(self):
        """
        Returns a list of chromosome names
        """
        namelist = [(v['num'],v['name']) for v in self.chromosomes.values()]
        return [x[1] for x in sorted(namelist)]

################################################################################
class GenRep(object):
    def __init__(self, url=None, root=None, config=None, section='genrep'):
        """Create an object to query a GenRep repository.

        GenRep is the in-house repository for sequence assemblies for the
        BBCF in Lausanne.  This is an object that wraps its use in Python
        in an idiomatic way.

        Create a GenRep object with the base URL to the GenRep system, and
        the root path of GenRep's files.  For instance::

            g = GenRep('genrep.epfl.ch', '/path/to/genrep/indices')

        To get an assembly from the repository, call the assembly
        method with either the integer assembly ID or the string assembly
        name.  This returns an Assembly object::

            a = g.assembly(3)
            b = g.assembly('mm9')

        You can also pass this to the Assembly call directly::

            a = Assembly(assembly='mm9',genrep=g)

        """
        if not(config is None):
            if url is None:
                url = config.get(section, 'genrep_url')
            if root is None:
                root = config.get(section, 'genrep_root')
        if url is None: url = default_url
        if root is None: root = default_root
        self.url = normalize_url(url)
        self.root = os.path.abspath(root)

    def assembly(self,assembly,intype=0):
        """ Backward compatibility """
        return Assembly( assembly=assembly, genrep=self, intype=intype )

    def is_up(self):
        """ Check if genrep webservice is available """
        try:
            urllib2.urlopen(self.url + "/nr_assemblies.json", timeout=2)
        except urllib2.URLError:
            return False
        return True

    def assemblies_available(self, assembly=None):
        """
        Returns a list of assemblies available on genrep, or tells if an
        assembly with name ``assembly`` is available.
        """
        request = urllib2.Request(self.url + "/assemblies.json")
        assembly_list = []
        for a in json.load(urllib2.urlopen(request)):
            name = a['assembly'].get('name')
            if name == None: continue
            if name == assembly: return True
            assembly_list.append(name)
        if assembly == None: return assembly_list

    def get_sequence(self, chr_id, coord_list):
        """Parses a slice request to the repository."""
        if len(coord_list) == 0:
            return []
        slices  = ",".join([",".join([str(y) for y in x]) for x in coord_list])
        url     = """%s/chromosomes/%i/get_sequence_part?slices=%s""" % (self.url, chr_id, slices)
        request = urllib2.Request(url)
        return urllib2.urlopen(request).read().split(',')

    def get_genrep_objects(self, url_tag, info_tag, filters = None, params = None):
        """
        Get a list of GenRep objets
        ... attribute url_tag: the GenrepObject type (plural)
        ... attribute info_tag: the GenrepObject type (singular)
        Optionals attributes:
        ... attribute filters: a dict that is used to filter the response
        ... attribute param: to add some parameters to the query
        from GenRep.
        Example:
        To get the genomes related to 'Mycobacterium leprae' species.
        First get the species with the right name:
        species = get_genrep_objects('organisms', 'organism', {'species':'Mycobacterium leprae'})[0]
        genomes = get_genrep_objects('genomes', 'genome', {'organism_id':species.id})
        """
        if not self.is_up(): return []
        if filters is None: filters = {}
        url = '%s/%s.json' % (self.url, url_tag)
        if params is not None:
            url += '?'
            url += '&'.join([k + '=' + v for k, v in params.iteritems()])
        infos = json.load(urllib2.urlopen(url))
        result = []
        for info in infos:
            obj = GenrepObject(info,info_tag)
            if not filters:
                result.append(obj)
            else:
                for k,v in filters.items():
                    if hasattr(obj,k) and getattr(obj,k) == v:
                        result.append(obj)
                        break
        return result

################################################################################
class GenrepObject(object):
    """
    Class wich will reference all different objects used by GenRep
    In general, you should never instanciate GenrepObject directly but
    call a method from the GenRep object.
    """
    def __init__(self, info, key):
        self.__dict__.update(info[key])

    def __repr__(self):
        return str(self.__dict__)

################################################################################

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
