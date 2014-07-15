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
import urllib2, json, os, sys, re, gzip, tarfile, bz2, sqlite3, itertools
from datetime import datetime

# Internal modules #
from bbcflib.common import normalize_url, unique_filename_in, fasta_length, fasta_composition, sam_faidx
from bbcflib.track import track, ensembl_to_ucsc, FeatureStream
from bbcflib.gfminer.common import shuffled as track_shuffle, split_field, map_chromosomes

# Constants #
default_url = 'http://bbcftools.epfl.ch/genrep/'
default_root = '/db/genrep'

################################################################################
class Assembly(object):
    """
    A representation of a GenRep assembly.
    To get an assembly from the repository, call the Assembly
    constructor with either the integer assembly ID or the string assembly
    name.  This returns an Assembly object::

        a = Assembly(3)
        b = Assembly('mm9')

    An Assembly has the following fields:

    .. attribute:: id

       An integer giving the assembly ID in GenRep.

    .. attribute:: name

       A string giving the name of the assembly in GenRep.

    .. attribute:: index_path

       The absolute path to the bowtie/bowtie2/SOAPsplice index for this assembly.

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

       All integers. ``intype`` is '0' for genomic data, '1' for exons, '2' for transcripts, '3' for junctions.

    .. attribute:: source_name

    .. attribute:: md5


    """
    def __init__(self, assembly=None, genrep=None, intype=0,
                 fasta=None, annot=None, ex=None, via='local',
                 bowtie2=False):
        if genrep is None: genrep = GenRep()
        self.genrep = genrep
        self.intype = int(intype)
        self.chromosomes = {}
        if fasta is not None:
            self.build_assembly(ex,assembly,fasta,annot,via,bowtie2=bowtie2)
        elif assembly is not None:
            self.set_assembly(assembly)
            if annot is not None and os.path.exists(annot):
                self.annot_path = self.gtf_to_sql(annot)
        else:
            raise ValueError("Need an assembly or a fasta file.")

    ### This is a hack for SNPs, should go in the database at some point
        _p = {"EB1_e_coli_k12":1,"MLeprae_TN":1,"mycoSmeg_MC2_155":1,
              "mycoTube_H37RV":1,"NA1000":1,"vibrChol1":1,"TB40-BAC4":1,
              "tbR25":1,"pombe":1,"sacCer2":1,"sacCer3":1,"ASM1346v1":1}
        self.ploidy = _p.get(self.name,2)

    def set_assembly(self, assembly):
        """
        Reset the Assembly attributes to correspond to *assembly*.

        :param assembly: integer giving the assembly ID, or a string giving the assembly name.
        """
        try: # Local db
            assembly_json = os.path.join(self.genrep.root,"nr_assemblies/info_json/%s.json" % assembly)
            assembly_info = json.load(open(assembly_json))[0]
            chromosomes_json = os.path.join(self.genrep.root,"nr_assemblies/info_json/%s_chromosomes.json" % assembly)
            chromosomes = json.load(open(chromosomes_json))
        except: # Online GenRep
            try:
                assembly = int(assembly)
                url = """%s/assemblies/%d.json""" % (self.genrep.url, assembly)
                assembly_info = json.load(urllib2.urlopen(urllib2.Request(url)))
                url2 = """%s/chromosomes.json?assembly_id=%d""" %(self.genrep.url, assembly)
                chromosomes = json.load(urllib2.urlopen(urllib2.Request(url2)))
            except:
                try:
                    url = """%s/assemblies.json?assembly_name=%s""" %(self.genrep.url, assembly)
                    assembly_info = json.load(urllib2.urlopen(urllib2.Request(url)))[0]
                    url2 = """%s/chromosomes.json?assembly_name=%s""" %(self.genrep.url, assembly)
                    chromosomes = json.load(urllib2.urlopen(urllib2.Request(url2)))
                except (IndexError, urllib2.URLError):
                    raise ValueError("URL not found: %s." % url)
        self._add_info(**dict((str(k),v) for k,v in assembly_info['assembly'].iteritems()))
        self.set_index_path()
        for c in chromosomes:
            chrom = dict((str(k),v) for k,v in c['chromosome'].iteritems())
            cnames = chrom.pop('chr_names')
            chrom['name'] = (str(x['chr_name']['value']) for x in cnames \
                             if x['chr_name']['assembly_id'] == self.id).next()
            self._add_chromosome(**chrom)
        return None

    def set_index_path(self,intype=None):
        if intype is not None: self.intype = int(intype)
        root = os.path.join(self.genrep.root,"nr_assemblies/bowtie")
        if self.intype == 1:
            root = os.path.join(self.genrep.root,"nr_assemblies/exons_bowtie")
        elif self.intype == 2:
            root = os.path.join(self.genrep.root,"nr_assemblies/cdna_bowtie")
        elif self.intype == 3:
            root = os.path.join(self.genrep.root,"nr_assemblies/soapsplice")
        if self.intype == 3:
            self.index_path = os.path.join(root,self.md5,self.md5+'.index')
        else:
            self.index_path = os.path.join(root,self.md5)

    def build_assembly(self, ex, assembly, fasta, annot, via, bowtie2=False):
        """
        Build an Assembly object from files.
        """
        from bbcflib.mapseq import bowtie_build
        if assembly is None:
            if fasta is not None:
                self.name = os.path.splitext( os.path.basename( fasta ) )[0]
            elif annot is not None:
                self.name = os.path.splitext( os.path.basename( fasta ) )[0]
            else:
                self.name = "_custom_"
            fasta_by_chrom = {}
        else:
            self.set_assembly(assembly)
            fasta_by_chrom = self.untar_genome_fasta()
            self.name = assembly
        if fasta is not None:
            self.fasta_origin = fasta
            chromosomes = {}
            add_chrom = self.untar_genome_fasta()
            fasta_files = list(set(add_chrom.values()))
            for f in fasta_files: chromosomes.update(fasta_length.nonblocking(ex,f,via=via).wait())
            fasta_by_chrom.update( add_chrom )
            fasta_files = list(set(fasta_by_chrom.values()))
            [g.wait() for g in [sam_faidx.nonblocking(ex,f,via=via) for f in fasta_files]]
            if len(fasta_files) > 0:
                self.index_path = bowtie_build.nonblocking(ex,fasta_files,bowtie2=bowtie2,via=via,memory=8).wait()
            self.stats_dict = {}
            for f in fasta_files:
                stats = fasta_composition(ex,f)
                for k,v in stats.iteritems():
                    self.stats_dict[k] = self.stats_dict.get(k,0)+v
            for k,chrom in self.chromosomes.iteritems():
                if chrom['name'] in chromosomes:
                    chrom_ac = str(k[0])+"_"+str(k[1])+"."+str(k[2]) if k[1] else str(k[0])
                    if chrom.get('synonyms'):
                        chromosomes[chrom['name']]['synonyms'] = chrom['synonyms']+","+chrom_ac
                    else:
                        chromosomes[chrom['name']]['synonyms'] = chrom_ac
            self.chromosomes.update(dict(((k,None,None),v)
                                         for k,v in chromosomes.iteritems()))
            self._fasta_by_chrom = fasta_by_chrom
            if hasattr(self,"annot_origin"):
                archive = tarfile.open(fasta)
                input_file = archive.extractfile(self.annot_origin)
                annot = unique_filename_in()
                with open(annot,"w") as outf:
                    [outf.write(line) for line in input_file]
        if annot is not None and os.path.exists(annot):
            self.annot_path = self.gtf_to_sql(annot)

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
        ``{"Ensembl":"http://ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000085692"}``.
        If *params* is a string, then it is assumed to be the *name* parameter with ``type=gene``.
        """
        if isinstance(params,basestring):
            params = {'name':str(params),'type':'gene'}
        request = "&".join(["md5=%s"%self.md5]+[k+"="+v for k,v in params.iteritems()])
        url = urllib2.urlopen(urllib2.Request(
                """%s/nr_assemblies/get_links/%s.json?%s""" %(self.genrep.url,self.nr_assembly_id,request)))
        return url.read()

    def fasta_from_regions(self, regions, out=None, path_to_ref=None, chunk=50000, shuffled=False, ex=None, intype=0):
        """
        Get a fasta file with sequences corresponding to the features in the
        bed or sqlite file.

        Returns a tuple `(out,size)` where `out` is the name of the output file (or a dict)
        and `size` is the total size of the extracted sequence.

        :param regions: (str or dict or list) bed or sqlite file name, or sequence of features.
            If *regions* is a dictionary {'chr': [[start1,end1],[start2,end2]]}
            or a list [['chr',start1,end1],['chr',start2,end2]],
            will simply iterate through its items instead of loading a track from file.
        :param out: (str, filehandle or dict) output file name or filehandle.
            If *out* is a (possibly empty) dictionary, will return the updated dictionary.
        :param path_to_ref: (str or dict) path to a fasta file containing the whole reference sequence,
            or a dictionary {chr_name: path} as returned by Assembly.untar_genome_fasta.
        :param chunk: (int) buffer size (length of the sequence kept in memory before writing).
        :param intype: (int) if 2, only transcribed sequences are returned (slices of mature RNAs).
            In this case, the fasta headers have the form
            ">assembly_id|transcript_id|genomic_coordinates".
            For each of the given *regions*, one sequence per intersecting cDNA sequence is reported. [0]
        :rtype: (str,int) or (dict,int)
        """
        if out is None: out = unique_filename_in()
        _is_filename = False
        if isinstance(out,basestring):
            _is_filename = True
            out = open(out,"w")
        if path_to_ref is None:
            path_to_ref = self.fasta_by_chrom

        def _push_slices(slices,start,end,name,cur_chunk,chrom,idx):
            """Add a feature to *slices*, and increment the buffer size *cur_chunk* by the feature's size."""
            if end > start:
                slices['coord'].append((start,end))
                slices['names'].append(name)
                cur_chunk += end-start
            return slices,cur_chunk

        def _flush_slices(slices,chrid,chrn,out,path_to_ref):
            """Write the content of *slices* to *out*."""
            names = slices['names']
            coord = slices['coord']
            seqs = self.genrep.get_sequence(chrid, coord, path_to_ref=path_to_ref, chr_name=chrn, ex=ex)
            if isinstance(out,file):
                for i,s in enumerate(seqs):
                    if s:
                        out.write(">%s|%s:%d-%d\n%s\n" % (names[i],chrn,coord[i][0],coord[i][1],s))
            else:
                out[chrn].extend(seqs)
            return {'coord':[],'names':[]}

        def _push_transcript_slices(slices,start,end,name,cur_chunk,chrom,idx):
            """Add a feature to *slices*, and increment the buffer size *cur_chunk* by the feature's size."""
            dbpath = self.sqlite_path
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            request = """SELECT DISTINCT transcript_id,exon_id,start,end from '%s'
                         WHERE end>=%d AND start<%d AND type='exon' ORDER BY transcript_id,start,end;""" \
                      % (chrom,start,end)
            cursor.execute(request)
            transcripts = dict((key,tuple((g[2],g[3]) for g in group)) \
                               for key,group in itertools.groupby(cursor, lambda x: x[0]))
            for tid,exon_coords in transcripts.iteritems():
                tstart = exon_coords[0][0] # start of first exon
                tend = exon_coords[-1][1]  # end of last exon
                for exon in exon_coords:
                    st = max(exon[0],start)
                    en = min(exon[1],end)
                    if en > st:
                        slices['names'].append((idx,'%s|%d-%d' % (tid,max(0,start-tstart),min(end-tstart,tend-tstart))))
                        slices['coord'].append((st,en))
                        cur_chunk += en-st
            return slices,cur_chunk

        def _flush_transcript_slices(slices,chrid,chrn,out,path_to_ref):
            """Write the content of *slices* to *out*."""
            nc = itertools.izip(slices['names'],slices['coord'])
            gb = itertools.groupby(nc, lambda x: x[0])
            for name,group in gb:
                coord = tuple(x[1] for x in group)
                s = self.genrep.get_sequence(chrid, coord, path_to_ref=path_to_ref, chr_name=chrn, ex=ex)
                s = ''.join(s)
                if isinstance(out,file) and s:
                    out.write(">%s|%s:%d-%d\n%s\n" % (self.name,chrn,coord[0][0],coord[-1][1],s))
                else:
                    out[chrn].append(s)
            return {'coord':[],'names':[]}

        if intype == 2:
            _push_slices = _push_transcript_slices
            _flush_slices = _flush_transcript_slices
        slices = {'coord':[],'names':[]}
        size = 0
        region_idx = 0
        if isinstance(regions,(list,tuple)): # convert to a dict
            reg_dict = {}
            for reg in regions:
                chrom = reg[0]
                reg_dict.setdefault(chrom, [])
                reg_dict[chrom].append(reg[1:])
            regions = reg_dict
        if isinstance(regions,dict):
            cur_chunk = 0
            for cid,chrom in self.chromosomes.iteritems():
                chrname = chrom['name']
                if not(chrname in regions): continue
                if isinstance(out,dict): out[chrname] = []
                if isinstance(path_to_ref,dict):
                    pref = path_to_ref.get(chrname)
                else:
                    pref = path_to_ref
                for row in regions[chrname]:
                    region_idx += 1
                    s = max(row[0],0)
                    e = min(row[1],chrom['length'])
                    slices,cur_chunk = _push_slices(slices,s,e,self.name,cur_chunk,chrname,region_idx)
                    if cur_chunk > chunk:
                        size += cur_chunk
                        slices = _flush_slices(slices,cid,chrname,out,pref)
                        cur_chunk = 0
                size += cur_chunk
                slices = _flush_slices(slices,cid,chrname,out,pref)
        else:
            with track(regions, chrmeta=self.chrmeta) as t:
                _f = [f for f in ["start","end","name"] if f in t.fields]
                cur_chunk = 0
                for cid,chrom in self.chromosomes.iteritems():
                    chrname = chrom['name']
                    if isinstance(path_to_ref,dict):
                        pref = path_to_ref.get(chrname)
                    else:
                        pref = path_to_ref
                    features = t.read(selection=chrname,fields=_f)
                    if shuffled:
                        features = track_shuffle( features, chrlen=chrom['length'],
                                                  repeat_number=1, sorted=False )
                    for row in features:
                        region_idx += 1
                        s = max(row[0],0)
                        e = min(row[1],chrom['length'])
                        name = re.sub('\s+','_',row[2]) if len(row)>2 else chrname
                        slices,cur_chunk = _push_slices(slices,s,e,name,cur_chunk,chrname,region_idx)
                        if cur_chunk > chunk: # buffer is full, write
                            size += cur_chunk
                            slices = _flush_slices(slices,cid,chrname,out,pref)
                            cur_chunk = 0
                    size += cur_chunk
                    slices = _flush_slices(slices,cid,chrname,out,pref)
        if _is_filename:
            out.close()
            out = out.name
        return (out,size)

    def motifs_available(self):
        return self.genrep.motifs_available(genome_id=self.genome_id)

    def get_motif_PWM(self, motif_name, output=None):
        return self.genrep.get_motif_PWM(self.genome_id, motif_name, output=output)

    def statistics(self, output=None, frequency=False, matrix_format=False, ex=None):
        """
        Return (di-)nucleotide counts or frequencies for an assembly, writes in file *output* if provided.
        Example of result::

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

        If *matrix_format* is True, *output* is like::

             >Assembly: sacCer2
             1   0.309798640038793   0.308714120881750   0.190593944221299   0.190893294858157
        """
        bases = "ACGT"
        if hasattr(self,"stats_dict"):
            stat = self.stats_dict
        else:
            request = urllib2.Request("%s/nr_assemblies/%d.json?data_type=counts" %(self.genrep.url, self.nr_assembly_id))
            stat = json.load(urllib2.urlopen(request))
        total = sum(float(stat[k]) for k in bases)
        if frequency:
            stat = dict((k,x/total) for k,x in stat.iteritems())
        if output == None:
            return stat
        else:
            with open(output, "w") as f:
                if matrix_format:
                    f.write(">Assembly: %s\n" % self.name)
                    f.write("1\t%s\t%s\t%s\t%s" %tuple(stat[x] for x in bases))
                    f.write("\n")
                else:
                    f.write("#Assembly: %s\n" % self.name)
                    [f.write("%s\t%s\n" % (x,stat[x])) for x in bases]
                    f.write("#\n")
                    [[f.write("%s\t%s\n" % (x+y,stat[x+y])) for y in bases] for x in bases]
            return output

    def fasta_path(self, intype=None, chromosome=None):
        """Return the path to the compressed fasta file, for the whole assembly or for a single chromosome."""
        if not intype:
            intype = self.intype
        if hasattr(self,"fasta_origin"):
            if chromosome is not None:
                if chromosome in self._fasta_by_chrom:
                    return self._fasta_by_chrom[chromosome]
                if chromosome[0] in self._fasta_by_chrom:
                    return self._fasta_by_chrom[chromosome[0]]
            return self.fasta_origin
        root = os.path.join(self.genrep.root,"nr_assemblies/fasta")
        path = os.path.join(root,self.md5+".fa.gz")
        if intype == 1:
            root = os.path.join(self.genrep.root,"nr_assemblies/exons_fasta")
            path = os.path.join(root,self.md5+".fa.gz")
        elif intype == 2:
            root = os.path.join(self.genrep.root,"nr_assemblies/cdna")
            path = os.path.join(root,self.md5+".fa.gz")
        elif chromosome is not None:
            chr_id = str(chromosome[0])
            if chromosome[1]: chr_id += "_"+str(chromosome[1])+"."+str(chromosome[2])
            root = os.path.join(self.genrep.root,"chromosomes/fasta")
            path = os.path.join(root,chr_id+".fa")
        return path

    def untar_genome_fasta(self, path_to_ref=None, convert=True):
        """Untar reference sequence fasta files.
        Returns a dictionary {chr_name: file_name}

        :param path_to_ref: (str) path to the fasta file of the reference sequence (possibly .tar).
        :param convert: (bool) True if chromosome names need conversion from id to chromosome name.
        """

        def _rewrite(inf,genomeRef,chrlist):
            newfa = unique_filename_in()
            outf = open(newfa,"w")
            for line in inf:
                headpatt = re.search(r'^>(\S+)\s',line)
                if headpatt:
                    chrom = headpatt.groups()[0]
                    if chrom in chrlist:
                        line = re.sub(chrom,chrlist[chrom],line)
                        chrom = chrlist[chrom]
                    else:
                        chrom = re.sub(r'[;|,/\\]','_',chrom)
                        line = ">"+chrom+"\n"
                    outf.write(line)
                    genomeRef[chrom] = newfa
                else: outf.write(line)
            outf.close()

        if hasattr(self,"_fasta_by_chrom"): return self._fasta_by_chrom
        if path_to_ref is None: path_to_ref = self.fasta_path()
        if not(os.path.exists(path_to_ref)):
            raise ValueError("Reference fasta archive not found: %s."%path_to_ref)
        if convert:
            # Create a dict of the form {'2704_NC_001144.4': chrXII}
            chrlist = dict((v['ac'],k) for k,v in self.chrmeta.iteritems())
        else:
            chrlist = {}
        genomeRef = {}
        is_gz = path_to_ref.endswith((".gz",".gzip"))
        is_bz2 = path_to_ref.endswith((".bz2",".bz"))
        if ".tar" in os.path.basename(path_to_ref) or ".tgz" in os.path.basename(path_to_ref):
            archive = tarfile.open(path_to_ref)
            for f in archive.getmembers():
                if f.isdir(): continue
                if os.path.splitext(f.name)[0][1:4] in ["sql","gtf"]:
                    self.annot_origin = f.name
                    continue
                input_file = archive.extractfile(f)
                _rewrite(input_file,genomeRef,chrlist)
                input_file.close()
            archive.close()
        elif is_gz:
            input_file = gzip.open(path_to_ref,mode='rb')
            _rewrite(input_file,genomeRef,chrlist)
            input_file.close()
        elif is_bz2:
            input_file = bz2.BZ2File(path_to_ref,mode='r')
            _rewrite(input_file,genomeRef,chrlist)
            input_file.close()
        else:
            input_file = open(path_to_ref)
            _rewrite(input_file,genomeRef,chrlist)
            input_file.close()
        return genomeRef

    def gtf_to_sql(self,gtf_path,sql_path=None):
        def _split_attributes(stream,field='attributes'):
            fi = stream.fields.index(field)
            for x in stream:
                for y in x[fi].split(';'):
                    header = re.search(r'\s*(\S+) \S*',y)
                    if header:
                        yield header.groups()[0].strip('"')

        def _exon_number(x):
            if x is None: x='1'
            return int(x)

        def _fix_exon_id(stream,count):
            exon_map = {}
            kix = [stream.fields.index(s) for s in ['chr','start','end','strand']]
            for x in stream:
                type = x[2]
                exon_id = x[12]
                if type == "exon" and not exon_id:
                    exkey = tuple(x[i] for i in kix)
                    if exkey not in exon_map:
                        exon_id = "gr_exid_%i"%count
                        count += 1
                        exon_map[exkey] = exon_id
                    x = x[:12]+(exon_map[exkey],)+x[13:]
                yield x

        _intypes = {'exon_number': _exon_number}
        std_outfields = ['gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number']
        sql_params = {'info': {'datatype':'relational'},
                      'outtypes': {'exon_number': 'integer'},
                      'chrmeta': self.chrmeta }
        gtf_read_fields = ['chr','source','name','start','end','strand','frame','attributes']
        sql_fields = ['chr','biotype','type','start','end','strand','frame']+std_outfields
        exon_count = 1

        if sql_path is None: sql_path = unique_filename_in()+".sql"
        if gtf_path is None:
            gtf_path = os.path.join(self.genrep.root,"nr_assemblies/gtf/%s_%i.gtf.gz"%(self.md5,0))
            for n in range(1,100):
                gtf_path_n = os.path.join(self.genrep.root,"nr_assemblies/gtf/%s_%i.gtf.gz"%(self.md5,n))
                if not os.path.exists(gtf_path_n): break
                gtf_path = gtf_path_n
        gtf = track(gtf_path, intypes=_intypes, chrmeta=self.chrmeta)
        if gtf.format == 'sql': return gtf_path
        all_fields = set(_split_attributes(gtf.read(fields=['attributes'])))
        new_fields = [f for f in all_fields if not(f in std_outfields)]
        xsplit = split_field(ensembl_to_ucsc(gtf.read(fields=gtf_read_fields)),
                             outfields=std_outfields+new_fields, infield='attributes',
                             header_split=' ', strip_input=True)
        outf = track(sql_path, fields=sql_fields+new_fields, **sql_params)
        outf.write(FeatureStream(_fix_exon_id(
                    map_chromosomes(xsplit,self.chromosomes),exon_count),
                                 sql_fields[:7]+xsplit.fields[7:]))
        outf.close()
        return os.path.abspath(sql_path)

    def get_sqlite_url(self):
        """Return the url of the sqlite file containing gene annotations."""
        return '%s/data/nr_assemblies/annot_tracks/%s.sql' %(self.genrep.url, self.md5)

    @property
    def sqlite_path(self):
        """Return the path to the sqlite file containing genes annotations."""
        if hasattr(self,"annot_path"): return self.annot_path
        if not hasattr(self,"md5"): return None
        root = os.path.join(self.genrep.root,"nr_assemblies/annot_tracks")
        return os.path.join(root,self.md5+".sql")

    @property
    def annotations_path(self):
        """Return the path an annotation file if available (e.g. microbiome)."""
        if hasattr(self,"annot_path"): return self.annot_path
        if hasattr(self,"md5"):
            root = os.path.join(self.genrep.root,"nr_assemblies/annot_txt")
            f = os.path.join(root,self.md5+".txt")
            if os.path.exists(f): return f
        return None

    def get_features_from_gtf(self,h,chr=None,method="dico"):
        """
        Return a dictionary *data* of the form
        ``{key:[[values],[values],...]}`` containing the result of an SQL request which
        parameters are given as a dictionary *h*. All [values] correspond to a line in the SQL.

        :param chr: (str, or list of str) chromosomes on which to perform the request. By default,
            every chromosome is searched.
        :param method: "dico" or "boundaries": ?

        Available keys for *h*, and possible values:

        * "keys":       "$,$,..."            (fields to `SELECT` and pass as a key of *data*)
        * "values":     "$,$,..."            (fields to `SELECT` and pass as respective values of *data*)
        * "conditions": "$:#,$:#,..."        (filter (SQL `WHERE`))
        * "uniq":       "whatever"           (SQL `DISTINCT` if specified, no matter what the -string- value is)
        * "at_pos":     "12,36,45,1124,..."  (to select only features overlapping this list of positions)

        where

        * $ holds for any column/field name in the database
        * # holds for any value in the database

        Available database fields:

        biotype, type, start, end, strand, frame, gene_id, gene_name, transcript_id, exon_id, exon_number.

        Note: giving several field names to "keys" permits to select unique combinations of these fields.
        The corresponding keys of *data* are a concatenation (by ';') of these fields.
        """
        data = {}
        if isinstance(chr,list): chromosomes = chr
        elif isinstance(chr,str): chromosomes = chr.split(',')
        else: chromosomes = [None]
        if not(method in ["dico","boundaries"]): return data
        # Sqlite3 request to /db/
        dbpath = self.sqlite_path
        if dbpath is None:
            raise ValueError("Annotations not available for this assembly.")
        if os.path.exists(dbpath) and not h.get('at_pos'):
            if chromosomes == [None]: chromosomes = self.chrnames
            db = sqlite3.connect(dbpath)
            cursor = db.cursor()
            dist = h.get('uniq','') and "DISTINCT "
            for chrom in chromosomes:
                if method == 'dico':
                    sql = "SELECT "+dist+",".join(h['keys'].split(',')+h['values'].split(','))
                    sql += " FROM '%s'" %chrom
                else:
                    h['names'] = ",".join([x for x in h['names'].split(',') if not(x == 'chr_name')])
                    sql = "SELECT "+dist+"MIN(start),MAX(end),"+h['names']
                    sql += " FROM '%s'" %chrom
                if h.get('conditions'):
                    sql += " WHERE "
                    conditions = []
                    for c in h['conditions'].split(','):
                        k,v = c.split(':')
                        if "|" in v:
                            conditions.append( "(%s in ('%s'))"%(k,"','".join(v.split("|"))) )
                        else:
                            conditions.append( "(%s='%s')"%(k,v) )
                    sql += " AND ".join(conditions)
                if method == 'boundaries':
                    sql += " GROUP BY "+h['names']
                cursor.execute(sql)
                for x in cursor:
                    if 'names' in h:
                        key = ';'.join([str(_) for _ in x[2:]]+[chrom])
                        values = x[:2]
                    else:
                        nkeys = len(h['keys'].split(','))
                        key = ';'.join([str(_) for _ in x[:nkeys]])
                        values = x[nkeys:]+(chrom,)
                    data.setdefault(key,[]).append(values)
        # Genrep url request
        else:
            for chrom in chromosomes:
                request = self.genrep.url+"/nr_assemblies/get_%s?md5=%s" %(method,self.md5)
                request += "&".join(['']+["%s=%s" %(k,v) for k,v in h.iteritems()])
                if chrom: request += "&chr_name="+chrom
                data.update(json.load(urllib2.urlopen(request)))
                for k,v in data.iteritems():
                    for i,x in enumerate(v):
                        data[k][i] = tuple(x)+(chrom,)
        return data

    def get_gene_mapping(self):
        """
        Return a dictionary ``{gene_id: Gene instance}``
        Note that the gene's length is not the sum of the lengths of its exons.
        """
        gene_mapping = {}
        h = {"keys":"gene_id", "values":"gene_name,start,end,strand,exon_id,transcript_id",
             "conditions":"type:exon", "uniq":"1"}
        for chrom in self.chrnames:
            resp = self.get_features_from_gtf(h,chrom)
            for k,v in resp.iteritems():
                start = min([x[1] for x in v])
                end = max([x[2] for x in v])
                name = str(v[0][0])
                strand = int(v[0][3])
                eid = [str(x[4]) for x in v]
                tid = list(set([str(x[5]) for x in v]))
                gene_mapping[str(k)] = Gene(exons=eid, transcripts=tid,
                    id=k, chrom=chrom, start=start, end=end, name=name, strand=strand)
            # Find gene lengths
            length = fend = 0
            for g,exons in resp.iteritems():
                exons.sort(key=lambda x:(x[1],x[2])) # sort w.r.t. start, then end
                for x in exons:
                    start = x[1]; end = x[2]
                    if start >= fend:
                        length += end-start # new exon
                        fend = end
                    elif end <= fend: pass  # embedded exon
                    else:
                        length += end-fend  # overlapping exon
                        fend = end
                gene_mapping[g].length = length
                length = 0
                fend = 0
        return gene_mapping

    def get_transcript_mapping(self):
        """Return a dictionary ``{transcript_id: Transcript instance}``"""
        transcript_mapping = {}
        h = {"keys":"transcript_id", "values":"gene_id,gene_name,start,end,strand,exon_id",
             "conditions":"type:exon", "uniq":"1"}
        for chrom in self.chrnames:
            resp = self.get_features_from_gtf(h,chrom)
            for k,v in resp.iteritems():
                start = min([x[2] for x in v])
                end = max([x[3] for x in v])
                length = sum([x[3]-x[2] for x in v])
                gid = str(v[0][0])
                gname = str(v[0][1])
                strand = int(v[0][4])
                eid = [str(x[5]) for x in v]
                transcript_mapping[str(k)] = Transcript(gene_id=gid, gene_name=gname, exons=eid,
                    id=k, chrom=chrom, start=start, end=end, strand=strand, length=length)
        return transcript_mapping

    def get_exon_mapping(self):
        """Return a dictionary ``{exon_id: Exon instance}``"""
        exon_mapping = {}
        h = {"keys":"exon_id", "values":"transcript_id,gene_id,gene_name,start,end,strand",
             "conditions":"type:exon", "uniq":"1"}
        for chrom in self.chrnames:
            resp = self.get_features_from_gtf(h,chrom)
            for k,v in resp.iteritems():
                start = int(v[0][3])
                end = int(v[0][4])
                tid = [str(x[0]) for x in v]
                gid = str(v[0][1])
                gname = str(v[0][2])
                strand = int(v[0][5])
                exon_mapping[str(k)] = Exon(gene_id=gid, gene_name=gname, transcripts=tid,
                    id=k, chrom=chrom, start=start, end=end, strand=strand, length=end-start)
        return exon_mapping

    def gene_coordinates(self,id_list):
        """Creates a BED-style stream from a list of gene ids."""
        _fields = ['chr','start','end','name','score','strand']
        _ids = "|".join(id_list)
        sort_list = []
        webh = { "names": "gene_id,gene_name,strand,chr_name",
                 "conditions": "gene_id:%s" %_ids,
                 "uniq":"1" }
        resp = self.get_features_from_gtf(webh,method='boundaries')
        for k,v in resp.iteritems():
            for vv in v:
                start,end = vv
                gene_id,gene_name,strand,chr_name = [str(y).strip() for y in k.split(';')]
                name = "%s|%s" %(gene_id,gene_name)
                sort_list.append((chr_name,start,end,name,0,int(strand)))
        sort_list.sort()
        return FeatureStream(sort_list,fields=_fields)

    def annot_track(self,annot_type='gene',chromlist=None,biotype=["protein_coding"]):
        """
        Return an iterator over all annotations of a given type in the genome.

        :param annot_type: (str) one of 'gene','transcript','exon','CDS'.
        :chrom_list: (list of str) return only features in the specified chromosomes.
        :biotype: (list of str, or None) return only features with the specified biotype(s).
            If None, all biotypes are selected.
        :rtype: track.FeatureStream
        """
        if self.sqlite_path is None:
            raise ValueError("Annotations not available for this assembly.")
        if chromlist is None: chromlist = self.chrnames
        elif isinstance(chromlist,basestring): chromlist = [chromlist]
        _fields = ['chr','start','end','name','strand']
        nmax = 4
        if annot_type == 'gene':
            flist = "gene_id,gene_name,strand"
            webh = { "keys": "gene_id",
                     "values": "start,end,"+flist,
                     "conditions": "type:exon",
                     "uniq":"1" }
        elif annot_type in ['CDS','exon']:
            flist = "exon_id,gene_id,gene_name,strand,frame"
            webh = { "keys": "exon_id,start,end",
                     "values": "start,end,"+flist,
                     "conditions": "type:"+annot_type,
                     "uniq":"1" }
            nmax = 5
            _fields += ['frame']
        elif annot_type == 'transcript':
            flist = "transcript_id,gene_name,strand"
            webh = { "keys": "transcript_id",
                     "values": "start,end,"+flist,
                     "conditions": "type:exon",
                     "uniq":"1" }
        else:
            raise TypeError("Annotation track type '%s' not implemented." %annot_type)
        def _query():
            for chrom in chromlist:
                sort_list = []
                if biotype is not None:
                    for bt in biotype:
                        wh = webh.copy()
                        wh["conditions"]+=",biotype:"+bt
                        resp = self.get_features_from_gtf(wh,chrom)
                        for k,v in resp.iteritems():
                            start = min([x[0] for x in v])
                            end = max([x[1] for x in v])
                            name = "|".join([str(y) for y in v[0][2:nmax]])
                            sort_list.append((start,end,name)+tuple(v[0][nmax:]))
                else:
                    wh = webh
                    if not os.path.exists(self.sqlite_path):
                        raise ValueError("Unavailable local database: try with `biotype=None`.")
                    resp = self.get_features_from_gtf(wh,chrom)
                    for k,v in resp.iteritems():
                        start = min([x[0] for x in v])
                        end = max([x[1] for x in v])
                        name = "|".join([str(y) for y in v[0][2:nmax]])
                        sort_list.append((start,end,name)+tuple(v[0][nmax:]))
                sort_list.sort()
                for k in sort_list: yield (chrom,)+k[:-1]
        return FeatureStream(_query(),fields=_fields)

    def gene_track(self,chromlist=None,biotype=["protein_coding"]):
        """Return a FeatureStream over all protein coding genes annotation in the genome:
        ('chr', start, end, 'gene_id|gene_name', strand)."""
        return self.annot_track(annot_type='gene',chromlist=chromlist,biotype=biotype)

    def exon_track(self,chromlist=None,biotype=["protein_coding"]):
        """Return a FeatureStream over all coding exons annotation in the genome:
        ('chr', start, end, 'exon_id|gene_id|gene_name', strand, phase)."""
        return self.annot_track(annot_type='exon',chromlist=chromlist,biotype=biotype)

    def transcript_track(self,chromlist=None,biotype=["protein_coding"]):
        """Return a FeatureStream over all protein coding transcripts annotation in the genome:
        ('chr', start, end, 'transcript_id|gene_name', strand)."""
        return self.annot_track(annot_type='transcript',chromlist=chromlist,biotype=biotype)

    @property
    def fasta_by_chrom(self):
        """Returns a dictionary of single chromosome fasta files."""
        return getattr(self,"_fasta_by_chrom",
                       dict((v['name'],self.fasta_path(chromosome=k))
                            for k,v in self.chromosomes.iteritems()))

    @property
    def chrmeta(self):
        """
        Return a dictionary of chromosome meta data of the type
        ``{'chr1': {'length': 249250621},'chr2': {'length': 135534747},'chr3': {'length': 135006516}}``
        """
        return dict([(v['name'], dict([('length', v['length']),
                                       ('ac',str(k[0])+"_"+str(k[1])+"."+str(k[2]) if k[1] else str(k[0]))]))
                     for k,v in self.chromosomes.iteritems()])

    @property
    def chrnames(self):
        """Return a list of chromosome names."""
        namelist = [(v.get('num',v['name']),v['name']) for v in self.chromosomes.values()]
        return [x[1] for x in sorted(namelist)]

################################################################################
class GenRep(object):
    """
    Create an object to query a GenRep repository.

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
    def __init__(self, url=None, root=None, config=None, section='genrep'):
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
        """Backward compatibility"""
        return Assembly( assembly=assembly, genrep=self, intype=intype )

    def is_up(self):
        """Check if genrep webservice is available"""
        try:
            urllib2.urlopen(self.url + "/nr_assemblies.json", timeout=2)
        except urllib2.URLError:
            return False
        return True

    def motifs_available(self, genome_id=None):
        """
        List motifs available in genrep, returns a list like (first number is genome id)::

            [('6 ABF1', 'Saccharomyces cerevisiae S288c - ABF1'),
             ('6 ABF2', 'Saccharomyces cerevisiae S288c - ABF2'),
             ('6 ACE2', 'Saccharomyces cerevisiae S288c - ACE2'), ...]

        """
        request = urllib2.Request(self.url + "/genomes.json")
        motif_list = []
        for g in json.load(urllib2.urlopen(request)):
            species = str(g['genome'].get('name')).strip()
            gid = g['genome'].get('id')
            if genome_id and gid != genome_id: continue
            source_url = g['genome'].get('motif_matrix_url')
            if not source_url: continue
            request2 = urllib2.Request("%s/genomes/%i/get_matrix.json" %(self.url, gid))
            for m in json.load(urllib2.urlopen(request2)):
                mname = m['motif']['name']
                motif_list.append(("%i %s"%(gid,mname), species+" - "+mname))
        motif_list.sort(key=lambda x: x[1])
        return motif_list

    def get_motif_PWM(self, genome_id, motif_name, output=None):
        """
        Retieves a motif PWM from its genome_id and name, and saves in the file named as output if not None.
        """
        request = urllib2.Request("%s/genomes/%i/get_matrix.json?gene_name=%s" %(self.url,genome_id,motif_name))
        motif = json.load(urllib2.urlopen(request))[0]['motif']
        if output is None: return motif
        with open(output,"w") as pwmfile:
            t = dict((alpha,n) for n,alpha in enumerate(motif['alphabet']))
            pwm = "\n".join([" ".join(["1"]+[str(x[t[y]]) for y in 'ACGT']) for x in motif['motif']])
            pwmfile.write(pwm)
        return output


    def assemblies_available(self, assembly=None, filter_valid=True):
        """
        Returns a list of tuples (assembly_name, species) available on genrep, or tells if an
        assembly with name *assembly* is available.
        """
        request = urllib2.Request(self.url + "/genomes.json")
        genome_list = {}
        assembly_list = {}
        try:
            for g in json.load(urllib2.urlopen(request)):
                species = str(g['genome'].get('name')).strip()
                gid = g['genome'].get('id')
                if species and gid:
                    genome_list[gid] = species
            request = urllib2.Request(self.url + "/assemblies.json")
            for a in json.load(urllib2.urlopen(request)):
                name = str(a['assembly'].get('name'))
                if name == None: continue
                if filter_valid and not a['assembly'].get('bbcf_valid'): continue
                species = genome_list.get(a['assembly'].get('genome_id'))
                info = "%s (%s)" %(species,name)
                if name == assembly: return info
                if species not in assembly_list: assembly_list[species] = []
                assembly_list[species].append((name,info))
            if assembly == None:
                return [x for k in sorted(assembly_list.keys()) for x in assembly_list[k]]
        except urllib2.URLError:
            return []

    def get_sequence(self, chr_id, coord_list, path_to_ref=None, chr_name=None, ex=None):
        """Parse a slice request to the repository.

        :param chr_id: tuple of the type ``(3066, u'NC_003279', 6)`` (keys of Assembly.chromosomes).
        :param coord_list: (list of (int,int)) sequences' (start,end) coordinates.
        :param path_to_ref: (str) path to a fasta file containing the whole reference sequence.
        :param ex: an optional bein execution to use the `sam_faidx` program.

        Fasta headers are assumed to be of the form ">3066_NC_003279.6 (...)".
        If *path_to_ref* is given and the header is different, give any random value to *chr_id*
        and set *chr_name* to be the fasta header. E.g. *chr_name='chrI'* if the fasta has ">chrI".
        """
        def _read_fasta(path,coord):
            a = 0
            with open(path) as f:
                sequences = []
                cl = iter(coord)
                start,end = cl.next()
                seq = ''
                line = True
                chrom = None
                for i,line in enumerate(f):
                    headpatt = re.search(r'^>(\S+)\s',line)
                    line = line.strip(' \t\r\n')
                    if headpatt:
                        chrom = headpatt.groups()[0]
                    if line.startswith('>'): continue
                    if not chrom in (chr_id[0], chr_name): continue
                    b = a+len(line)
                    while start <= b:
                        start = max(a,start)
                        if end <= b:
                            seq += line[start-a:end-a]
                            sequences.append(seq.upper())
                            seq = ''
                            try:
                                start,end = cl.next()
                            except StopIteration:
                                f.close()
                                return sequences
                        elif end > b:
                            seq += line[start-a:]
                            a = b
                            break
                    a = b
            return sequences

        if len(coord_list) == 0:
            return []
        for k,c in enumerate(coord_list):
            if c[1] < c[0]: coord_list[k] = (c[1],c[0]) # end < start
            elif c[1] == c[0]: coord_list.pop(k)        # end = start
        coord_list = sorted(coord_list)
        if isinstance(chr_id,tuple):
            chrname = "%s_%s.%s"%chr_id if chr_id[1] else str(chr_id[0])
        else:
            chrname = chr_name
            chr_id = (chr_name,'',0)
        if path_to_ref and os.path.exists(path_to_ref):
            locus = ["%s:%i-%i"%(chrname,start+1,end) for start,end in coord_list]
            seq_all = []
            try:
                if ex:
                    for _ln in range(0, len(locus), 100):
                        seq_all += sam_faidx(ex,path_to_ref,locus[_ln:_ln+100])
                else:
                    from bein import execution
                    with execution(None) as ex:
                        for _ln in range(0, len(locus), 100):
                            seq_all += sam_faidx(ex,path_to_ref,locus[_ln:_ln+100])
                return seq_all
            except:  # raw text fasta without index
                return _read_fasta(path_to_ref,coord_list)
        else:
            slices  = ','.join([','.join([str(y) for y in x]) for x in coord_list])
            url     = """%s/chromosomes/%i/get_sequence_part?slices=%s""" % (self.url, chr_id[0], slices)
            request = urllib2.Request(url)
            return urllib2.urlopen(request).read().split(',')

    def get_genrep_objects(self, url_tag, info_tag, filters = None, params = None):
        """
        Get a list of GenRep objets.

        :param url_tag: the GenrepObject type (plural)
        :param info_tag: the GenrepObject type (singular)
        :param filters: a dict that is used to filter the response
        :param params: to add some parameters to the query from GenRep.

        Example:

        To get the genomes related to 'Mycobacterium leprae' species::

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

    .. method:: __repr__()
    """
    def __init__(self, info, key):
        self.__dict__.update(info[key])

    def __repr__(self):
        return str(self.__dict__)

class GenomicObject(object):
    def __init__(self, id='',chrom='',start=0,end=0,name='',score=0.0,strand=0,length=0,seq=''):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.length = length
        self.seq = seq  # sequence

class Gene(GenomicObject):
    def __init__(self, exons=[],transcripts=[], **args):
        GenomicObject.__init__(self, **args)
        self.gene_id = self.id           # to make the structure similar to Exon and Transcript
        self.gene_name = self.name       # same
        self.exons = exons               # list of exons contained
        self.transcripts = transcripts   # list of transcripts contained

class Exon(GenomicObject):
    def __init__(self, gene_id='',gene_name='',transcripts=[], **args):
        GenomicObject.__init__(self, **args)
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.transcripts = transcripts   # list of transcripts it is contained in

class Transcript(GenomicObject):
    def __init__(self, gene_id='',gene_name='',exons=[], **args):
        GenomicObject.__init__(self, **args)
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.exons = exons               # list of exons it contains


################################################################################

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
