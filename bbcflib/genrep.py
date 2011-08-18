"""
======================
Module: bbcflib.genrep
======================

This module provides an interface to GenRep repositories.  It
provides two classes. ``GenRep`` connects to a GenRep repository and
handles all queries.  A query via the ``GenRep`` object returns an
``Assembly``, giving the information on a particular entry in GenRep.

The primary GenRep repository on VITAL-IT has URL
``http://bbcftools.vital-it.ch/genrep/`` and root directory
``/db/genrep``.  To connect to
this GenRep repository and fetch an ``Assembly`` named ``ce6``, we
would write::

    g = GenRep(url='http://bbcftools.vital-it.ch/genrep/',root='/db/genrep')
    g.assembly('ce6')

Assemblies in GenRep are also assigned unique integer IDs.  The unique
integer ID for assembly ``ce6`` is 14.  We can use these IDs anywhere
we would use the name, so the third line in the prevous code could
equally well be written::

    g.assembly(14)

``GenRep`` objects can also be created from ``ConfigParser`` objects.
Instead of a URL and root directory, we pass keyword arguments
``config`` with the ``ConfigParser`` and optionally ``section`` to
choose what section of the configuration file to use.  If no
``section`` is specified, it defaults to "genrep".  The two fields
read from that section are

  * genrep_url
  * genrep_root

With a ``ConfigParser``, the previous code would look like::

    c = ConfigParser(...)
    ... fill the ConfigParser
    g = GenRep(config=c) # or GenRep(config=c,section='genrep')
    g.assembly('ce6')
"""

# Built-in modules #
import urllib2, json, os, re
from datetime import datetime
# Internal modules #
from .common import normalize_url
# Other modules #
from bein import unique_filename_in

################################################################################
class GenRep(object):
    def __init__(self, url='http://bbcftools.vital-it.ch/genrep/', root='/db/genrep', intype=0, config=None, section='genrep'):
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
            b = g.assembly('mus')

        """
        if (url == None or root == None) and config == None:
            raise TypeError("GenRep requires either a 'url' and 'root', or a 'config'")
        elif config != None:
            self.root = os.path.abspath(config.get(section, 'genrep_root'))
            self.url = normalize_url(config.get(section, 'genrep_url'))
            if url != None:
                self.url = normalize_url(url)
            if root != None:
                self.root = os.path.abspath(root)
        else:
            self.url = normalize_url(url)
            self.root = os.path.abspath(root)
        self.intype = intype

    def is_up(self):
        """ Check if genrep webservice is available """
        try:
            urllib2.urlopen(self.url + "/nr_assemblies.json", timeout=2)
        except urllib2.URLError:
            return False
        return True

    def query_url(self, method, assembly):
        """Assemble a URL to call *method* for *assembly* on the repository."""
        if isinstance(assembly, basestring):
            return urllib2.Request("""%s/%s.json?assembly_name=%s""" % (self.url, method, assembly))
        elif isinstance(assembly, int):
            return urllib2.Request("""%s/%s.json?assembly_id=%d""" % (self.url, method, assembly))
        else:
            raise ValueError("Argument 'assembly' to must be a " + \
                                 "string or integer, got " + str(assembly))

    def get_sequence(self, chr_id, coord_list):
        """Parses a slice request to the repository."""
        if len(coord_list) == 0:
            return []
        slices  = ",".join([",".join([str(y) for y in x]) for x in coord_list])
        url     = """%s/chromosomes/%i/get_sequence_part?slices=%s""" % (self.url, chr_id, slices)
        request = urllib2.Request(url)
        return urllib2.urlopen(request).read().split(',')

    def get_sequence_straight(self, input_list):
        """ each element in input_list are lists containing: the chromosome name, the start and the end position."""    
        h_data = {}
        h_result = {}
        
        for e in input_list:
            if (h_data.get(e[0]) == None):
                h_data[e[0]]=[] 
            h_data[e[0]].append([e[1],e[2]])
        for chr in h_data.keys():
            h_data[chr].sort()
            res_list=get_sequence(chr, h_data[chr])
            for i in range(0,len(res_list)):
                h_results[chr]=[h_data[chr][i], res_list[i]]
        return h_results

    def fasta_from_data(self, chromosomes, data_path, out=None, chunk=50000):
        """Get a fasta file with sequences corresponding to the features in the
        bed or sqlite file.

        Returns the name of the file and the total sequence size.
        """
        from .track import load
        if out == None:
            out = unique_filename_in()
        def push_slices(slices,start,end,name,cur_chunk):
            if end>start:
                slices['coord'].append([start,end])
                slices['names'].append(name)
                cur_chunk += end-start
            return slices,cur_chunk
        def flush_slices(slices,chrid,chrn,out):
            names = slices['names']
            coord = slices['coord']
            with open(out,"a") as f:
                for i,s in enumerate(self.get_sequence(chrid,coord)):
                    f.write(">"+names[i]+"|"+chrn+":"+str(coord[i][0])+"-"+str(coord[i][1])+"\n"+s+"\n")
            return {'coord':[],'names':[]}
        slices = {'coord':[],'names':[]}
        chrlist = dict((v['name'], {'length': v['length']}) for v in chromosomes.values())
        size = 0
        with load(data_path, chrmeta=chrlist) as t:
            cur_chunk = 0
            for cid,chrom in chromosomes.iteritems():
                for row in t.read(selection=chrom['name'], fields=["start","end","name"]):
                    s = max(row[0],0)
                    e = min(row[1],chrom['length'])
                    name = re.sub('\s+','_',row[2])
                    slices,cur_chunk= push_slices(slices,s,e,name,cur_chunk)
                    if cur_chunk > chunk:
                        size += cur_chunk
                        slices = flush_slices(slices,cid[0],chrom['name'],out)
                        cur_chunk = 0
                size += cur_chunk
                slices = flush_slices(slices,cid[0],chrom['name'],out)
        return (out,size)

    def assembly(self, assembly):
        """Get an Assembly object corresponding to *assembly*.

        *assembly* may be an integer giving the assembly ID, or a
        string giving the assembly name.
        """
        if isinstance(assembly, basestring):
            assembly_info = json.load(urllib2.urlopen(self.query_url('assemblies', assembly)))[0]
        elif isinstance(assembly, int):
            assembly_info = json.load(urllib2.urlopen("""%s/assemblies/%d.json""" % (self.url, assembly)))
        else:
            raise ValueError("Argument 'assembly' must be a string or integer, got " + str(assembly))

        root = os.path.join(self.root,"nr_assemblies/bowtie")
        if self.intype == 1:
            root = os.path.join(self.root,"nr_assemblies/exons_bowtie")
        a = Assembly(assembly_id = int(assembly_info['assembly']['id']),
                     assembly_name = assembly_info['assembly']['name'],
                     index_path = os.path.join(root,str(assembly_info['assembly']['md5'])),
                     bbcf_valid = assembly_info['assembly']['bbcf_valid'],
                     updated_at = datetime.strptime(assembly_info['assembly']['updated_at'],
                                                    '%Y-%m-%dT%H:%M:%SZ'),
                     nr_assembly_id = int(assembly_info['assembly']['nr_assembly_id']),
                     genome_id = int(assembly_info['assembly']['genome_id']),
                     source_name = assembly_info['assembly']['source_name'],
                     md5 = assembly_info['assembly']['md5'],
                     source_id = int(assembly_info['assembly']['source_id']),
                     created_at = datetime.strptime(assembly_info['assembly']['created_at'],
                                                    '%Y-%m-%dT%H:%M:%SZ'))
        chromosomes = json.load(urllib2.urlopen(self.query_url('chromosomes', assembly)))
        for c in chromosomes:
            name_dictionary = dict([ (x['chr_name']['assembly_id'],
                                      x['chr_name']['value'])
                                     for x in c['chromosome']['chr_names']])
            a.add_chromosome(c['chromosome']['id'],
                             c['chromosome']['refseq_locus'],
                             c['chromosome']['refseq_version'],
                             name_dictionary[a.id],
                             c['chromosome']['length'])
        return a

    def assemblies_available(self):
        """
        Returns a list of assemblies available on genrep
        """
        request         = urllib2.Request(self.url + "/assemblies.json")
        assembly_info   = json.load(urllib2.urlopen(request))
        assembly_list   = [a['assembly']['name'] for a in assembly_info]
        return [ a for a in assembly_list if a is not None ]

    def is_available(self, assembly):
        """
        Returns a list of assemblies available on genrep
        """
        return assembly in self.assemblies_available()

    def statistics(self, assembly, output=None, frequency=False, matrix_format=False):
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

        if matrix_format is True output is for example:
         >Assembly: sacCer2
        1   0.309798640038793   0.308714120881750   0.190593944221299   0.190893294858157
        """
        request = urllib2.Request("%s/nr_assemblies/%d.json?data_type=counts" % (self.url, assembly.nr_assembly_id))
        stat    = json.load(urllib2.urlopen(request))
        total   = float(stat["A"] + stat["T"] + stat["G"] + stat["C"])
        if frequency:
            stat = dict((k,x/total) for k,x in stat.iteritems())
        else:
            stat.update
        if output == None:
            return stat
        else:
            with open(output, "w") as f:
                if matrix_format:
                    f.write(">Assembly: %s\n" % assembly.name)
                    f.write("%s\t%s\t%s\t%s" %(stat["A"],stat["C"],stat["G"],stat["T"]))
                    f.write("\n")
                else:
                    f.write("#Assembly: %s\n" % assembly.name)
                    [f.write("%s\t%s\n" % (x,stat[x])) for x in ["A","C","G","T"]]
                    f.write("#\n")
                    [[f.write("%s\t%s\n" % (x+y,stat[x+y])) for y in ["A","C","G","T"]] for x in ["A","C","G","T"]]
            return output

    def fasta_path(self, assembly, chromosome=None):
        """
        Returns the path to the compressed fasta file, for the whole assembly or for a single chromosome.
        """
        root = os.path.join(self.root,"nr_assemblies/fasta")
        path = os.path.join(root,assembly.md5+".tar.gz")
        if chromosome != None:
            chr_id = str(chromosome[0])+"_"+str(chromosome[1])+"."+str(chromosome[2])
            root = os.path.join(self.root,"chromosomes/fasta")
            path = os.path.join(root,chr_id+".fa.gz")
        elif self.intype == 1:
            root = os.path.join(self.root,"nr_assemblies/exons_fasta")
            path = os.path.join(root,assembly.md5+".fa.gz")
        return path

    def get_genrep_objects(self, url_tag, info_tag, filters = None):
        """
        Get a list of GenRep objets
        ... attribute url_tag : the GenrepObject type (plural)
        ... attribute info_tag : the GenrepObject type (singular)
        Optionals attributes :
        ... attribute filters : a dict that is used to filter the response
        from GenRep.
        Example :
        To get the genomes related to 'Mycobacterium leprae' species.
        First get the species with the right name :
        species = get_genrep_objects('organisms', 'organism', {'species':'Mycobacterium leprae'})[0]
        genomes = get_genrep_objects('genomes', 'genome', {'organism_id':species.id})
        """
        if filters is None:
            filters = {}
        url = """%s/%s.json""" % (self.url, url_tag)
        infos = json.load(urllib2.urlopen(url))
        result = []
        # get objects
        for info in infos:
            obj = GenrepObject(info,info_tag)
            if not filters :
                result.append(obj)
            else : # filter
                for k,v in filters.items():
                    if hasattr(obj,k):
                        if getattr(obj,k)==v:
                            result.append(obj)
        return result

################################################################################
class Assembly(object):
    def __init__(self, assembly_id, assembly_name, index_path,
                 bbcf_valid, updated_at, nr_assembly_id, genome_id,
                 source_name, md5, source_id, created_at):
        """A representation of a GenRep assembly.

        In general, Assembly objects should always be created by calls to
        a GenRep object.

        An Assembly has the following fields:

        .. attribute:: id

        An integer giving the assembly ID in GenRep.

        .. attribute:: name

        A string giving the name of the assembly in GenRep.

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

        All integers.

        .. attribute:: source_name

        .. attribute:: md5

        """
        self.id = int(assembly_id)
        self.name = assembly_name
        self.chromosomes = {}
        self.index_path = os.path.abspath(index_path)
        self.bbcf_valid = bbcf_valid
        self.updated_at = updated_at
        self.nr_assembly_id = nr_assembly_id
        self.genome_id = genome_id
        self.source_name = source_name
        self.md5 = md5
        self.source_id = source_id
        self.created_at = created_at

    def add_chromosome(self, chromosome_id, refseq_locus, refseq_version, name, length):
        self.chromosomes[(chromosome_id, refseq_locus, refseq_version)] = \
            {'name': name, 'length': length}

    @property
    def chrmeta(self):
        """ Returns a dictionary of chromosome meta data looking something like:

            {'chr1': {'length': 249250621},
             'chr2': {'length': 135534747},
             'chr3': {'length': 135006516},

        """
        return dict([(v['name'], dict([('length', v['length'])])) for v in self.chromosomes.values()])


class GenrepObject(object):
    """
    Class wich will reference all different objects used by GenRep
    In general, you should never instanciate GenrepObject directly but
    call a method from the GenRep object.
    """
    def __init__(self, info, key):
        self.__dict__.update(info[key])



#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
