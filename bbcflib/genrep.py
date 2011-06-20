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
``/scratch/frt/yearly/genrep/nr_assemblies/bowtie``.  To connect to
this GenRep repository and fetch an ``Assembly`` named ``ce6``, we
would write::

    g = GenRep('http://bbcftools.vital-it.ch/genrep/',
               '/scratch/frt/yearly/genrep/nr_assemblies/bowtie')
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

.. autoclass:: GenRep

.. autoclass:: Assembly
"""
import sqlite3
import urllib2
import json
import os
from ConfigParser   import ConfigParser
from datetime       import datetime
from decimal        import Decimal, getcontext

from common import normalize_url

class GenRep(object):
    """Create an object to query a GenRep repository.

    GenRep is the in-house repository for sequence assemblies for the
    BBCF in Lausanne.  This is an object that wraps its use in Python
    in an idiomatic way.

    Create a GenRep object with the base URL to the GenRep system, and
    the root path of GenRep's files.  For instance,

        g = GenRep('genrep.epfl.ch','/path/to/genrep/indices')

    To get an assembly from the repository, call the assembly
    method with either the integer assembly ID or the string assembly
    name.  This returns an Assembly object.

        a = g.assembly(3)
        b = g.assembly('mus')
    """
    def __init__(self, url=None, root=None, config=None, section='genrep'):
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

    def query_url(self, method, assembly):
        """Assemble a URL to call *method* for *assembly* on the repository."""
        if isinstance(assembly, str):
            return """%s/%s.json?assembly_name=%s""" % (self.url, method, assembly)
        elif isinstance(assembly, int):
            return """%s/%s.json?assembly_id=%d""" % (self.url, method, assembly)
        else:
            raise ValueError("Argument 'assembly' to must be a " + \
                                 "string or integer, got " + str(assembly))

    def get_chromosome(self, chromosome):
        """
        Give a chromosome id or directly a dictionary from json file (genrep format)
        return a Chromosome object
        """
        chromosomeInfo  = None

        if isinstance(chromosome,int): # is an id
            chromosomeInfo = json.load(urllib2.urlopen("""%s/chromosomes/%d.json""" % (self.url, chromosome)))
        elif isinstance(chromosome,dict):
            chromosomeInfo = chromosome
        else:# is not a tuple (json)
            TypeError(u"chromosome type must be an integer for id or a dict from json! Type submited was: "+unicode(type(chromosome)))

        return chromosomeInfo

    def get_chromosome_sequence(self, chromosome):
        """
        This methode take an integer (id) or a Chromosome object
        return full chromosome sequence
        """
        ## @warning this method could be too overload the server get instead the compressed fasta file from assembly_to_compressed_fasta method
        chromosomeId    = None
        if isinstance(chromosome, Chromosome):  # is a Chromosome object
            chromosomeId = chromosome.id
        elif not isinstance(chromosome, int):   # is an integer
            chromosomeId    = chromosome
            chromosome      = self.get_chromosome(chromosomeId)
        else:                                   # is not an integer or a Chromosome object
            raise TypeError(u"chromosome type must be an integer or a Chromosome object! Type submited was: "+unicode(type(chromosome)))
        return self.get_sequence( chromosomeId, [[0, chromosome.length -1]] )

    def get_sequence(self, chr_id, coord_list):
        """Parses a slice request to the repository."""
        if len(coord_list) == 0:
            return []
        slices = ",".join([",".join([str(y) for y in x]) for x in coord_list])
        url = """%s/chromosomes/%i/get_sequence_part?slices=%s""" % (self.url, chr_id, slices)
        return urllib2.urlopen(url).read().split(',')

    def fasta_from_bed(self, chromosomes, out=None, bed=None, sql=None, chunk=50000):
        """Get a fasta file with sequences corresponding to the features in the
        bed or sqlite file.

        Returns the name of the file and the total sequence size.
        """
        def push_slices(slices,start,end,name,cur_chunk):
            if end>start:
                slices['coord'].append([start,end])
                slices['names'].append(name)
                cur_chunk += end-start
            return slices,cur_chunk
        def flush_slices(slices,cid,out):
            names=slices['names']
            coord=slices['coord']
            chr = chr_names[cid]
            with open(out,"a") as f:
                for i,s in enumerate(self.get_sequence(cid,coord)):
                    f.write(">"+names[i]+" "+chr+":"+str(coord[i][0])+"-"+str(coord[i][1])+"\n"+s+"\n")
            return {'coord':[],'names':[]}
        def set_feature_name(features_names, feature_name):
            num             = 0
            isSearchingName = True
            if feature_name in features_names:
                while isSearchingName:
                    num +=1
                    tmp_name = feature_name+"_"+str(num)
                    if tmp_name not in features_names:
                        isSearchingName = False
                        feature_name    = tmp_name
            features_names.add(feature_name)
            return features_names, feature_name

        slices      = {'coord':[],'names':[]}
        chr_names   = dict((c[0],cn['name']) for c,cn in chromosomes.iteritems())
        chr_len     = dict((c[0],cn['length']) for c,cn in chromosomes.iteritems())
        size        = 0
        if bed != None:
            if out == None:
                out = bed+".fa"
            chr_ids = dict((cn['name'],c[0]) for c,cn in chromosomes.iteritems())
            cur_chunk = 0
            cur_chr = 0
            with open(bed,"r") as f:
                for l in f:
                    row             = l.rstrip('\n').split('\t')
                    cid             = chr_ids[row[0]]
                    s               = max(int(row[1]),0)
                    e               = min(int(row[2]),chr_len[cid])
                    name            = ''
                    features_names  = set()
                    if len(row)>3:
                        name = row[3]
                    else:
                        name = "feature"
                    features_names, name    = set_feature_name(features_names, name)
                    slices,cur_chunk= push_slices(slices,s,e,name,cur_chunk)
                    if (cur_chr>0 and cur_chr != cid) or cur_chunk > chunk:
                        size += cur_chunk
                        slices = flush_slices(slices,cur_chr,out)
                        cur_chunk = 0
                    cur_chr = cid
            size += cur_chunk
            slices = flush_slices(slices,cur_chr,out)
        elif sql != None:
            if out == None:
                out = sql+".fa"
            cur_chunk       = 0
            connection      = sqlite3.connect( sql )
            cur             = connection.cursor()
            features_names  = set()
            for k in chromosomes.keys():
                cur.execute('select start,end,name from "'+chr_names[k]+'"')
                connection.commit()
                for row in cur:
                    s               = max(row[0],0)
                    e               = min(row[1],chr_len[cid])
                    features_names, name            = set_feature_name(features_names, row[2])
                    slices,cur_chunk= push_slices(slices,s,e,name,cur_chunk)
                    if cur_chunk > chunk:
                        size        += cur_chunk
                        slices      = flush_slices(slices,k,f)
                        cur_chunk   = 0
                size += cur_chunk
                slices = flush_slices(slices,k,f)
                cur.close()
        else:
            raise TypeError("fasta_from_bed requires either a 'sqlite' or a 'bed' file.")
        return (out,size)

    def assembly(self, assembly):
        """Get an Assembly object corresponding to *assembly*.

        *assembly* may be an integer giving the assembly ID, or a
        string giving the assembly name.
        """
        if isinstance(assembly, str):
            assembly_info = json.load(urllib2.urlopen(self.query_url('assemblies', assembly)))[0]
        elif isinstance(assembly, int):
            assembly_info = json.load(urllib2.urlopen("""%s/assemblies/%d.json""" % (self.url, assembly)))
        else:
            raise ValueError("Argument 'assembly' must be a string or integer, got " + str(assembly))

        a = Assembly(assembly_id = int(assembly_info['assembly']['id']),
                     assembly_name = assembly_info['assembly']['name'],
                     index_path = os.path.join(self.root,
                                               str(assembly_info['assembly']['md5'])),
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

    def assembly_to_compressed_fasta(self, assembly):
        """
        Fields:
        - assembly could be an assembly object or an assembly id or r an assembly name
        return the path to compressed fasta file (tar.gz)
        """
        a       = None
        path    = None

        try:
            if isinstance(assembly, Assembly):
                a = assembly
            else:
                a = self.assembly(assembly)
        except TypeError:
            raise TypeError(u"Assembly type must be an integer for id or an tring for name or a Assembly object! Type submited was: "+unicode(type(assembly)))

        url         = self.url + "/" + "data/nr_assemblies/fasta/" + a.md5 + ".tar.gz"
        with urllib2.urlopen(url) as webFile:
            tempfilename = tempfile.NamedTemporaryFile(suffix=".tar.gz")
            with open(tempfilename.name, 'w') as tempfile:
                tempfile.write(webFile.read())
            path = tempfilename.name
        return path

    def assembly_statistic(self, assembly):
        """
        Return statistic about an assembly
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
        """
        a               = None
        genrepObject    = "nr_assemblies"
        query           = None
        genrepFormat    = "json"
        parameter       = "data_type=counts"

        try:
            if isinstance(assembly, Assembly):
                a = assembly
            else:
                a = self.assembly(assembly)
        except TypeError:
            raise TypeError(u"Assembly type must be an integer for id or an tring for name or a Assembly object! Type submited was: "+unicode(type(assembly)))

        return json.load(urllib2.urlopen("""%s/nr_assemblies/%d.json?data_type=counts""" % (self.url, a.nr_assembly_id)))

    def assembly_statistic_to_file(self, assembly, output):
        getcontext().prec   = 15
        statistic           = self.assembly_statistic( assembly )
        total               = Decimal(statistic["A"] + statistic["T"] + statistic["G"] + statistic["C"])
        name                = (isinstance(assembly, Assembly)) and assembly.name or assembly
        with open(os.path.expanduser(output), "w") as f:
                f.write(u">Assembly: %s\n" %name)
                f.write(u"1\t%s\t%s\t%s\t%s" % ( statistic["A"] / total, statistic["T"] / total,statistic["G"] / total, statistic["C"] / total ) )

    def assemblies_available(self):
        """
        Return list of assemblies available on genrep
        """
        assembly_info   = json.load(urllib2.urlopen(self.url + "/assemblies.json"))
        assembly_list   = [self._assembly(a) for a in assembly_info]
        return [ a for a in assembly_list if a is not None ]


class Assembly(object):
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
    def __init__(self, assembly_id, assembly_name, index_path,
                 bbcf_valid, updated_at, nr_assembly_id, genome_id,
                 source_name, md5, source_id, created_at):
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
