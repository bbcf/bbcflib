import urllib2
import json

class GenRep(object):
    """Create an object to query a GenRep repository.

    GenRep is the in-house repository for sequence assemblies for the
    BBCF in Lausanne.  This is an object that wraps its use in Python
    in an idiomatic way.

    Create a GenRep object with the base URL to the GenRep system, and
    the root path of GenRep's files.  For instance,

    >>> g = GenRep('genrep.epfl.ch','/path/to/genrep/indices')

    To get an assembly from the repository, call the get_assembly
    method with either the integer assembly ID or the string assembly
    name.  This returns an Assembly object.

    >>> a = GenRep.get_assembly(3)
    >>> b = GenRep.get_assembly('mus')
    """
    def __init__(self, url, root):
        self.root = os.path.abspath(root)
        self.url = normalize_url(url)

    def query_url(self, method, assembly):
        """Assemble a URL to call *method* for *assembly* on the repository."""
        if isinstance(assembly, int):
            return """%s/%s.json?assembly_name=%s""" % (self.url, method, assembly)
        elif isinstance(assembly, str):
            return """%s/%s.json?assembly_id=%d""" % (self.url, method, assembly)
        else:
            raise ValueError("Argument 'assembly' to index_path must be a " + \
                                 "string or integer, got " + str(assembly))

    def get_assembly(self, assembly):
        """Get an Assembly object corresponding to *assembly*.

        *assembly* may be an integer giving the assembly ID, or a
        string giving the assembly name.
        """
        assembly_info = json.load(urllib2.urlopen(self.query_url('assemblies', assembly)))
        assembly = Assembly(assembly_id = str(assembly_info[0]['assembly']['id']),
                            assembly_name = assembly_info[0]['assembly']['name'],
                            index_path = os.path.join(self.root, 
                                                      str(assembly_info[0]['assembly']['md5'])))
        chromosomes = json.load(urllib2.urlopen(self.query_url('chromosomes', assembly)))
        for c in chromosomes:
            name_dictionary = dict([ (x['chr_name']['assembly_id'],
                                      x['chr_name']['value'])
                                     for x in c['chromosome']['chr_names']])
            assembly.add_chromosome(c['chromosome']['id'],
                                    c['chromosome']['refseq_locus'],
                                    c['chromosome']['refseq_version'],
                                    name_dictionary[assembly.id],
                                    c['chromosome']['length'])
        return assembly
                                    
            


class Assembly(object):
    """A representation of a GenRep assembly.

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

    In general, Assembly objects should always be created by calls to
    a GenRep object.
    """
    def __init__(self, assembly_id, assembly_name, index_path):
        self.id = assembly_id
        self.name = assembly_name
        self.chromosomes = {}
        self.index_path = os.path.abspath(index_path)

    def add_chromosome(self, chromosome_id, refseq_locus, refseq_version, name, length):
        self.chromosomes[(chromosome_id, refseq_locus, refseq_version)] = \
            {'name': name, 'length': length}

        
        
        
def normalize_url(url):
    """Produce a fixed form for an HTTP URL.

    Make sure the URL begins with http:// and does *not* end with a /.

    >>> normalize_url('http://www.google.com')
    'http://www.google.com'
    >>> normalize_url('http://www.google.com/')
    'http://www.google.com'
    >>> normalize_url('www.google.com/')
    'http://www.google.com'
    """
    url = url.lower()
    if not(url.startswith("http://")):
        url = "http://" + url
    if url.endswith("/"):
        url = url[:-1]
    return url

if __name__ == '__main__':
    import doctest
    doctest.testmod()
