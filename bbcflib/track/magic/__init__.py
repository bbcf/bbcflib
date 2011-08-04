"""
===============================
Subpackage: bbcflib.track.magic
===============================

This function uses the python-magic project at ftp://ftp.astron.com/pub/file/ to guess formats.
The python-magic is part of file project (file unix command)

"""

# Built-in modules #
from pkg_resources import resource_filename

###############################################################################
def guess_file_format(path):
    # Link between descriptions and file extension #
    known_formats = {
        'SQLite 3.x database':                          'sql',
        'SQLite database (Version 3)':                  'sql',
        'Hierarchical Data Format (version 5) data':    'hdf5',
        'BED Document, ':                               'bed',
        'PSL Document, ':                               'psl',
        'GFF Document, ':                               'gff',
        'GTF Document, ':                               'gtf',
        'WIG Document, ':                               'wig',
        'MAF Document, ':                               'maf',
    }
    # Try import #
    try: import magic
    except ImportError: return ''
    # Try version #
    if not hasattr(magic, 'open'): return ''
    # Add our definitions #
    mime = magic.open(magic.NONE)
    mime.load(file = resource_filename(__name__, 'magic_data'))
    filetype = mime.file(path)
    # Otherwise try standard definitions #
    if not filetype in known_formats:
        mime = magic.open(magic.NONE)
        mime.load()
        filetype = mime.file(path)
    # Try the conversion dict #
    return known_formats.get(filetype, '')
