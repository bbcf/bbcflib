"""
===================================
Submodule: bbcflib.track.track_util
===================================

Usefull stuff.
"""

# Built-in modules #
import os, sys, random

# Internal modules #
from . import Track, new
from . import formats

###########################################################################
def determine_format(path):
    '''Try to guess the format of a track given its path. Returns a three letter extension'''
    # Try magic first #
    filetype = magic_format(path)
    # Then try the extension #
    if not filetype: filetype = os.path.splitext(path)[1][1:]
    # If still nothing, raise exception #
    if not filetype:
        raise Exception("The format of the path '" + path + "' cannot be determined. Please specify a format or add an extension.")
    # Synonyms #
    if filetype == 'db': filetype = 'sql'
    # Return the format #
    return filetype

def magic_format(path):
    # List of names to three letter extension #
    known_format_extensions = {
        'SQLite 3.x database':                       'sql',
        'SQLite database (Version 3)':               'sql',
        'Hierarchical Data Format (version 5) data': 'hdf5',
    }
    # Try import #
    try:
        import magic
    except ImportError:
        return ''
    # Try usage #
    try:
        mime = magic.Magic(magic.NONE)
    except AttributeError:
        return ''
    # Let the user customize magic #
    if os.path.exists('magic'): mime.load(file='magic')
    else: mime.load()
    # Does the file even exist ? #
    try: filetype = mime.file(path)
    except IOError: return ''
    # Try the conversion dict #
    return known_format_extensions.get(filetype, filetype)

#-----------------------------------------------------------------------------#
def import_implementation(format):
    '''Try to import the implementation of a given format'''
    if not hasattr(sys.modules[__package__].formats, format):
        __import__(    __package__ + '.formats.' + format)
    return sys.modules[__package__ + '.formats.' + format]

#-----------------------------------------------------------------------------#
def join_read_queries(track, selections, fields):
    '''Join read results when selection is a list'''
    def _add_chromsome(sel, data):
        if type(sel) == str: chrom = (sel,)
        else:                chrom = (sel['chr'],)
        for f in data: yield chrom + f
    for sel in selections:
        for f in _add_chromsome(sel, track.read(sel, fields)): yield f

#-----------------------------------------------------------------------------#
def make_cond_from_sel(selection):
    '''Make an SQL condition string from a selection dictionary'''
    query = ""
    if "start" in selection and "end" in selection:
        if selection.get('inclusion') == 'strict':
            query = "start < " + str(selection['end'])   + " and " + "start >= " + str(selection['start']) + \
               " and end   > " + str(selection['start']) + " and " + "end <= "   + str(selection['end'])
        else:
            query = "start < " + str(selection['end'])   + " and " + "end > " + str(selection['start'])
    if "score" in selection:
        if query: query += " and "
        statements  = [i for i in selection["score"].split(" ") if i != ""]
        for i in statements:
            try: number = float(i)
            except ValueError: symbol = i
        query += "score " + locals().get("symbol", "=") + " " + str(number)
    return query

###########################################################################
def shuffle_track(track_path, random_track_path, repeat_number=1):
    with new(random_track_path, "sql", name="random_track") as random_track:
        with Track(track_path, format="sql") as track:
            random_track.meta_chr   = track.meta_chr
            random_track.meta_track = track.meta_track
            number                  = 0
            for i in range(repeat_number):
                features_list = []
                for chromosome in track.meta_chr:
                    data            = track.read(chromosome["name"])
                    for information in data:
                        distance        = information[0] - information[1]
                        random_start    = random.randint(0, chromosome["length"] - distance)
                        random_end      = random_start + distance
                        feature_name    = "random%d" %(number)
                        features_list.append((random_start, random_end, feature_name, "+", None))
                        number          += 1
                    random_track.write(chromosome["name"], features_list)

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
