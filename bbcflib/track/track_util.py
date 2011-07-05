"""
===================================
Submodule: bbcflib.track.track_util
===================================

Usefull stuff.
"""

# Built-in modules #
import sys, random

# Internal modules #
from . import Track, new

###########################################################################
def determine_format(path):
    '''Try to guess the format of a track given its path. Returns a three letter extension'''
    # List of names to three letter extension #
    known_format_extensions = {
        'SQLite 3.x database':                       'sql',
        'SQLite database (Version 3)':               'sql',
        'Hierarchical Data Format (version 5) data': 'hdf5',
    }
    # Get extension #
    extension = os.path.splitext(path)[1][1:]
    # If no extension found then try magic #
    if not extension:
        try:
            import magic
        except ImportError:
            raise Exception("The format of the track '" + path + "' cannot be determined")
        # Let the user customize magic #
        if os.path.exists('magic'):
            m = magic.Magic('magic')
        else:
            m = magic.Magic()
        # Does the file even exist ? #
        try:
            filetype = m.from_file(path)
        except IOError:
            raise Exception("The format of the path '" + path + "' cannot be determined. Please specify a format or add an extension.")
        # Check the result of magic #
        try:
            extension = known_format_extensions[filetype]
        except KeyError:
            raise Exception("The format of the track '" + path + "' resolves to " + filetype + " which is not supported at the moment.")
    # Synonyms #
    if extension == 'db': extension = 'sql'
    # Return the format #
    return extension

def import_implementation(format):
    '''Try to import the implementation of a given format'''
    try:
        if not hasattr(sys.modules[__package__], format):
            __import__(__package__ + '.' + format)
        return sys.modules[__package__ + '.' + format]
    except (ImportError, AttributeError):
        raise Exception("The format '" + format + "' is not supported at the moment")

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
def random_track(track_path, random_track_path, repeat_number=1):
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
