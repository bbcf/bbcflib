"""
This script will take a genomic track file and convert
the chromosome names from the 'chr4' format to the 'chrIV'
format
"""

import shutil
from bbcflib import track
from bbcflib.track.track_collection import track_collections

old = track_collections['Yeast']['RP genes']['path_sql']
new = 'gdv_rpgenes.sql'
shutil.copyfile(old, new)
with track.load(new) as t:
    t.integer_to_roman()