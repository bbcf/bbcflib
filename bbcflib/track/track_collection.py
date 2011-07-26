"""
=========================================
Submodule: bbcflib.track.track_collection
=========================================

Creates a large dictionary of different tracks found in the "extras" directory (mainly small validation tracks). This is then used by the unittests.
"""

# Built-in modules #
import os

# Internal modules #
from . import Track

###########################################################################
# Tracks path #
if os.environ.has_key('TRACKSPATH'):
    tracks_path = os.environ['TRACKSPATH']
    if not tracks_path.endswith('/'): tracks_path += '/'
else:
    try:
        tracks_path = 'extras/tracks/'
        tracks_path = os.path.abspath('/'.join(os.path.realpath(__file__).split('/')[:-1]) + '../../../' + tracks_path) + '/'
    except NameError:
        pass

# Tracks collection #
yeast_chr_file = tracks_path + 'chr/yeast.chr'
track_collections = {
'Validation': {
  1:   {'path': tracks_path+'qual/bed/validation1.bed',         'type':'qualitative', ' fields':Track.qualitative_fields[:4]},
  2:   {'path': tracks_path+'qual/bed/validation2.bed',         'type':'qualitative',  'fields':Track.qualitative_fields    },
  3:   {'path': tracks_path+'qual/bed/validation3.bed',         'type':'qualitative',  'fields':Track.qualitative_fields    },
  4:   {'path': tracks_path+'qual/bed/validation4.bed',         'type':'qualitative',  'fields':Track.qualitative_fields    },
  },
'Scores': {
  1:   {'path': tracks_path+'quan/wig/scores1.wig',             'type':'quantitative', 'fields':Track.quantitative_fields   },
  2:   {'path': tracks_path+'quan/wig/scores2.wig',             'type':'quantitative', 'fields':Track.quantitative_fields   },
  3:   {'path': tracks_path+'quan/wig/scores3.wig',             'type':'quantitative', 'fields':Track.quantitative_fields   },
  4:   {'path': tracks_path+'quan/wig/scores4.wig',             'type':'quantitative', 'fields':Track.quantitative_fields   },
    },
'Random': {
  1:   {'path': tracks_path+'qual/bed/random1.bed',             'type':'qualitative',  'fields':Track.qualitative_fields    },
  2:   {'path': tracks_path+'qual/bed/random2.bed',             'type':'qualitative',  'fields':Track.qualitative_fields    },
  3:   {'path': tracks_path+'qual/bed/random3.bed',             'type':'qualitative',  'fields':Track.qualitative_fields    },
  4:   {'path': tracks_path+'qual/bed/random4.bed',             'type':'qualitative',  'fields':Track.qualitative_fields    },
    },
'Signals': {
  1:   {'path': tracks_path+'quan/bedgraph/scores1.bedGraph',   'type':'quantitative', 'fields':Track.quantitative_fields   },
  'A': {'path': tracks_path+'quan/bedgraph/testA.bedGraph',     'type':'quantitative', 'fields':Track.quantitative_fields   },
  'B': {'path': tracks_path+'quan/bedgraph/testB.bedGraph',     'type':'quantitative', 'fields':Track.quantitative_fields   },
    },
'Binary': {
  1:   {'path': tracks_path+'quan/bigWig/scores1.bigWig',       'type':'quantitative', 'fields':Track.quantitative_fields   },
  'A': {'path': tracks_path+'quan/bigWig/testA.bigWig',         'type':'quantitative', 'fields':Track.quantitative_fields   },
  'B': {'path': tracks_path+'quan/bigWig/testB.bigWig',         'type':'quantitative', 'fields':Track.quantitative_fields   },
    },
'Peaks': {
  'Rap1': {'path':tracks_path+'quan/wig/rap1.wig',              'type':'quantitative', 'fields':Track.quantitative_fields   },
  'Pol2': {'path':tracks_path+'quan/wig/pol2.wig',              'type':'quantitative', 'fields':Track.quantitative_fields   },
    },
'Special': {
  'Corrupted': {'path':tracks_path+'quan/sql/corrupted.sql',    'type':'quantitative', 'fields':Track.quantitative_fields   },
    },
'Yeast': {
  'All genes':  {'path':tracks_path+'qual/bed/all_yeast_genes.bed',   'type':'qualitative', 'fields':Track.qualitative_fields},
  'Ribi genes': {'path':tracks_path+'qual/bed/ribosome_genesis.bed',  'type':'qualitative', 'fields':Track.qualitative_fields},
  'RP genes':   {'path':tracks_path+'qual/bed/ribosome_proteins.bed', 'type':'qualitative', 'fields':Track.qualitative_fields},
    },
}

# Modify tracks collection #
for col_name, col in sorted(track_collections.items()):
    for track_num, track in sorted(col.items()):
        # Make the SQL path #
        old_path = track['path']
        old_name = old_path.split('/')[-1]
        new_dir  = '/'.join(old_path.split('/')[:-2]) + '/' + 'sql' + '/'
        new_name = '.'.join(old_name.split('.')[:-1]  + ['sql'])
        track['path_sql'] = new_dir + new_name
        # Make a name #
        track['name'] = col_name + ' ' + str(track_num)

# Extra test tracks #
parser_tests = []
for f in ['bed']:
    parser_tests += [(path + t, True)   for path in [tracks_path + 'qual/' +f+ '/should_pass/'] for t in os.listdir(path) if t.endswith('.'+f)]
    parser_tests += [(path + t, False)  for path in [tracks_path + 'qual/' +f+ '/should_fail/'] for t in os.listdir(path) if t.endswith('.'+f)]
for f in ['wig', 'bedgraph']:
    parser_tests += [(path + t, True)   for path in [tracks_path + 'quan/' +f+ '/should_pass/'] for t in os.listdir(path) if t.endswith('.'+f)]
    parser_tests += [(path + t, False)  for path in [tracks_path + 'quan/' +f+ '/should_fail/'] for t in os.listdir(path) if t.endswith('.'+f)]

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
