# General Modules #
import os

# Specific module #
from ..track import Track

# Tracks path #
tracks_path = 'extras/tracks/' 
try:
    tracks_path = os.path.abspath('/'.join(os.path.realpath(__file__).split('/')[:-1]) + '../../../' + tracks_path) + '/'
except NameError:
    pass

# Tracks collection #
yeast_chr_file = tracks_path + 'chr/yeast.chr' 
track_collections = {
'Validation': {
  1: {'path':tracks_path+'qual/bed/validation1.bed', 'type':'qualitative', ' fields':Track.qualitative_fields[:4], 'chrfile':yeast_chr_file},
  2: {'path':tracks_path+'qual/bed/validation2.bed', 'type':'qualitative',  'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  3: {'path':tracks_path+'qual/bed/validation3.bed', 'type':'qualitative',  'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  4: {'path':tracks_path+'qual/bed/validation4.bed', 'type':'qualitative',  'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  },
'Scores': {
  1: {'path':tracks_path+'quan/wig/scores1.wig',     'type':'quantitative', 'fields':Track.quantitative_fields,    'chrfile':yeast_chr_file},
  2: {'path':tracks_path+'quan/wig/scores2.wig',     'type':'quantitative', 'fields':Track.quantitative_fields,    'chrfile':yeast_chr_file},
  3: {'path':tracks_path+'quan/wig/scores3.wig',     'type':'quantitative', 'fields':Track.quantitative_fields,    'chrfile':yeast_chr_file},
  4: {'path':tracks_path+'quan/wig/scores4.wig',     'type':'quantitative', 'fields':Track.quantitative_fields,    'chrfile':yeast_chr_file},
    },
'Random': {
  1: {'path':tracks_path+'qual/bed/random1.bed', 'type':'qualitative', ' fields':Track.qualitative_fields},
  2: {'path':tracks_path+'qual/bed/random2.bed', 'type':'qualitative',  'fields':Track.qualitative_fields},
  3: {'path':tracks_path+'qual/bed/random3.bed', 'type':'qualitative',  'fields':Track.qualitative_fields},
  4: {'path':tracks_path+'qual/bed/random4.bed', 'type':'qualitative',  'fields':Track.qualitative_fields},
    },
'Yeast': {
  'All genes':  {'path':tracks_path+'qual/bed/all_yeast_genes.bed',   'type':'quantitative',
                 'fields':Track.qualitative_fields, 'chrfile':yeast_chr_file},
  'Ribi genes': {'path':tracks_path+'qual/bed/ribosome_genesis.bed',  'type':'quantitative',
                 'fields':Track.qualitative_fields, 'chrfile':yeast_chr_file},
  'RP genes':   {'path':tracks_path+'qual/bed/ribosome_proteins.bed', 'type':'quantitative',
                 'fields':Track.qualitative_fields, 'chrfile':yeast_chr_file},
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
parser_tests  = [(path + b, True)  for path in [tracks_path + 'qual/bed/should_pass/'] for b in os.listdir(path) if b.endswith(".bed")]
parser_tests += [(path + b, False) for path in [tracks_path + 'qual/bed/should_fail/'] for b in os.listdir(path) if b.endswith(".bed")]
parser_tests += [(path + w, True)  for path in [tracks_path + 'quan/wig/should_pass/'] for w in os.listdir(path) if w.endswith(".wig")]
parser_tests += [(path + w, False) for path in [tracks_path + 'quan/wig/should_fail/'] for w in os.listdir(path) if w.endswith(".wig")]

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------#
