# Genreral Modules #
import os

# Tracks path #
tracks_path = 'extras/tracks/' 
try:
    tracks_path = os.path.abspath('/'.join(__file__.split('/')[:-1]) + '../../../' + tracks_path) + '/'
except NameError:
    pass

# Tracks collection #
yeast_chr_file = tracks_path + 'chr/yeast.chr' 
track_collections = {
'Validation': {
  '1': {'path':tracks_path+'qual/bed/validation1.bed', 'type':'qualitative', 'fields':Track.qualitative_fields[:4], 'chrfile':yeast_chr_file},
  '2': {'path':tracks_path+'qual/bed/validation2.bed', 'type':'qualitative', 'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  '3': {'path':tracks_path+'qual/bed/validation3.bed', 'type':'qualitative', 'fields':Track.qualitative_fields,     'chrfile':yeast_chr_file},
  },
}

# Modify tracks collection #
for col_name, col in sorted(track_collections.items()):
    for track_name, track in sorted(col.items()):
        # Make the SQL path #
        old_path = track['path']
        old_name = old_path.split('/')[-1]
        new_dir  = '/'.join(old_path.split('/')[:-2]) + '/' + 'sql' + '/'
        new_name = '.'.join(old_name.split('.')[:-1]  + ['sql'])
        track['path_sql'] = new_dir + new_name        
