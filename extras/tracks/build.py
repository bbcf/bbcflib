'''
execfile('extras/tracks/build.py')
create_tracks()
create_wigs()
'''

# Genreral Modules #
import os, subprocess

# Specific module #
from bbcflib.track.common import terminal_colors
from bbcflib.track import load
from bbcflib.track.track_collection import track_collections, tracks_path, yeast_chr_file

# Same randomness #
import random
random.seed(0)

# Create tracks #
def create_tracks():
    for col_name, col in sorted(track_collections.items()):
        if col_name in ['Binary']: continue
        for track_num, track in sorted(col.items()):
            print terminal_colors['txtylw'] + "Creating track '" + track['name'] + "'" + terminal_colors['end']
            if os.path.exists(track['path_sql']): os.remove(track['path_sql'])
            if col_name == 'Random':
                with load('/dev/null', 'random', 'Test random track ' + str(track_num)) as t:
                    t.size = 1000.0*(float(track_num)/2.0)
                    t.convert(track['path_sql'])
            else:
                with load(track['path'], chrfile=yeast_chr_file) as t:
                    t.convert(track['path_sql'])

# Special wig tracks #
def create_wigs():
    from bbcflib.track.format_wig import random_track
    # Pol2
    random.seed(0)
    if os.path.exists(tracks_path + 'quan/wig/pol2.wig'): os.remove(tracks_path + 'quan/wig/pol2.wig')
    with open(tracks_path + 'quan/wig/pol2.wig', 'w') as file: file.writelines(random_track('fixed'))
    if os.path.exists(tracks_path + 'quan/sql/pol2.sql'): os.remove(tracks_path + 'quan/sql/pol2.sql')
    with load(tracks_path + 'quan/wig/pol2.wig', chrfile=yeast_chr_file) as t: t.convert(tracks_path + 'quan/sql/pol2.sql')
    # Rap1
    random.seed(0)
    if os.path.exists(tracks_path + 'quan/wig/rap1.wig'): os.remove(tracks_path + 'quan/wig/rap1.wig')
    with open(tracks_path + 'quan/wig/rap1.wig', 'w') as file: file.writelines(random_track('variable'))
    if os.path.exists(tracks_path + 'quan/sql/rap1.sql'): os.remove(tracks_path + 'quan/sql/rap1.sql')
    with load(tracks_path + 'quan/wig/rap1.wig', chrfile=yeast_chr_file) as t: t.convert(tracks_path + 'quan/sql/rap1.sql')

# Special bigWig tracks #
def create_bigwigs():
    for track_num, d in sorted(track_collections['Binary'].items()):
        print terminal_colors['txtylw'] + "Creating track '" + d['name'] + "'" + terminal_colors['end']
        proc = subprocess.Popen(['bedGraphToBigWig', d['from'], yeast_chr_file, d['path']], stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr: raise Exception("The tool bedGraphToBigWig exited with message: " + '"' + stderr.strip('\n') + '"')
