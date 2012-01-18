# Built-in modules #
import os, re, json, shutil, gzip, tarfile, pickle, urllib

# Internal modules #
from . import frontend, genrep, daflims
from bbcflib.common import get_files, cat, set_file_descr, merge_sql, gzipfile, unique_filename_in
from .track import Track, new

# Other modules #
from bein import program, ProgramFailed, MiniLIMS
from bein.util import add_pickle, touch, split_file, count_lines


@program
def pileup(ex,job,minilims=None,hts_url=None,via='lsf'):
 
