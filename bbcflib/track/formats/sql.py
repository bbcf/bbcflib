"""
====================================
Submodule: bbcflib.track.formats.sql
====================================

Implementation of the SQL format.
"""

# Built-in modules #
import sqlite3

# Internal modules #
from .. import Track
from ..track_util import join_read_queries, make_cond_from_sel
from ..extras.sql import TrackExtras
from ..common import natural_sort

################################################################################
class TrackFormat(Track, TrackExtras):
    special_tables = ['attributes', 'chrNames', 'types']

    def __init__(self, path, format=None, name=None, chrmeta=None, datatype=None, readonly=False, empty=False):
        # Set the path #
        self.path = path
        # Prepare the connection #
        self.connection = sqlite3.connect(self.path)
        self.cursor = self.connection.cursor()
        self.modified = False
        # Get the track meta data #
        self.attributes = self.attributes_read()
        self.attributes.modified = False
        # Get the chrom meta data #
        if chrmeta: self.chrmeta = chrmeta
        else:
            self.chrmeta = self.chrmeta_read()
            self.chrmeta.modified = False
        # Check for missing attributes #
        if name: self.name = name
        if datatype:
            if not self.datatype: self.datatype = datatype
            if self.datatype != datatype: raise Exception("You cannot change the datatype of the track '" + self.path + "'")
        # Set chromosome list #
        self.all_chrs = self.chrs_from_tables
        self.all_chrs.sort(key=natural_sort)

    def unload(self, datatype=None, value=None, traceback=None):
        if self.attributes.modified: self.attributes_write()
        if self.chrmeta.modified: self.chrmeta_write()
        self.make_missing_indexes()
        self.connection.commit()
        self.cursor.close()
        self.connection.close()

    def commit(self):
        self.connection.commit()

    def make_missing_indexes(self):
        if self.readonly: return
        try:
            for chrom in self.chrs_from_tables:
                self.cursor.execute(    "create index IF NOT EXISTS '" + chrom + "_range_idx' on '" + chrom + "' (start,end)")
                if 'score' in self.get_fields_of_table(chrom):
                    self.cursor.execute("create index IF NOT EXISTS '" + chrom + "_score_idx' on '" + chrom + "' (score)")
                if 'name' in self.get_fields_of_table(chrom):
                    self.cursor.execute("create index IF NOT EXISTS '" + chrom + "_name_idx' on '" +  chrom + "' (name)")
        except sqlite3.OperationalError as err:
            raise Exception("The index creation on the database '" + self.path + "' failed with error: " + str(err))

    #--------------------------------------------------------------------------#
    @property
    def fields(self):
        if self.chrs_from_tables: return self.get_fields_of_table(self.chrs_from_tables[0])
        else:                     return []

    @property
    def all_tables(self):
        self.cursor.execute("select name from sqlite_master where type='table'")
        return [x[0].encode('ascii') for x in self.cursor.fetchall()]

    @property
    def chrs_from_names(self):
        self.cursor.execute("select name from chrNames")
        return [x[0].encode('ascii') for x in self.cursor.fetchall()]

    @property
    def chrs_from_tables(self):
        return [x for x in self.all_tables if x not in self.special_tables and not x.endswith('_idx')]

    def get_fields_of_table(self, table):
        return [x[1] for x in self.cursor.execute('pragma table_info("' + table + '")').fetchall()]

    #--------------------------------------------------------------------------#
    def attributes_read(self):
        if not 'attributes' in self.all_tables: return {}
        self.cursor.execute("select key, value from attributes")
        return dict(self.cursor.fetchall())

    def attributes_write(self):
        if self.readonly: return
        if not 'attributes' in self.all_tables: self.cursor.execute('create table attributes (key text, value text)')
        if not self.attributes: self.cursor.execute('delete from attributes')
        for k in self.attributes.keys(): self.cursor.execute('insert into attributes (key,value) values (?,?)', (k, self.attributes[k]))

    def chrmeta_read(self):
        if not 'chrNames' in self.all_tables: return {}
        self.cursor.execute("pragma table_info(chrNames)")
        column_names = [x[1] for x in self.cursor.fetchall()]
        all_rows = self.cursor.execute("select * from chrNames").fetchall()
        return (column_names, all_rows)

    def chrmeta_write(self):
        if self.readonly: return
        if not 'chrNames' in self.all_tables: self.cursor.execute('create table chrNames (name text, length integer)')
        if not self.chrmeta: self.cursor.execute('delete from chrNames')
        for r in self.chrmeta.rows: self.cursor.execute('insert into chrNames (' + ','.join(r.keys()) + ') values (' + ','.join(['?' for x in r.keys()])+')', tuple(r.values()))

    @property
    def chrmeta(self):
        return self._chrmeta

    @chrmeta.setter
    def chrmeta(self, value):
        self._chrmeta(value)

    @property
    def attributes(self):
        return self._attributes

    @attributes.setter
    def attributes(self, value):
        self._attributes(value)

    @property
    def datatype(self):
        return self.attributes.get('datatype')

    @datatype.setter
    def datatype(self, value):
        if value not in ['quantitative', 'qualitative']:
            raise Exception("The datatype you are trying to use is invalid: '" + str(value) + "'.")
        self.attributes['datatype'] = value

    @property
    def name(self):
        return self.attributes.get('name', 'Unamed')

    @name.setter
    def name(self, value):
        self.attributes['name'] = value

    #--------------------------------------------------------------------------#
    def read(self, selection=None, fields=None, order='start,end', cursor=False):
        # Default selection #
        if not selection:
            selection = self.chrs_from_tables
        # Case list of things #
        if isinstance(selection, list) or isinstance(selection, tuple):
            return join_read_queries(self, selection, fields)
        # Case chromsome name #
        if isinstance(selection, basestring):
            if selection not in self.chrs_from_tables: return ()
            if not fields: fields = self.get_fields_of_table(selection)
            sql_request = "select " + ','.join(fields) + " from '" + selection + "'"
        # Case selection dictionary #
        if isinstance(selection, dict):
            chrom = selection['chr']
            if chrom not in self.chrs_from_tables: return ()
            if not fields: fields = self.get_fields_of_table(chrom)
            sql_request = "select " + ','.join(fields) + " from '" + chrom + "' where " + make_cond_from_sel(selection)
        # Ordering #
        order_by = 'order by ' + order
        # Return the results #
        if cursor: cur = self.connection.cursor()
        else:      cur = self.cursor
        return cur.execute(sql_request + ' ' + order_by)

    def write(self, chrom, data, fields=None):
        self.modified = True
        if self.readonly: return
        # Default fields #
        if self.datatype == 'quantitative': fields = Track.quantitative_fields
        if fields        == None:           fields = Track.qualitative_fields
        # Maybe create the table #
        if chrom not in self.chrs_from_tables:
            columns = ','.join([field + ' ' + Track.field_types.get(field, 'text') for field in fields])
            self.cursor.execute('create table "' + chrom + '" (' + columns + ')')
        # Execute the insertion
        sql_command = 'insert into "' + chrom + '" (' + ','.join(fields) + ') values (' + ','.join(['?' for x in range(len(fields))])+')'
        try:
            self.cursor.executemany(sql_command, data)
        except sqlite3.OperationalError as err:
            raise Exception("The command '" + sql_command + "' on the database '" + self.path + "' failed with error: " + str(err))

    def remove(self, chrom=None):
        self.modified = True
        if self.readonly: return
        if not chrom:
            chrom = self.chrs_from_tables
        if isinstance(chrom, list):
            for ch in chrom:
                self.remove(ch)
        else:
            self.cursor.execute("DROP table '" + chrom + "'")

    def count(self, selection=None):
        # Default selection #
        if not selection:
            selection = self.chrs_from_tables
        # Case multi-chromosome #
        if isinstance(selection, list) or isinstance(selection, tuple):
            return sum([self.count(s) for s in selection])
        # Case chromsome name #
        if isinstance(selection, basestring):
            if selection not in self.chrs_from_tables: return 0
            sql_request = "select COUNT(*) from '" + selection + "'"
        # Case span dictionary #
        if isinstance(selection, dict):
            chrom = selection['chr']
            if chrom not in self.chrs_from_tables: return 0
            sql_request = "select COUNT(*) from '" + chrom + "' where " + make_cond_from_sel(selection)
        # Return the results #
        return self.cursor.execute(sql_request).fetchone()[0]

    def ucsc_to_ensembl(self):
        '''Converts all entries of a track from the UCSC standard to the Ensembl standard.

           ``ucsc_to_ensembl`` returns nothing.
        '''
        for chrom in self.chrs_from_tables: self.cursor.execute("update '" + chrom + "' set start=start+1")

    def ensembl_to_ucsc(self):
        '''Converts all entries of a track from the Ensembl standard to the UCSC standard.

           ``ensembl_to_ucsc`` returns nothing.
        '''
        for chrom in self.chrs_from_tables: self.cursor.execute("update '" + chrom + "' set start=start-1")

    #--------------------------------------------------------------------------#
    @staticmethod
    def create(path):
        connection = sqlite3.connect(path)
        cursor = connection.cursor()
        connection.commit()
        cursor.close()
        connection.close()

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
