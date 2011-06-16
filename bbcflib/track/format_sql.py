"""
===================================
Submodule: bbcflib.track.format_sql
===================================

Implementation of the SQL format.
"""

# General modules #
import sqlite3

# Internal modules #
from ..track import *

###########################################################################   
class GenomicFormat(Track):
    special_tables = ['attributes', 'chrNames', 'types']

    def load(self):
        # Prepare the connection #
        self.connection = sqlite3.connect(self.path)
        self.cursor = self.connection.cursor()
        # Test the file #
        try:
            self.datatype
            self.fields
            self.meta_chr    
            self.meta_track
            self.all_tables
            self.chrs_from_names
            self.chrs_from_tables
        except sqlite3.OperationalError:
            raise Exception("The track '" + self.path + "' doesn't seem to have the correct format")
        # Set attributes #
        self.all_chrs = self.chrs_from_tables

    def unload(self, datatype=None, value=None, trackback=None):
        self.make_missing_indexes()
        self.connection.commit()
        self.cursor.close()
        self.connection.close()

    def commit(self):
        self.connection.commit()

    #-----------------------------------------------------------------------------#
    @property
    def datatype(self):
        self.cursor.execute("select value from attributes where key='datatype'")
        return self.cursor.fetchall()[0][0].encode('ascii')

    @property
    def fields(self):
        if self.chrs_from_tables: return self.get_fields(self.chrs_from_tables[0])
        else:                     return []

    @property
    def meta_chr(self): 
        self.cursor.execute("pragma table_info(chrNames)")
        chrkeys = [x[1] for x in self.cursor.fetchall()]
        self.cursor.execute("select * from chrNames")
        return [dict([(k, e[i]) for i, k in enumerate(chrkeys)]) for e in self.cursor.fetchall()]

    @meta_chr.setter
    def meta_chr(self, data):
        if self.readonly: return
        if data == []: self.cursor.execute('delete from chrNames')  
        for x in data: self.cursor.execute('insert into chrNames (' + ','.join(x.keys()) + ') values (' + ','.join(['?' for y in range(len(x.keys()))])+')', tuple([x[k] for k in x.keys()]))

    @property
    def meta_track(self): 
        self.cursor.execute("select key, value from attributes")
        return dict(self.cursor.fetchall())

    @meta_track.setter
    def meta_track(self, data): 
        if self.readonly: return
        if data == {}: self.cursor.execute('delete from attributes') 
        for k in data.keys(): self.cursor.execute('insert into attributes (key,value) values (?,?)',(k,data[k]))

    @property
    def name(self):
        if self._name: return self.name
        if 'name' in self.meta_track: return self.meta_track['name']
        return 'Unamed'

    @name.setter
    def name(self, value): 
        self._name = value

    #-----------------------------------------------------------------------------#
    def read(self, selection=None, fields=None, order='start,end', cursor=False):
        # Default selection #
        if not selection:
            selection = self.chrs_from_tables
        # Case multi-chromosome #
        if type(selection) == list:
            return join_read_queries(self, selection, fields)
        # Case chromsome name #
        if type(selection) == str:
            if selection not in self.chrs_from_tables: return ()
            if not fields: fields = self.get_fields(selection)
            sql_request = "select " + ','.join(fields) + " from '" + selection + "'"
        # Case span dictionary #
        if type(selection) == dict:
            chrom = selection['chr']
            if chrom not in self.chrs_from_tables: return ()
            if not fields: fields = self.get_fields(chrom)
            sql_request = "select " + ','.join(fields) + " from '" + chrom + "' where " + make_cond_from_sel(selection)
        # Ordering #
        order_by = 'order by ' + order
        # Return the results #
        if cursor: cur = self.connection.cursor()
        else:      cur = self.cursor
        return cur.execute(sql_request + ' ' + order_by)

    def write(self, chrom, data, fields=None):
        if self.readonly: return
        # Default fields #
        if self.datatype == 'quantitative': fields = Track.quantitative_fields
        if fields        == None:           fields = Track.qualitative_fields
        # Get types #
        columns = ','.join([field + ' ' + Track.field_types.get(field, 'text') for field in fields])
        # Maybe create the table #
        if chrom not in self.chrs_from_tables:
            self.cursor.execute('create table "' + chrom + '" (' + columns + ')')
        # Execute the insertion
        sql_command = 'insert into "' + chrom + '" values (' + ','.join(['?' for x in range(len(fields))])+')'
        try: 
            self.cursor.executemany(sql_command, data)
        except sqlite3.OperationalError as err:
            raise Exception("The command '" + sql_command + "' on the database '" + self.path + "' failed with error: " + str(err))

    def remove(self, chrom=None):
        if self.readonly: return
        if not chrom:
            chrom = self.chrs_from_tables
        if type(chrom) == list:
            for ch in chrom:
                self.remove(ch)
        else:
            self.cursor.execute("DROP table '" + chrom + "'")

    def count(self, selection=None):
        # Default selection #
        if not selection:
            selection = self.chrs_from_tables
        # Case multi-chromosome #
        if type(selection) == list:
            return sum([self.count(s) for s in selection])
        # Case chromsome name #
        if type(selection) == str:
            if selection not in self.chrs_from_tables: return 0
            sql_request = "select COUNT(*) from '" + selection + "'"
        # Case span dictionary #
        if type(selection) == dict:
            chrom = selection['chr']
            if chrom not in self.chrs_from_tables: return 0
            sql_request = "select COUNT(*) from '" + chrom + "' where " + make_cond_from_sel(selection)
        # Return the results #
        return self.cursor.execute(sql_request).fetchone()[0]

    #-----------------------------------------------------------------------------#
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

    #-----------------------------------------------------------------------------#
    def get_fields(self, chrom):
        return [x[1] for x in self.cursor.execute('pragma table_info("' + chrom + '")').fetchall()]

    def make_missing_indexes(self):
        if self.readonly: return
        try: 
            for chrom in self.chrs_from_tables:
                self.cursor.execute(    "create index IF NOT EXISTS '" + chrom + "_range_idx' on '" + chrom + "' (start,end)")
                if 'score' in self.get_fields(chrom):
                    self.cursor.execute("create index IF NOT EXISTS '" + chrom + "_score_idx' on '" + chrom + "' (score)")
                if 'name' in self.get_fields(chrom):
                    self.cursor.execute("create index IF NOT EXISTS '" + chrom + "_name_idx' on '" +  chrom + "' (name)")
        except sqlite3.OperationalError as err:
            raise Exception("The index creation on the database '" + self.path + "' failed with error: " + str(err))

    #-----------------------------------------------------------------------------#
    @staticmethod
    def create(path, datatype, name):
        connection = sqlite3.connect(path)
        cursor = connection.cursor()
        cursor.execute('create table chrNames (name text, length integer)') 
        cursor.execute('create table attributes (key text, value text)')
        if datatype == 'quantitative':
            cursor.execute('insert into attributes (key,value) values ("datatype","quantitative")') 
        if datatype == 'qualitative':
            cursor.execute('insert into attributes (key,value) values ("datatype","qualitative")') 
        connection.commit()
        cursor.close()
        connection.close()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------# 
