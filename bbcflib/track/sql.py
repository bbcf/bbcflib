"""
============================
Submodule: bbcflib.track.sql
============================

Implementation of the SQL format.
"""

import sqlite3
from ..track import *

###########################################################################   
class GenomicFormat(Track):
    special_tables = ['attributes', 'chrNames']

    def load(self):
        # Prepare the connection #
        self.connection = sqlite3.connect(self.path)
        self.cursor = self.connection.cursor()
        # Test the file #
        try:
            self.type
            self.fields
            self.meta_chr    
            self.meta_track
            self.all_tables
            self.chrs_from_names
            self.chrs_from_tables
        except sqlite3.OperationalError:
            raise Exception("The sql track '" + self.path + "' doesn't seem to have the correct format")
        # Set attributes #
        self.format   = "sql"
        self.all_chrs = self.chrs_from_tables

    def unload(self, type, value, trackback):
        if not type: self.connection.commit() 
        self.cursor.close()
        self.connection.close()

    #-----------------------------------------------------------------------------#
    @property
    def type(self):
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
        for x in data: self.cursor.execute('insert into chrNames (' + ','.join(x.keys()) + ') values (' + ','.join(['?' for y in range(len(x.keys()))])+')', tuple([x[k] for k in x.keys()]))

    @property
    def meta_track(self): 
        self.cursor.execute("select key, value from attributes")
        return dict(self.cursor.fetchall())

    @meta_track.setter
    def meta_track(self, data): 
        for k in data.keys(): self.cursor.execute('insert into attributes (key,value) values (?,?)',(k,data[k]))

    #-----------------------------------------------------------------------------#
    def read(self, selection=None, fields=None):
        # Default fields #
        if not fields:
            fields = self.fields
        # Default selection #
        if not selection:
            selection = self.chrs_from_tables
        # Special case multi-chromosome #
        if type(selection) == list:
            return join_read_queries(self, selection, fields)
        # Other cases #
        if type(selection) == str:
            if selection not in self.chrs_from_tables: return ()
            sql_request = "select " + ','.join(fields) + " from '" + selection + "'"
        if type(selection) == dict:
            chrom = selection['chr']
            if chrom not in self.chrs_from_tables: return ()
            if selection.get('inclusion') == 'strict':
                condition = "start < " + str(selection['end'])   + " and " + "start >= " + str(selection['start']) + \
                  " and " + "end   > " + str(selection['start']) + " and " + "end <= "   + str(selection['end'])
            else:
                condition = "start < " + str(selection['end']) + " and " + "end > " + str(selection['start'])
            sql_request = "select " + ','.join(fields) + " from '" + chrom + "' where " + condition
        # Return the results #
        order_by = 'order by start,end'
        return self.cursor.execute(sql_request + ' ' + order_by)

    def write(self, chrom, data, fields):
        # Default fields #
        if fields == None:
            fields = self.default_fields
        # Get types #
        columns = ','.join([field + ' ' + Track.field_types[field] for field in fields])
        # Maybe create the table #
        if chrom not in self.chrs_from_tables:
            self.cursor.execute('create table "' + chrom + '" (' + columns + ')')
        # Execute the insertion
        sql_command = 'insert into "' + chrom + '" values (' + ','.join(['?' for x in range(len(fields))])+')' 
        self.cursor.executemany(sql_command, data)
    
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
        return [x for x in self.all_tables if x not in self.special_tables]

    #-----------------------------------------------------------------------------#
    def get_fields(self, chrom):
        return [x[1] for x in self.cursor.execute('pragma table_info("' + chrom + '")').fetchall()]

    def make_indexes(self):
        for chrom in self.chrs_from_tables:
            self.cursor.execute(    "create index '" + chrom + "_range_idx' on '" + chrom + "' (start,end)")
            if 'score' in self.fields:
                self.cursor.execute("create index '" + chrom + "_score_idx' on '" + chrom + "' (score)")
            if 'name' in self.fields:
                self.cursor.execute("create index '" + chrom + "_name_idx' on '" +  chrom + "' (name)")

###########################################################################   
def create(path, type, name):
    connection = sqlite3.connect(path)
    cursor = connection.cursor()
    cursor.execute('create table chrNames (name text, length integer)') 
    cursor.execute('create table attributes (key text, value text)')
    if type == 'quantitative':
        cursor.execute('insert into attributes (key,value) values ("datatype","quantitative")') 
    if type == 'qualitative':
        cursor.execute('insert into attributes (key,value) values ("datatype","qualitative")') 
    connection.commit()
    cursor.close()
    connection.close()

#-----------------------------------------#
# This code was written by Lucas Sinclair #
# lucas.sinclair@epfl.ch                  #
#-----------------------------------------# 
