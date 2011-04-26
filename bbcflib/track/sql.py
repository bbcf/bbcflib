"""
============================
Submodule: bbcflib.track.sql
============================

Implementation of the SQL format.
"""

from ..track import Track

class GenomicFormat(Track):
    special_tables = ['attributes', 'chrNames']

    def load():
        # Prepare the connection #
        self.connection = sqlite3.connect(self.location)
        self.cursor = self.connection.cursor()
        # Test the file #
        try:
            self.type
            self.all_tables
            self.chrs_from_names
            self.meta_chr    
            self.meta_track
        except sqlite3.OperationalError:
            raise Exception("The sql track '" + self.path + "' doesn't seem to have the correct format")
        # Set attributes #
        self.format   = "sql"
        self.all_chrs = self.chrs_from_tables

    #-----------------------------------------------------------------------------#
    @property
    def type(self):
        self.cursor.execute("select value from attributes where key='datatype'")
        return self.cursor.fetchall()[0][0].encode('ascii')

    @property
    def fields(self):
        if self.chrs_from_tables: return get_fields(self.chrs_from_tables[0])
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

    def get_fields(self, chr):
        self.fields = [x[1] for x in self.cursor.execute('pragma table_info("' + self.all_chrs[0] + '")').fetchall()]
    
    #-----------------------------------------------------------------------------#
    def read():
        pass

    def write():
        pass 
