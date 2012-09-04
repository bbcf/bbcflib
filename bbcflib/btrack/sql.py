from bbcflib.btrack import *
import sqlite3, warnings

_sql_types = {'start':        'integer',
              'end':          'integer',
              'score':        'real',
              'strand':       'integer',
              'thick_start':  'integer',
              'thick_end':    'integer',
              'block_count':  'integer',
              'frame':        'integer'}

class SqlTrack(Track):
    """
    Track class for sqlite3 files (extension ".sql" or ".db").
    Additional attributes:

    .. attribute:: readonly

        If True, tables will not be updated to reflect, e.g. the chrmeta or info attributes.

    .. attribute:: connection

       The sqlite3 file connection.

    .. attribute:: cursor

       The sqlite3 connection cursor.

    .. attribute:: types

       The field types as defined in the sqlite3 tables.

    """
    def __init__(self,path,**kwargs):
        self.readonly = kwargs.get('readonly',False)
        self.connection = sqlite3.connect(path)
#        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()
        kwargs['format'] = 'sql'
        Track.__init__(self,path,**kwargs)
        self.fields = self._get_fields(fields=self.fields)
        self.types = dict((k,v) for k,v in _sql_types.iteritems() if k in self.fields)
        if isinstance(kwargs.get('types'),dict): self.types.update(kwargs["types"])

    def open(self):
        self._prepare_db()

    def close(self):
        self.connection.commit()
        self.cursor.close()

    @property
    def tables(self):
        """Returns the complete list of SQL tables."""
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        return [x[0].encode('ascii') for x in self.cursor.fetchall()]

###### attributes #######
    def _fix_attributes(self,info=None):
        if info is None:
            info = self.info
        try:
            sql_command = "CREATE TABLE IF NOT EXISTS 'attributes' ('key' TEXT, 'value' TEXT)"
            self.cursor.execute(sql_command)
            self.cursor.execute("SELECT key,value FROM 'attributes'")
            table_attr = dict((x[0],x[1]) for x in self.cursor.fetchall())
            for k,v in info.iteritems():
                if k in table_attr:
                    if v != table_attr[k]:
                        sql_command = "UPDATE 'attributes' SET value='%s' WHERE key='%s'"%(v,k)
                        self.cursor.execute(sql_command)
                else:
                    sql_command = "INSERT INTO 'attributes' (key,value) VALUES ('%s','%s')"%(k,v)
                    self.cursor.execute(sql_command)
        except sqlite3.OperationalError as err:
            raise Exception("Sql error: %s\n on file %s, with\n%s"%(err,self.path,sql_command))
        return True

###### chrNames #######
    def _fix_chrmeta(self,chrmeta=None):
        if chrmeta is None: chrmeta = self.chrmeta
#        chrmeta = self._get_chrmeta( chrmeta )
        try:
            sql_command = "CREATE TABLE IF NOT EXISTS 'chrNames' ('name' TEXT, 'length' INTEGER)"
            self.cursor.execute(sql_command)
            self.cursor.execute("SELECT name,length FROM 'chrNames'")
            table_chr = dict((x[0].encode('ascii'), {'length': x[1]}) for x in self.cursor.fetchall())
            for k,v in chrmeta.iteritems():
                if not(k in table_chr):
                    sql_command = "INSERT INTO 'chrNames' (name,length) VALUES ('%s',%i)"%(k,v['length'])
                    self.cursor.execute(sql_command)
                elif int(v['length']) != int(table_chr[k]['length']):
                    sql_command = "UPDATE 'chrNames' SET length=%i WHERE name='%s'"%(v['length'],k)
                    self.cursor.execute(sql_command)
        except sqlite3.OperationalError as err:
            raise Exception("Sql error: %s\n on file %s, with\n%s"%(err,self.path,sql_command))
        return True

###### chrom tables #######
    def _fix_datatables(self,fields=None):
        if fields is None:
            fields = self.fields
        try:
            fields_qry = ','.join(['"%s" %s'%(f,self.types.get(f,'text')) for f in fields if f != 'chr'])
            for chrom in self.chrmeta:
                sql_command = "CREATE TABLE IF NOT EXISTS '%s' (%s)"%(chrom,fields_qry)
                self.cursor.execute(sql_command)
                self.cursor.execute("SELECT * FROM '%s' LIMIT 1" % chrom)
                table_fields = [x[0] for x in self.cursor.description]
                for field in fields:
                    if field == 'chr' or field in table_fields: continue
                    sql_command = "ALTER TABLE '%s' ADD '%s' %s"%(chrom,field,self.types.get(field,'text'))
                    self.cursor.execute(sql_command)
                    sql_command = "CREATE INDEX IF NOT EXISTS '%s' ON '%s' (start,end)"%(chrom+"_range_idx",chrom)
                    self.cursor.execute(sql_command)
                if 'score' in fields:
                    sql_command = "CREATE INDEX IF NOT EXISTS '%s' ON '%s' (score)"%(chrom+"_score_idx",chrom)
                    self.cursor.execute(sql_command)
                if 'name' in fields:
                    sql_command = "CREATE INDEX IF NOT EXISTS '%s' ON '%s' (name)"%(chrom+"_name_idx",chrom)
                    self.cursor.execute(sql_command)
        except sqlite3.OperationalError as err:
            raise Exception("Sql error: %s\n on file %s, with\n%s"%(err,self.path,sql_command))
        return True

    def _prepare_db(self):
        status = not(self.readonly) and \
                 self._fix_attributes() and \
                 self._fix_chrmeta() and \
                 self._fix_datatables()
        self.connection.commit()
        return status

    def _get_info(self,info=None):
        if info: return info
        if 'attributes' in self.tables:
            self.cursor.execute("SELECT key,value FROM 'attributes'")
            return dict((x[0].encode('ascii'),x[1]) for x in self.cursor.fetchall())
        return {}

    def _get_fields(self,fields=None):
        if fields: return fields
        for chrom in self.chrmeta.keys():
            if chrom in self.tables:
                self.cursor.execute("SELECT * FROM '%s' LIMIT 1" % chrom)
                return [x[0].encode('ascii') for x in self.cursor.description]
        return ['chr','start','end']

    def _get_chrmeta(self,chrmeta=None):
        _chrmeta = Track._get_chrmeta(self,chrmeta)
        if _chrmeta: return _chrmeta
        try:
            sql_command = 'SELECT * FROM "chrNames"'
            query = self.cursor.execute(sql_command)
            columns = dict((x[0],n) for n,x in enumerate(query.description))
            return dict((x[columns['name']].encode('ascii'),
                         dict((k.encode('ascii'),x[n]) for k,n in columns.iteritems() if k != 'name'))
                        for x in query)
        except (sqlite3.OperationalError, sqlite3.ProgrammingError) as err:
            raise Exception("Sql error: %s\n on file %s, with\n%s" % (err,self.path,sql_command))

################################ Read ##########################################
    def _make_selection(self,selection):
        query = []
        for k,v in selection.iteritems():
            if isinstance(v,tuple):
                query.append(str(k)+" >= "+str(v[0]))
                query.append(str(k)+" <= "+str(v[1]))
            else:
                if isinstance(v,basestring):
                    str_val = "'%s'" % v
                else:
                    str_val = str(v)
                query.append(str(k)+" = "+str_val)
        return " AND ".join(query)

    def _read(self, fields, selection, order, add_chr):
        cursor = self.connection.cursor()
        chr_idx = 0
        if isinstance(selection,FeatureStream):
            chr_idx = selection.fields.index('chr')
            start_idx = selection.fields.index('start')
            end_idx = selection.fields.index('end')
        for sel in selection:
            chrom = sel[chr_idx]
            if add_chr: qfields = "'"+chrom+"' as chr,"+fields
            else: qfields = fields
            sql_command = "SELECT %s FROM '%s'" % (qfields, chrom)
            if isinstance(selection,FeatureStream):
                sql_command += " WHERE end>%s AND start<%s" % (sel[start_idx],sel[end_idx])
            elif sel[1]:
                sql_command += " WHERE %s" % self._make_selection(sel[1])
            sql_command += " ORDER BY %s" % order
            try:
                cursor.execute(sql_command)
                for x in cursor: yield x
            except sqlite3.OperationalError as err:
                raise Exception("Sql error: %s\n on file %s, with\n%s" % (err,self.path,sql_command))
        cursor.close()


    def read(self, selection=None, fields=None, order='start,end', **kw):
        """
        :param selection: list of dict of the type
            `[{'chr':'chr1','start':(12,24)},{'chr':'chr3','end':(25,45)},...]`,
            where tuples represent ranges, or a FeatureStream.
        :param fields: (list of str) list of field names (columns) to read.
        :param order: (str, comma-separated) fields with respect to which the result must
            be sorted. ['start,end']
        """
        if selection is None:
            selection = sorted(self.chrmeta.keys())
        if isinstance(selection, (basestring,dict)):
            selection = [selection]
        if isinstance(selection, (list, tuple)):
            sel2 = []
            for x in selection:
                if isinstance(x,basestring):
                    sel2.append( (str(x),{}) )
                elif isinstance(x,dict):
                    if 'chr' in x:
                        sel2.append( (x['chr'],dict((k,v) for k,v in x.iteritems()
                                                    if not(k=='chr')) ) )
                    else:
                        sel2.extend([(chrom,x) for chrom in sorted(self.chrmeta.keys())])
            selection = sel2
        elif not(selection is None or isinstance(selection,FeatureStream)):
            raise TypeError("Unexpected selection type: %s." % type(selection))
        if fields:
            add_chr = 'chr' in fields
            _fields = [f for f in fields if f in self.fields and f != 'chr']
            query_fields = ','.join(_fields)
            if add_chr: _fields = ['chr']+_fields
        else:
            query_fields = "*"
            _fields = ['chr']+[f for f in self.fields if f != 'chr']
            add_chr = True
        return FeatureStream(self._read(query_fields,selection,order,add_chr), _fields)

################################ Write ##########################################
    def _clip(self):
        for chrom,val in self.chrmeta.iteritems():
            sql_command = "UPDATE '%s' SET start=0 WHERE start<0" %(chrom)
            self.cursor.execute(sql_command)
            sql_command = "UPDATE '%s' SET end=%i WHERE end>%i" %(chrom,val['length'],val['length'])
            self.cursor.execute(sql_command)
            sql_command = "DELETE FROM '%s' WHERE end<=start" %(chrom)
            self.cursor.execute(sql_command)
        self.connection.commit()

    def write(self, source, fields=None, chrom=None, **kw):
        if not(self._prepare_db()):
            raise IOError("Cannot write database %s, readonly is %s."%(self.path,self.readonly))
        if hasattr(source, 'fields'):
            srcfields = source.fields
        elif hasattr(source, 'description'):
            srcfields = [x[0] for x in source.description]
        elif not(fields is None):
            srcfields = fields
        else:
            srcfields = self.fields
        if fields is None:
            fields = self.fields
        else:
            fields = [f for f in fields if f in self.fields]

        try:
            chr_idx = 0
            if 'chr' in srcfields:
                chr_idx = srcfields.index('chr')
            fields_left = [srcfields.index(f) for f in fields if f != 'chr' and f in srcfields]
            fields_list = ','.join([srcfields[n] for n in fields_left])
            pholders = ','.join(['?' for f in fields_left])
            def _sub(x): return [str(x[f]) for f in fields_left]
            def _sqlc(x,y,z): return "INSERT INTO '%s' (%s) VALUES (%s)" %(x,y,z)
            if chrom is None:
                if not 'chr' in srcfields:
                    raise Exception("Need a chromosome name in the source fields or in the arguments.")
                for row in source:
                    sql_command = _sqlc(row[chr_idx],fields_list,pholders)
                    self.cursor.execute(sql_command, _sub(row))
            else:
                sql_command = _sqlc(chrom,fields_list,pholders)
                if 'chr' in srcfields:
                    self.cursor.executemany(sql_command,
                                            (_sub(row) for row in source
                                             if str(row[chr_idx])==chrom))
                else:
                    self.cursor.executemany(sql_command,
                                            (_sub(row) for row in source))
            self.connection.commit()
            if kw.get('clip'): self._clip()
        except (sqlite3.OperationalError, sqlite3.ProgrammingError) as err:
            raise Exception("Sql error: %s\n on file %s, with\n%s"%(err,self.path,sql_command))

