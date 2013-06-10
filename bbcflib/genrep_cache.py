'''
Build a local database that store queries to GenRep. Once stored, result are fetched from it
instead of from GenRep web-service.

Use it :
from bbcflib import genrep_cache as genrep
print genrep.Assembly('mm9').chrnames

Warning :
You can only call Assembly and all methods belonging to it.
Direct calls to GenRep() or GenRepObject() are not supported

Functioning :
On import it create the configuration file `.genrepcache.conf` in your home directory with
parameters :
    db_name = Database name.
    time_limit = Time limit when a cached response must be reloaded (in days).
    assembly_cache = List of methods you want to cache.
    unique_ids = Important attributes in your class that will determine the uniqueness of the child object
    connector = Default connector : change it to what you want and supported by sqlalchemy
        (see http://docs.sqlalchemy.org/en/latest/core/engines.html#supported-databases).

Then on method call, if the method is not in `assembly_cache` parameter, your call is redirected
to the original GenRep library, else the response is fetched from the local database (if query already logged
and not outdated)

'''


import genrep, hashlib, datetime, os, ConfigParser, json, pickle

from sqlalchemy import create_engine, and_, MetaData, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Text, PickleType
from sqlalchemy.types import DateTime
from sqlalchemy.orm import sessionmaker






'''
Directory where genrep_cache db will be installed (if connector is sqlite)
'''
home_dir =  os.getenv('HOME') or os.getenv('USERPROFILE', '.')

'''
Configuration file name.
'''
conf_file = '.genrepcache.conf'

'''
Configuration path.
'''
conf_path = os.path.join(home_dir, conf_file)


default_conf = {'db_name' : 'genrepcache.sqlite',
                'time_limit' : 1,
                'assembly_cache' : ['chrnames', 'fasta_path', 'get_sqlite_url', 'sqlite_path'],
                'unique_ids' : ['intype', 'nr_assembly_id', 'source_id'],
                }

default_conf['connector'] = 'sqlite:///%s' % os.path.join(home_dir, default_conf['db_name'])

'''
Section in the config file.
'''
section = 'GenRep cache'



def dec(obj, func, _property=False):
    '''
    Hook to fetch response from cache if there is one.
    '''
    def wrapper(*args, **kw):
        uid = obj.unique_id()
        if _property : args = ()
        return cached(obj, func, uid, args, kw, _property)
    return wrapper


def cached_ass(*args, **kw):
    '''
    Cache assembly
    '''
    session = get_session()
    ass = session.query(Ass).filter(and_(Ass.args == str(args),
                                     Ass.kw == str(kw))).first()
    if ass is not None:
        now = datetime.datetime.now()
        if timedelta > (now - ass.date):
            session.close()
            return pickle.loads(str(ass.assembly))
    else :
        ass = Ass()

    assembly = genrep.Assembly(*args, **kw)
    ass.args = str(args)
    ass.kw = str(kw)
    ass.assembly = pickle.dumps(assembly)
    session.add(ass)
    session.commit()
    session.close()
    return assembly

def cached(obj, func, uid, args, kw, _property):
    '''
    Look if the method is already cached and not outdated
    then return the response
    '''
    session = get_session()
    store = session.query(Store).filter(and_(Store.uid == uid,
                                     Store.args == str(args),
                                     Store.kw == str(kw))).first()
    # look if the response is here / not outdated
    if store is not None:
        now = datetime.datetime.now()
        if timedelta > (now - store.date):
            session.close()
            return store.response
    else :
        store = Store()

    # get the response from genRep and store it
    if _property :
        response = getattr(obj.gr, func)
    else :
        response = getattr(obj.gr, func)(*args, **kw)
    store.uid = uid
    store.args = str(args)
    store.kw = str(kw)
    store.response = str(response)

    session.add(store)
    session.commit()
    session.close()
    return response

class Assembly(object):
    '''
    Assembly object that reflect the `original`.
    '''
    def __init__(self, *args, **kw):
        '''
        Initialize the decoration system.
        Each method in `assembly_cache` will be decorated by the
        `dec` method. Others will be redirected to the Assembly decorated.
        '''
        self.uid = None
        self.gr = cached_ass(*args, **kw)


        for func in assembly_cache:
            if hasattr(self.gr, func):
                v = self.gr.__class__.__dict__.get(str(func))
                if v is None :
                    ## instance
                    raise AttributeError('Only function or property can be decorated : "%s" is an attribute' % func)
                ## class
                elif isinstance(v, property):
                    prop = property(dec(self, str(func), _property=True))
                    setattr(self.__class__, str(func), prop)
                else :
                    setattr(self, str(func), dec(self, func))
            else :
                raise AttributeError('%r object has no attribute %r' % (type(self.gr).__name__, func))


    def __getattr__(self, name):
        '''
        Redirect attribute search on the decorated parameter.
        '''
        if hasattr(self.gr, name):
            return getattr(self.gr, name)
        raise AttributeError('%r object has no attribute %r' % (type(self).__name__, name))




    def unique_id(self):
        '''
        Generate an unique id with the `unique_ids attributes`.
        '''
        if self.uid is None:
            attrs = ''
            for att in unique_ids:
                attrs += str(getattr(self.gr, att))
            self.uid = hashlib.sha1(attrs).hexdigest()
        return self.uid



'''
Database configuration
'''

Base = declarative_base()

class Store(Base):
    '''
    Data store.
    '''
    __tablename__ = 'grc'
    uid = Column(Text, primary_key=True)
    args = Column(Text, primary_key=True)
    kw = Column(Text, primary_key=True)
    response = Column(Text)
    date = Column(DateTime, default=datetime.datetime.now, nullable=False)

class Ass(Base):
    '''
    Data store.
    '''
    __tablename__ = 'ass'
    args = Column(Text, primary_key=True)
    kw = Column(Text, primary_key=True)
    assembly = Column(Text)
    date = Column(DateTime, default=datetime.datetime.now, nullable=False)


metadata = MetaData()
grc = Table('grc', metadata,
     Column('uid', Text, primary_key=True),
     Column('args', Text, primary_key=True),
     Column('kw', Text, primary_key=True),
     Column('response', Text),
     Column('date', DateTime))

ass = Table('ass', metadata,
    Column('args', Text, primary_key=True),
    Column('kw', Text, primary_key=True),
    Column('assembly', PickleType),
    Column('date', DateTime))


def init(connector, metadata, tables, echo=False):
    '''
    Initialization of the engine.
    '''
    engine = create_engine(connector, echo=echo)
    metadata.create_all(bind=engine, tables=tables, checkfirst=True)
    Session = sessionmaker(bind=engine)
    Session.configure(bind=engine)
    return Session


def configuration(config):
    if not config.read(conf_path):
        config.add_section(section)
        for k, v in default_conf.iteritems():
            config.set(section, k, json.dumps(v))
        with open(conf_path, 'w') as cf:
            config.write(cf)
        configuration(config)


config = ConfigParser.ConfigParser()
configuration(config)


nb = config.getint(section, 'time_limit')
connector = json.loads(config.get(section, 'connector'))
assembly_cache = json.loads(config.get(section, 'assembly_cache'))
json.loads(config.get(section, 'assembly_cache'))

unique_ids = json.loads(config.get(section, 'unique_ids'))
timedelta = datetime.timedelta(days=nb)
DBSession = init(connector, metadata, [grc, ass], echo=False)

def get_session():
    return DBSession()



if __name__ == '__main__':
    session = get_session()
    a = Assembly('mm9')
#    print 'CHRNAMES'
#    print a.chrnames
#    print 'FASTA PATH'
#    print a.fasta_path('chr1')
    print 'AGAIN'
    print a.chromosomes




