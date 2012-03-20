import genrep, hashlib

#from sqlalchemy import create_engine, MetaData, Table
#from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy import Column, Integer, String, Float, Text
#
#
#Base = declarative_base()
#
#class Store(Base):
#    __tablename__ = 'grc'
#    hash = Column(Text, primary_key=True)
#    args = Column(Text, primary_key=True)
#    response = Column(Text)
#    Session = sessionmaker(bind=engine)                                                                              
#    Session.configure(bind=engine)                                                                                   
#    session = Session()                                                                                              
#    return session

#def init(connector, echo=False):
#    return create_engine(connector, echo=echo)

'''
List of methods you want to cache
'''
assembly_cache = ['chrnames']

'''
Important attributes in your classe that will determine the uniquiness of the child object. 
'''
unique_ids = ['assembly', 'intype']

def cache_it(obj, func):
    return getattr(obj, func)


        
class Assembly(object):
    def __init__(self, *args, **kw):
        '''
        Initialize the decoration system.
        Each method in `assembly_cache` will be decorated by the 
        `chache_it` method. Others will be redirected to the Aeembly decorated.
        '''
        self.gr = genrep.Assembly(*args, **kw)
        for func in assembly_cache:
            if hasattr(self.gr, func):
                setattr(self, func, cache_it(self.gr, func))
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
        if self.uid is None:
            attrs = ''
            for att in unique_ids:
                attrs += str(getattr(self.gr, att))
            self.uid = hashlib.sha1(attrs).hexdigest()
        return self.uid
    