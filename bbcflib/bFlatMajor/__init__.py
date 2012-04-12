__all__ = ['bFlatMajorGroup']

import os, sys
from bbcflib import genrep
from bbcflib import btrack as track
from bbcflib.common import unique_filename_in

###############################################################################
class bFlatMajorGroup(object):
    def __init__(self,members):
        self.__dict__ = members

    def __call__(self,fct,*args,**kwargs):
        if fct in self.__dict__.keys():
            return eval(fct)(*args,**kwargs)
        else:
            raise ValueError("No such function in %s: %s." %(self.__module__,fct))

    def loadable(self,fct):
        return self.__dict__.get(fct,'')

###############################################################################
_function_map = {}
for module in ['stream','numeric','figure']:
    __import__('bbcflib.bFlatMajor.'+module)
    smod = sys.modules['bbcflib.bFlatMajor.'+module]
    for k in getattr(smod, module)().__dict__.keys():
        _function_map[k] = module
###############################################################################
def run(**kwargs):
    funct = kwargs.pop("operation")
    module = _function_map[funct]
    output = kwargs.pop("output","./")
    if os.path.isdir(output):
        output = os.path.join(output,unique_filename_in(output)+".sql")
        format = "sql"
    else:
        format = os.path.splitext(output)[1][1:] or "sql"
    if format in ['gz','gzip']:
        format = os.path.splitext(output.strip("."+format))[1][1:]+"."+format
    __import__('bbcflib.bFlatMajor.'+module)
    smod = sys.modules['bbcflib.bFlatMajor.'+module]
    trackSet = {}
    for targ in getattr(smod, module)().loadable(funct):
        trackSet[targ] = [track.track(t) for t in kwargs[targ].split(",")]
    assembly = None
    if 'assembly' in kwargs:
        assembly = kwargs.pop('assembly')
    if assembly:
        chrmeta = genrep.Assembly(assembly).chrmeta
    else:
        chrmeta = trackSet[targ][0].chrmeta
        if 'chromosome' in kwargs:
            chrom = kwargs.pop('chromosome')
            chrmeta = {chrom: chrmeta.get(chrom,{})}
    chr = chrmeta.keys()[0]
    files = None
    for targ in getattr(smod, module)().loadable(funct):
        kwargs[targ] = [t.read(selection=chr) for t in trackSet[targ]]
    funct_output = getattr(smod, funct)(**kwargs)
    if isinstance(funct_output,list):
        files = []
        for n,stream in enumerate(funct_output):
            outf = "%s_%i.%s" %(output.strip(format),n,format)
            files.append(outf)
            fields = stream.fields
            track.track(outf,chrmeta=chrmeta,fields=fields).write(stream,chrom=chr)
        for chr in chrmeta.keys()[1:]:
            for targ in getattr(smod, module)().loadable(funct):
                kwargs[targ] = [t.read(selection=chr) for t in trackSet[targ]]
            funct_output = getattr(smod, funct)(**kwargs)
            for n,stream in enumerate(funct_output):
                track.track(files[n],chrmeta=chrmeta).write(stream,chrom=chr)
    else:
        files = output
        fields = funct_output.fields
        track.track(files,chrmeta=chrmeta,fields=fields).write(funct_output,chrom=chr)
        for chr in chrmeta.keys()[1:]:
            for targ in getattr(smod, module)().loadable(funct):
                kwargs[targ] = [t.read(selection=chr) for t in trackSet[targ]]
            funct_output = getattr(smod, funct)(**kwargs)
            track.track(files,chrmeta=chrmeta).write(funct_output,chrom=chr)
    return files
    
