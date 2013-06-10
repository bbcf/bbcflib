"""
This packages provides algorithms working on :func:`FeatureStream <bbcflib.track.FeatureStream>`.
It is divided into three major groups depending on the algorithm's return type:

* :mod:`bbcflib.gfminer.stream` returns :func:`FeatureStream <bbcflib.track.FeatureStream>` objects,
* :mod:`bbcflib.gfminer.numeric` returns numeric arrays or matrices,
* :mod:`bbcflib.gfminer.figure` returns figure files (pdf, png, etc.).

Algorithms can be used via a direct import::

    from bbcflib.gfminer import stream
    for chrom in track1.chrmeta.keys():
        catstream = stream.concatenate([track1.read(chrom),track2.read(chrom)])

or via the global :func:`run <bbcflib.gfminer.run>` function. The latter will take file names as input parameters, while a direct call requires :func:`FeatureStream <bbcflib.track.FeatureStream>` to be created beforehands.

Most functions in :mod:`gfminer <bbcflib.gfminer>` take one or more lists of :func:`FeatureStream <bbcflib.track.FeatureStream>` as parameters, plus additional algorithm-specific parameters.

"""

__all__ = ['gfminerGroup']

import os, sys
from bbcflib.track import track
from bbcflib.common import unique_filename_in

_here = 'bbcflib.gfminer.'
_module_list = ['stream','numeric','figure']
###############################################################################
class gfminerGroup(object):
    def __init__(self,members):
        self.__dict__ = members

    def __call__(self,fct,*args,**kwargs):
        if fct in self.__dict__.keys():
            return eval(fct)(*args,**kwargs)
        else:
            raise ValueError("No such function in %s: %s." %(self.__module__,fct))

    def loadable(self,fct):
        return self.__dict__.get(fct,["trackList"])

###############################################################################
def run(**kwargs):
    """
    Wrapper function to execute any operation contained in this package, directly from
    file inputs. Arguments are:

    :param operation: (str) the name of the function to be called.
    :param output: (str) a filename or a directory to write the results into.
    :param assembly: (str) a genome assembly identifier if needed.
    :param chromosome: (str) a chromosome name if operation must be restricted to a single chromsome.
    :param ...: additional parameters passed to `operation`.

    Example::

        run(operation="score_by_feature",
            output="score_output.bed", chromosome="chr1",
            trackScores="density_file.sql", trackFeatures="genes.sql")
    """
    from bbcflib import genrep
    def _map(fct):
        for module in _module_list:
            __import__(_here+module)
            smod = sys.modules[_here+module]
            if hasattr(getattr(smod, module)(),fct): return module
        return None
    funct = kwargs.pop("operation",'None')
    module = _map(funct)
    if module is None:
        raise ValueError("No such operation %s." %funct)
    output = kwargs.pop("output","./") or "./"
    if os.path.isdir(output):
        output = os.path.join(output,unique_filename_in(output)+".sql")
        format = "sql"
    else:
        format = os.path.splitext(output)[1][1:] or "sql"
    if format in ['gz','gzip']:
        format = os.path.splitext(output.strip("."+format))[1][1:]+"."+format
    smod = sys.modules[_here+module]
    trackSet = {}
    for targ in getattr(smod, module)().loadable(funct):
        trackSet[targ] = [track(t) for t in kwargs[targ].split(",")]
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
    info = None
    if 'datatype' in kwargs: info = {'datatype': kwargs.pop('datatype')}
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
            track(outf,chrmeta=chrmeta,fields=fields,
                  info=info).write(stream,chrom=chr)
        for chr in chrmeta.keys()[1:]:
            for targ in getattr(smod, module)().loadable(funct):
                kwargs[targ] = [t.read(selection=chr) for t in trackSet[targ]]
            funct_output = getattr(smod, funct)(**kwargs)
            for n,stream in enumerate(funct_output):
                track(files[n],chrmeta=chrmeta).write(stream,chrom=chr,mode='append')
    else:
        files = output
        fields = funct_output.fields
        track(files,chrmeta=chrmeta,fields=fields,
              info=info).write(funct_output,chrom=chr)
        for chr in chrmeta.keys()[1:]:
            for targ in getattr(smod, module)().loadable(funct):
                kwargs[targ] = [t.read(selection=chr) for t in trackSet[targ]]
            funct_output = getattr(smod, funct)(**kwargs)
            track(files,chrmeta=chrmeta).write(funct_output,chrom=chr,mode='append')
    return files

