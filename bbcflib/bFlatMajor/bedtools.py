from bein import *
from bbcflib.common import unique_filename_in

@program
def bedtools(tool, args=None):
    if args is None: args = []
    if isinstance(args,basestring):
        args = args.split(" ")
    if isinstance(args,dict):
        args2 = []
        for k,v in args.iteritems():
            k = str(k)
            if not(k.startswith("-")): k = "-"+k
            if not(isinstance(v,list)): v = [str(v)]
            args2.extend([k]+v)
        args = args2
    return {"arguments": ["bedtools",tool]+args, "return_value": None}

def _outfile(kw):
    return kw.pop('outfile',unique_filename_in())

def _wait(wait,outfile,future):
    if wait:
        future.wait()
        return outfile
    else:
        return (future,outfile)

def annotateBed(ex,intervals,files,wait=True,via='local',**kw):
    outfile = _outfile(kw):
    if not(isinstance(files,(list,tuple))): files = [files]
    kw.update({"-i":intervals,"-files": list(files)})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"annotate",kw,via=via,stdout=outfile))

def bamToBed(ex,bamfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bamfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"bamtobed",kw,via=via,stdout=outfile))

def bamToFastq(ex,bamfile,wait=True,via='local',**kw):
    outfile = kw.pop('fq',kw.pop('-fq',_outfile(kw)))
    kw.update({"-i": bamfile, "-fq": outfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"bamtofastq",kw,via=via))

def bed12ToBed6(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"bed12tobed6",kw,via=via,stdout=outfile))

def bedpeToBam(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"bedpetobam",kw,via=via,stdout=outfile))

def bedToBam(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"bedtobam",kw,via=via,stdout=outfile))
 
def bedToIgv(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"igv",kw,via=via,stdout=outfile))

def closestBed(ex,afile,bfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-a": afile,"-b": bfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"closest",kw,via=via,stdout=outfile))
 
def clusterBed(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"cluster",kw,via=via,stdout=outfile))

def complementBed(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"complement",kw,via=via,stdout=outfile))
 
def coverageBed(ex,bfile,afile=None,bamfile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(afile is None):
        kw["-a"] = afile
    elif not(bamfile is None):
        kw["-abam"] = bamfile
    else: 
        raise ValueError("Need either a bed or a bam in coverageBed.")
    kw["-b"] = bfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"coverage",kw,via=via,stdout=outfile))

def expandCols(ex,infile,column,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":infile,"-c":str(column)})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"expand",kw,via=via,stdout=outfile))
 
def fastaFromBed(ex,bedfile,fasta,wait=True,via='local',**kw):
    outfile = kw.pop('fo',kw.pop('-fo',_outfile(kw)))
    kw.update({"-fi": fasta,"-bed": bedfile,"-fo":outfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"getfasta",kw,via=via))

def flankBed(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"flank",kw,via=via,stdout=outfile))

def genomeCoverageBed(ex,genomefile,bedfile=None,bamfile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(bedfile is None):
        kw["-i"] = bedfile
    elif not(bamfile is None):
        kw["-ibam"] = bamfile
    else:
        raise ValueError("Need either a bed or a bam in genomeCoverageBed.")
    kw["-g"] = genomefile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"genomecov",kw,via=via,stdout=outfile))

def getOverlap(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"overlap",kw,via=via,stdout=outfile))
 
def groupBy(ex,bedfile,groupcol,opcol,operation,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if isinstance(groupcol,(list,tuple)):
        groupcol = ",".join([str(x) for x in groupcol])
    kw.update({"-i":bedfile,"-g":str(groupcol):"-c",str(opcol):"-o",operation})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"groupby",kw,via=via,stdout=outfile))
 
def intersectBed(ex,bfile,afile=None,bamfile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(afile is None):
        kw["-a"] = afile
    elif not(bamfile is None):
        kw["-abam"] = bamfile
    else: 
        raise ValueError("Need either a bed or a bam in intersectBed.")
    kw["-b"] = bfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"intersect",kw,via=via,stdout=outfile))

def linksBed(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"links",kw,via=via,stdout=outfile))
 
def mapBed(ex,afile,bfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-a":afile,"-b":bfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"map",kw,via=via,stdout=outfile))
 
def maskFastaFromBed(ex,bedfile,fasta,wait=True,via='local',**kw):
    outfile = kw.pop('fo',kw.pop('-fo',_outfile(kw)))
    kw.update({"-fi":fasta,"-bed":bedfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"maskfasta",kw,via=via))
 
def mergeBed(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"merge",kw,via=via,stdout=outfile))
  
def multiBamCov(ex,bedfile,bamfiles,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(isinstance(bamfiles,(list,tuple))):
        bamfiles = [bamfiles]
    kw.update({"-bed":bedfile,"-bams":list(bamfiles)})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"multicov",kw,via=via,stdout=outfile))

def multiIntersectBed(ex,bedfiles,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(isinstance(bedfiles,(list,tuple))):
        bedfiles = [bedfiles]
    kw["-i"] = list(bedfiles)
    future = 
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"multiinter",kw,via=via,stdout=outfile))
 
def nucBed(ex,bedfile,fasta,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-fi":fasta,"-bed":bedfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"nuc",kw,via=via,stdout=outfile))
 
def pairToBed(ex,bfile,afile=None,bamfile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(afile is None):
        kw["-a"] = afile
    elif not(bamfile is None):
        kw["-abam"] = bamfile
    else: 
        raise ValueError("Need either a bed or a bam in pairToBed.")
    kw["-b"] = bfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"pairtobed",kw,via=via,stdout=outfile))
 
def pairToPair(ex,afile,bfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-a":afile,"-b":bfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"pairtopair",kw,via=via,stdout=outfile))

def randomBed(ex,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-g"] = genomefile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"random",kw,via=via,stdout=outfile))

def shuffleBed(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"shuffle",kw,via=via,stdout=outfile))
 
def slopBed(ex,bedfile,genomefile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-i":bedfile,"-g":genomefile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"slop",kw,via=via,stdout=outfile))
 
def sortBed(ex,bedfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw["-i"] = bedfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"sort",kw,via=via,stdout=outfile))
 
def subtractBed(ex,afile,bfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    kw.update({"-a":afile,"-b":bfile})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"subtract",kw,via=via,stdout=outfile))
 
def tagBam(ex,bedfiles,labels,bamfile,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(isinstance(bedfiles,(list,tuple))):
        bedfiles = [bedfiles]
    if not(isinstance(labels,(list,tuple))):
        labels = [labels]
    if len(labels) < len(bedfiles):
        raise ValueError("Need one label per file in tagBam.")
    kw.update({"-i":bamfile,"-files":list(bedfiles),"-labels":list(labels)})
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"tag",kw,via=via,stdout=outfile))
 
def unionBedGraphs(ex,bgfiles,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(isinstance(bgfiles,(list,tuple))):
        bgfiles = [bgfiles]
    kw["-i"] = list(bgfiles)
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"unionbedg",kw,via=via,stdout=outfile))
 
def windowBed(ex,bfile,afile=None,bamfile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(afile is None):
        kw["-a"] = afile
    elif not(bamfile is None):
        kw["-abam"] = bamfile
    else: 
        raise ValueError("Need either a bed or a bam in windowBed.")
    kw["-b"] = bfile
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"window",kw,via=via,stdout=outfile))
 
def windowMaker(ex,bedfile=None,genomefile=None,wait=True,via='local',**kw):
    outfile = _outfile(kw)
    if not(bedfile is None):
        kw["-b"] = bedfile
    elif not(genomefile is None):
        kw["-g"] = genomefile
    else: 
        raise ValueError("Need either a bed or a genome in windowMaker.")
    return _wait(wait,outfile,
                 bedtools.nonblocking(ex,"makewindows",kw,via=via,stdout=outfile))
