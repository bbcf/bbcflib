from bein import *
from bbcflib.common import unique_filename_in

@program
def bedtools(tool, args=None):
    return {"arguments": ["bedtools",tool]+args, "return_value": None}

def annotateBed(ex,intervals,files,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(isinstance(files,(list,tuple))):
        files = [files]
    args += ["-i",intervals,"-files"]+list(files)
    future = bedtools.nonblocking(ex,"annotate",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def bamToBed(ex,bamfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bamfile]
    future = bedtools.nonblocking(ex,"bamtobed",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def bamToFastq(ex,bamfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    elif 'fq' in kw:
        outfile = kw.pop('fq')
    else:
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bamfile,"-fq",outfile]
    future = bedtools.nonblocking(ex,"bamtofastq",args,via=via)
    future.wait()
    return outfile

def bed12ToBed6(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"bed12tobed6",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def bedpeToBam(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"bedpetobam",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def bedToBam(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"bedtobam",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def bedToIgv(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"igv",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def closestBed(ex,afile,bfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-a",afile,"-b",bfile]
    future = bedtools.nonblocking(ex,"closest",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def clusterBed(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"cluster",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def complementBed(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"complement",args,via=via,stdout=outfile)
    future.wait()
    return outfile 
 
def coverageBed(ex,bfile,afile=None,bamfile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(afile is None):
        args += ["-a",afile]
    elif not(bamfile is None):
        args += ["-abam",bamfile]
    else: 
        raise ValueError("Need either a bed or a bam in coverageBed.")
    args += ["-b",bfile]
    future = bedtools.nonblocking(ex,"coverage",args,via=via,stdout=outfile)
    future.wait()
    return outfile  

def expandCols(ex,infile,column,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = ["-i",infile,"-c",str(column)]
    future = bedtools.nonblocking(ex,"expand",args,via=via,stdout=outfile)
    future.wait()
    return outfile  
 
def fastaFromBed(ex,bedfile,fasta,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    elif 'fo' in kw:
        outfile = kw.pop('fo')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-fi",fasta,"-bed",bedfile,"-fo",outfile]
    future = bedtools.nonblocking(ex,"getfasta",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def flankBed(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"flank",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def genomeCoverageBed(ex,genomefile,bedfile=None,bamfile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(bedfile is None):
        args += ["-i",bedfile]
    elif not(bamfile is None):
        args += ["-ibam",bamfile]
    else:
        raise ValueError("Need either a bed or a bam in genomeCoverageBed.")
    args += ["-g",genomefile]
    future = bedtools.nonblocking(ex,"genomecov",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def getOverlap(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"overlap",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def groupBy(ex,bedfile,groupcol,opcol,operation,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if isinstance(groupcol,(list,tuple)):
        groupcol = ",".join([str(x) for x in groupcol])
    args += ["-i",bedfile,"-g",str(groupcol),"-c",str(opcol),"-o",operation]
    future = bedtools.nonblocking(ex,"groupby",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def intersectBed(ex,bfile,afile=None,bamfile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(afile is None):
        args += ["-a",afile]
    elif not(bamfile is None):
        args += ["-abam",bamfile]
    else: 
        raise ValueError("Need either a bed or a bam in intersectBed.")
    args += ["-b",bfile]
    future = bedtools.nonblocking(ex,"intersect",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def linksBed(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()+".html"
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,outfile]
    future = bedtools.nonblocking(ex,"links",args,via=via)
    future.wait()
    return outfile
 
def mapBed(ex,afile,bfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-a",afile,"-b",bfile]
    future = bedtools.nonblocking(ex,"map",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def maskFastaFromBed(ex,bedfile,fasta,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    elif 'fo' in kw:
        outfile = kw.pop('fo')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-fi",fasta,"-bed",bedfile,"-fo",outfile]
    future = bedtools.nonblocking(ex,"maskfasta",args,via=via,stdout=outfile)
    future.wait()
    return outfile 
 
def mergeBed(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"merge",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
multiBamCov(ex,bedfile,bamfiles,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(isinstance(bamfiles,(list,tuple))):
        bamfiles = [bamfiles]
    args += ["-bed",bedfile,"-bams"]+list(bamfiles)
    future = bedtools.nonblocking(ex,"multicov",args,via=via,stdout=outfile)
    future.wait()
    return outfile

def multiIntersectBed(ex,bedfiles,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(isinstance(bedfiles,(list,tuple))):
        bedfiles = [bedfiles]
    args += ["-i"]+list(bedfiles)
    future = bedtools.nonblocking(ex,"multiinter",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def nucBed(ex,bedfile,fasta,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    elif 'fo' in kw:
        outfile = kw.pop('fo')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-fi",fasta,"-bed",bedfile]
    future = bedtools.nonblocking(ex,"nuc",args,via=via,stdout=outfile)
    future.wait()
    return outfile 
 
def pairToBed(ex,bfile,afile=None,bamfile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(afile is None):
        args += ["-a",afile]
    elif not(bamfile is None):
        args += ["-abam",bamfile]
    else: 
        raise ValueError("Need either a bed or a bam in pairToBed.")
    args += ["-b",bfile]
    future = bedtools.nonblocking(ex,"pairtobed",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def pairToPair(ex,afile,bfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-a",afile,"-b",bfile]
    future = bedtools.nonblocking(ex,"pairtopair",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def randomBed(ex,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-g",genomefile]
    future = bedtools.nonblocking(ex,"random",args,via=via,stdout=outfile)
    future.wait()
    return outfile 

def shuffleBed(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"shuffle",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def slopBed(ex,bedfile,genomefile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile,"-g",genomefile]
    future = bedtools.nonblocking(ex,"slop",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def sortBed(ex,bedfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-i",bedfile]
    future = bedtools.nonblocking(ex,"sort",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def subtractBed(ex,afile,bfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    args += ["-a",afile,"-b",bfile]
    future = bedtools.nonblocking(ex,"subtract",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def tagBam(ex,bedfiles,labels,bamfile,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(isinstance(bedfiles,(list,tuple))):
        bedfiles = [bedfiles]
    if not(isinstance(labels,(list,tuple))):
        labels = [labels]
    if len(labels) < len(bedfiles):
        raise ValueError("Need one label per file in tagBam.")
    args += ["-i",bamfile,"-files"]+list(bedfiles)+["-labels"]+list(labels)
    future = bedtools.nonblocking(ex,"tag",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def unionBedGraphs(ex,bgfiles,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(isinstance(bgfiles,(list,tuple))):
        bgfiles = [bgfiles]
    args += ["-i"]+list(bgfiles)
    future = bedtools.nonblocking(ex,"unionbedg",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def windowBed(ex,bfile,afile=None,bamfile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(afile is None):
        args += ["-a",afile]
    elif not(bamfile is None):
        args += ["-abam",bamfile]
    else: 
        raise ValueError("Need either a bed or a bam in windowBed.")
    args += ["-b",bfile]
    future = bedtools.nonblocking(ex,"window",args,via=via,stdout=outfile)
    future.wait()
    return outfile
 
def windowMaker(ex,bedfile=None,genomefile=None,via='local',**kw):
    if 'outfile' in kw:
        outfile = kw.pop('outfile')
    else: 
        outfile = unique_filename_in()
    args = []
    for k,v in kw.iteritems():
        args.extend(["-"+str(k),str(v)])
    if not(bedfile is None):
        args += ["-b",bedfile]
    elif not(genomefile is None):
        args += ["-g",genomefile]
    else: 
        raise ValueError("Need either a bed or a genome in windowMaker.")
    future = bedtools.nonblocking(ex,"makewindows",args,via=via,stdout=outfile)
    future.wait()
    return outfile 
