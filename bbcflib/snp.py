# Built-in modules #
import os, re, json, shutil, gzip, tarfile, pickle, urllib

# Internal modules #
from . import frontend, genrep, daflims
from bbcflib.common import get_files, cat, set_file_descr, merge_sql, gzipfile, unique_filename_in
from .track import Track, new

# Other modules #
from bein import program, ProgramFailed, MiniLIMS
from bein.util import add_pickle, touch, split_file, count_lines


@program
def gunzip_untar(file,workingDirectory):
    print file
    return {"arguments": ["tar","xvfz",file,"-C",workingDirectory],
            "return_value": None}


def untar_cat(ex,path):
    archive = tarfile.open(path)
    archive.extractall()
    allfiles=[]
    for f in archive.getmembers():
        if not f.isdir():
            allfiles.append(f.name)
    genomeRef=cat(allfiles)
    archive.close()
    return genomeRef

@program
def sam_pileup(job,bamfile,refGenome,via='lsf'):

    if(job.assembly_id == 'MLeprae_TN' or job.assembly_id == 'MSmeg_MC2_155' or job.assembly_id == 'MTb_H37Rv' or job.assembly_id == 'NA1000' or job.assembly_id == 'TB40-BAC4' ):
        ploidy=1
    else:
        ploidy=2
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
             "return_value": None}

def parse_pileupFile(ex,job,listPileupFile,allSNPpos,via='lsf',minCoverage=80,minSNP=10):
    formatedPileupFilename=unique_filename_in()
    pos=open(allSNPpos,'rb')
    posSNP=pos.readlines()
    pos.close()

    for p in listPileupFile:
        sample=open(p,'rb')
        cpt=0
        flag=0
        
        for position in posSNP:
            if position.strip("\n") in sample:
                    info=l.split("\t")
                    if int(info[7])<minSNP:
                         string="* "
                    else:
                        string=""
                    if re.search(r'[ACGT]',info[3]):
                        string+=info[3]+"\n"
                        print string
                    else:
                        print "c'est plus complique"        
        sample.close()
      
    return formatedPileupFilename


def posAllUniqSNP(ex,job,listPileupFile,minCoverage=80):
    allSNPpos=unique_filename_in()
    file=open(allSNPpos,'wb')
    #file.write("position\treference\n")
    d={}
    for p in listPileupFile:
        f=open(p,'rb')
        for l in f:
            data=l.split("\t")
            code=data[8].split()
            cpt=0
            for c in code:
                if c == "." or c == ",":
                    cpt+=1
            if (cpt*100/int(data[7]))<int(minCoverage):
                #file.write("%d\t%s\n" % (int(data[1]),data[2])) 
                d[data[1]]=data[2]

        f.close() 
    
    for pos in sorted(d.keys()):
        file.write("%s\t%s\n" % (pos,d[pos]))

    file.close()
    return allSNPpos
    
