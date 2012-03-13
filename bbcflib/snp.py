# Built-in modules #
import os, re, json, shutil, gzip, tarfile, pickle, urllib

# Internal modules #
from bbcflib import frontend, genrep, daflims
from bbcflib.common import get_files, cat, set_file_descr, merge_sql, gzipfile, unique_filename_in
from bbcflib.btrack import Track, new

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
        minSNP=10
        minCoverage=80
    else:
        ploidy=2
        minSNP=20
        minCoverage=40
    return {"arguments": ["samtools","pileup","-B","-cvsf",refGenome,"-N",str(ploidy),bamfile],
             "return_value": [minCoverage,minSNP]}

def parse_pileupFile(ex,job,dictPileupFile,allSNPpos,via='lsf',minCoverage=80,minSNP=10):
    formatedPileupFilename=unique_filename_in()
    pos=open(allSNPpos,'rb')
    posSNP=pos.readlines()
    pos.close()

    allSample={}
    iupac={'M':['A','a','C','c'],'Y':['T','t','C','c'],'R':['A','a','G','g'],'S':['G','g','C','c'],'W':['A','a','T','t'],'K':['T','t','G','g']}

    for p in dictPileupFile:
        sample=open(p,'rb')
        cpt=0
        flag=0
        s=sample.readlines()
        sample.close()
        allSample[dictPileupFile[p]]={}
        for position in posSNP:
            position=position.strip("\n")
            flag=1
            for line in s:
                if re.search(position,line):
                    flag=0
                    info=line.split("\t")
                    if int(info[7])<minSNP:
                        string="* "
                    else:
                        string=""
                    if re.search(r'[ACGT]',info[3]):
                        string+=info[3]
                        allSample[dictPileupFile[p]][position]=string
                    else:
                        snp=0
                        snp2=0

                        snp=info[8].count(iupac[info[3]][0])+info[8].count(iupac[info[3]][1])
                        snp2=info[8].count(iupac[info[3]][2])+info[8].count(iupac[info[3]][3])
                        if ((snp*100.)/int(info[7]))+((snp2*100.)/int(info[7])) > minCoverage :
                            if info[2] == iupac[info[3]][0]:
                                allSample[dictPileupFile[p]][position]=string+str(snp2*100./int(info[7]))+" % "+iupac[info[3]][2]+" / "+str((int(info[7])-snp2)*100./int(info[7]))+" % "+iupac[info[3]][0] 
                            elif info[2] == iupac[info[3]][2]:
                                allSample[dictPileupFile[p]][position]=string+str(snp*100./int(info[7]))+" % "+iupac[info[3]][0]+" / "+str((int(info[7])-snp)*100./int(info[7]))+" % "+iupac[info[3]][2]
                            else:
                                allSample[dictPileupFile[p]][position]=string+str(snp*100./int(info[7]))+" % "+iupac[info[3]][0]+" / "+str(snp2*100./int(info[7]))+" % "+iupac[info[3]][2]
                        else:
                            allSample[dictPileupFile[p]][position]="-"
         
            if flag==1:
                allSample[dictPileupFile[p]][position]="-"
    
    firstSample=allSample[allSample.keys()[0]]

    with open(formatedPileupFilename,'w') as outfile:
        outfile.write("position\treference\t"+"\t".join(dictPileupFile.values())+"\n")
        
        for p in sorted(firstSample):
            nbNoSnp=0
            for s in allSample:
                nbNoSnp+=allSample[s][p].count("-")

            if nbNoSnp!=len(allSample.keys()):
                outfile.write(str(p)+"\t"+"\t".join([allSample[s][p] for s in allSample])+"\n")

    return formatedPileupFilename

def synonymous(ex,job,allSnp):
    allCodon=unique_filename_in()
    file=open(allSnp,'rb')
    outfile=open(allCodon,'wb')
    outfile.write(file.readline())
    
    return allCodon
    

def posAllUniqSNP(ex,job,dictPileupFile,minCoverage=80):
    allSNPpos=unique_filename_in()
    file=open(allSNPpos,'wb')
    #file.write("position\treference\n")
    d={}
    for p in dictPileupFile.keys():
        f=open(p,'rb')
        for l in f:
            data=l.split("\t")
            cpt=data[8].count(".")+data[8].count(",")
            if (cpt*100./int(data[7]))<int(minCoverage):
                d[int(data[1])]=data[2]

        f.close() 
    
    for pos in sorted(d.keys()):
        file.write("%s\t%s\n" % (pos,d[pos]))

    file.close()
    return allSNPpos
    
