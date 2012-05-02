"""
======================
Module: bbcflib.createlib.py
======================

Functions used for the creation of a library in a 4c-seq analysis.
"""

from bein import program
import bbcflib.btrack as track
from bbcflib import genrep
from bbcflib.common import cat, set_file_descr, unique_filename_in, coverageBed
import os, json, re
import tarfile

GlobalLibPath="/archive/epfl/bbcf/data/genomes/4cLibraries"
GlobalRepbasePath="/archive/epfl/bbcf/data/genomes/repeats"

@program
def getRestEnzymeOccAndSeq(assembly_or_fasta, prim_site, sec_site, l_seg, l_type='typeI'):
    '''
    Creates segments and fragments files of the new library from the genome sequence (via a call to getRestEnzymeOccAndSeq.pl)
    The genome sequence (assembly_or_fasta) can be given either as a genrep assembly or as a fasta file.
    '''
    if isinstance(assembly_or_fasta,genrep.Assembly):
        fasta_path=assembly_or_fasta.fasta_path()
        tar = tarfile.open(fasta_path)
        tar.extractall()
        allfiles=[]
        for finfo in tar.getmembers():
            if not finfo.isdir():
                allfiles.append(finfo.name)
        fasta_file=cat(allfiles)
        tar.close()
    else:
        fasta_file=assembly_or_fasta

    segFile = unique_filename_in()
    fragFile = unique_filename_in()
    logFile = unique_filename_in()
    outfiles=[segFile, fragFile, logFile, fasta_file]
#	script_path='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_SAM/'
    progname = (l_type=='typeI') and "getRestEnzymeOccAndSeq.pl" or "getRestEnzymeOccAndSeq_typeII.pl"
    options=[progname,
             "-i",fasta_file,"-m",prim_site,"-s",sec_site,
             "-l",l_seg,"-o",segFile,"-f",fragFile,"-x",logFile]
    return {'arguments': options, 'return_value':outfiles}

def parse_fragFile(fragfile):
    '''
    Parse fragment file to create segment info bed file and fragment bed file
    '''
    segInfoBedFile=unique_filename_in()
    fragmentBedFile=unique_filename_in()
    segmentBedFile=unique_filename_in()
    o=open(segInfoBedFile,'w')
    obed=open(fragmentBedFile,'w')
    oseg=open(segmentBedFile,'w')
    with open(fragfile,'r') as f:
        s = f.next()
        for s in f:
            if re.search('FragIsNotValid',s): continue
            s=s.strip().split('\t')
            chrom="chr"+s[1]
            fragmentInfo='|'.join(['',s[0],chrom+':'+str(int(s[2])+1)+'-'+s[3],
                                   'indexOfSecondRestSiteOcc='+s[10],
                                   'status='+s[-1],'length='+str(int(s[3])-int(s[2])),
                                   0,0,0,0])
            o.write('\t'.join([chrom,s[5],s[6],'type=startSegment'+fragmentInfo])+'\n')
            o.write('\t'.join([chrom,s[8],s[9],'type=endSegment'+fragmentInfo])+'\n')
            row = [chrom,s[2],s[3],'frag'+s[0]]
            obed.write('\t'.join(row)+'\n')
            row[1:3]=[s[5],s[6]]
            obed.write('\t'.join(row)+'_startSeq\n')
            row[1:3]=[s[8],s[9]]
            obed.write('\t'.join(row)+'_endSeq\n')
    o.close()
    obed.close()
    oseg.close()
    return([segInfoBedFile,fragmentBedFile,segmentBedFile])


def coverageInRepeats(ex,infile,genomeName='mm9',repeatsPath=GlobalRepbasePath,via='lsf'):
    '''
    Completes the segment info bed file with the coverage in repeats of each segment.
    For now, works only for mm9, hg19 and dm3.
    '''
    repeatsFile=os.path.join(repeatsPath,genomeName,genomeName+'_rmsk.bed')
    if not(os.path.exists(repeatsFile)):
        print("coverage in repeats not calculated as file "+repeatsFile+" does not exist.")
        return(infile)

    tmpfile = unique_filename_in()
    coverageBed.nonblocking(ex,repeatsFile,infile,via=via,stdout=tmpfile).wait()
    resfile = unique_filename_in()+".bed"
    with open(resfile,'w') as o:
        with open(tmpfile,'r') as f:
            for s in f:
                s=s.strip().split('\t')
                s_split=s[3].split('|')
                infos='|'.join(s_split[0:(len(s_split)-4)]+s[4:8])
                o.write('\t'.join(s[0:3]+[infos])+'\n')
    return resfile

def getEnzymeSeqId(enzyme_id,byId=False,enzymes_dict=None,libpath=GlobalLibPath):
    '''
    Returns the restriction site corresponding to a given enzyme id (from existing enzymes)
    '''
    if not(isinstance(enzymes_dict,list)):
        with open(os.path.join(libpath,'enzymes.json')) as g:
            enzymes_dict=json.load(g)
    if byId:
        src = 'site'
        trg = 'id'
        default = 0
    else:
        src = 'id'
        trg = 'site'
        default = None
    for x in enzymes_dict:
        for v in x.values():
            if v[src]==enzyme_id:
                return v[trg]
    return default

def lib_exists(params,path=GlobalLibPath,returnType="id"):
    '''
    Return id or filename corresponding to the library described in params.
    '''
    if not(isinstance(libs_dict,dict)):
        with open(os.path.join(path,'libraries.json')) as f:
            libs_dict = json.load(f)
    with open(os.path.join(path,'enzymes.json')) as g:
        enzymes_dict=json.load(g)
    for lib in libs_dict:
        enz1=getEnzymeSeqId(lib['library']['enzyme1_id'],False,enzymes_dict)
        enz2=getEnzymeSeqId(lib['library']['enzyme2_id'],False,enzymes_dict)
		#temporary...
        if not 'type' in lib['library']: lib['library']['type']='typeI'
        if params['species']==lib['library']['assembly_name'] and params['primary']==enz1 and params['secondary']==enz2 and int(params['length'])==int(lib['library']['segment_length']) and params['type']==lib['library']['type']:
            if returnType=="id":
                return lib['library']['id']
            else:
                return lib['library']['filename']

    if returnType == "id" :
        return 0
    else:
        return None

# *** main call to create the library
def createLibrary(ex,fasta_allchr,params,via='local'):
    '''
    Main call to create the library
    '''
    if len(params['primary'])<2 or len(params['secondary'])<2:
        print('Some parameters are missing, cannot create the library')
        print('primary='+params['primary']+" ; "+'secondary='+params['secondary'])
        return [None,None,None,None]
    libfiles=getRestEnzymeOccAndSeq(ex,fasta_allchr,params['primary'],params['secondary'],
                                    params['length'], params['type'])
    bedfiles=parse_fragFile(libfiles[1])
    resfile=coverageInRepeats(ex,bedfiles[0],params['species'],via=via)
    resfile_sql=resfile+".sql"
    track.convert((resfile,'bed'),(resfile_sql,'sql'),assembly=params['species'])
    infos_lib={'assembly_name':params['species'],
               'enzyme1_id':getEnzymeSeqId(params['primary'],True),
               'enzyme2_id':getEnzymeSeqId(params['secondary'],True),
               'segment_length':params['length'],
               'type':params['type'],
               'filename':resfile}
    return [libfiles,bedfiles,resfile,infos_lib,resfile_sql]

def get_libForGrp(ex,group,fasta_or_assembly,new_libraries, job_id, grpId, lib_dir=None):
#wd_archive="/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/" #temporary: will be /scratch/cluster/monthly/htsstation/4cseq/job.id
#os.path.split(ex.remote_working_directory)[0]
    def _libfile(id_lib):
        with open(os.path.join(GlobalLibPath,'libraries.json')) as f: libs_dict = json.load(f)
            #id_lib=13
        for lib in libs_dict:
            if lib['library']['id']==int(id_lib):
                return lib['library']['filename']
        return None

    def _paramsFile(paramsfile):
        '''Returns a dictionary with the parameters required for the creation of a new library'''
        paramslib={'name': 'myLibrary',
                   'length': '30',
                   'type': 'typeI'}
        with open(paramsfile) as f:
            for s in f:
                s=s.strip().split('=')
                key = None
                if   re.search('Library name',s[0],re.I) and len(s[1])>1:   key='name'
                elif re.search('Genome name',s[0],re.I):                    key='species'
                elif re.search('Primary',s[0],re.I):                        key='primary'
                elif re.search('Secondary',s[0],re.I):                      key='secondary'
                elif re.search('Segment length',s[0],re.I) and len(s[1])>0: key='length'
                elif re.search('Type',s[0],re.I) and len(s[1])>1:           key='type'
                if key: paramslib[key]=s[1]
        return paramslib

    if lib_dir is None: lib_dir = os.path.split(ex.remote_working_directory)[0]
    libfile = group.get('library_param_file',False)
    if str(group['library_param_file']).lower() in ['1','true','on','t']: libfile = True
    if libfile:
        library_filename = os.path.join(lib_dir,'group_' + group['name'] + "_paramsFileLibrary.txt")
        paramslib=_paramsFile(library_filename);
        lib_id=lib_exists(paramslib)
        ex_libfile=lib_exists(paramslib,new_libraries,returnType="filename")
        if lib_id == 0 and ex_libfile == None :
            libfiles=createLibrary(ex,fasta_or_assembly,paramslib);
            reffile=libfiles[4]
            ex.add(libfiles[2],description=set_file_descr("new_library.bed",groupId=grpId,step="library",type="bed"))
            ex.add(reffile,description=set_file_descr("new_library.sql",groupId=grpId,step="library",type="sql",view='admin'))
            new_libraries.append({'library':libfiles[3]})
        elif lib_id > 0 :
            reffile=_libfile(lib_id)+".sql"
        else:
            reffile=ex_libfile+".sql"
    elif 'library_id' in group and group['library_id']> 0 and not str(group['library_id'])=="":
        reffile=_libfile(group['library_id'])
        if reffile is None:
            raise TypeError("No valid parameter passed for the library.")
        if not(os.path.exists(reffile) or os.path.exists(reffile+'.bed.gz')):
            raise TypeError("library file ("+reffile+") is not valid")
        if not os.path.exists(reffile):
            reffile=reffile+'.bed.gz'
    elif 'library_file_url' in group and group['library_file_url'] != "" :
        reffile=group['library_file_url']
    else:
        reffile=None
    return reffile
