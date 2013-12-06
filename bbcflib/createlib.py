"""
=========================
Module: bbcflib.createlib
=========================

Functions to create and manage a restrition fragments library in a 4c-seq analysis.
"""

from bein import program
from bein.util import touch
from bbcflib import genrep
from bbcflib.track import track
from bbcflib.common import cat, set_file_descr, unique_filename_in, coverageBed, gzipfile
from bbcflib.gfminer.common import sorted_stream
import os, json, re, urllib2, time, shutil

#GlobalLibPath="/archive/epfl/bbcf/data/genomes/4cLibraries"
GlobalHtsUrl="http://htsstation.epfl.ch"
GlobalRepbasePath="/archive/epfl/bbcf/data/genomes/repeats"

@program
def getRestEnzymeOccAndSeq(fasta_file, prim_site, sec_site, l_seg, l_type='typeI'):
    """
    Creates segments and fragments files of the new library from the genome sequence
    (via a call to getRestEnzymeOccAndSeq.pl).
    """
    segFile = unique_filename_in()
    fragFile = unique_filename_in()
    logFile = unique_filename_in()
#	script_path='/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_SAM/'
    progname = (l_type=='typeI') and "getRestEnzymeOccAndSeq.pl" or "getRestEnzymeOccAndSeq_typeII.pl"
    options=["-i",fasta_file,"-m",prim_site,"-s",sec_site,
             "-l",l_seg,"-o",segFile,"-f",fragFile,"-x",logFile]
    return {'arguments': [progname]+options, 'return_value': [ segFile, fragFile, logFile ]}

def parse_fragFile(fragfile,chrom_dict={}):
    """
    Parse fragment file to create segment info bed file and fragment bed file
    """
    segInfoBedFile=unique_filename_in()
    fragmentBedFile=unique_filename_in()
    o=open(segInfoBedFile,'w')
    obed=open(fragmentBedFile,'w')
    with open(fragfile,'r') as f:
        s = f.next()
        for s in f:
            if re.search('FragIsNotValid',s): continue
            s=s.strip().split('\t')
            chrom=chrom_dict.get(s[1],s[1])
            fragmentInfo='|'.join(['',s[0],chrom+':'+str(int(s[2])+1)+'-'+s[3],
                                   'indexOfSecondRestSiteOcc='+s[10],
                                   'status='+s[-1],'length='+str(int(s[3])-int(s[2])),
                                   '0','0','0','0'])
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
    return([segInfoBedFile,fragmentBedFile])


def coverageInRepeats(ex, infile, genomeName='mm9', repeatsPath=GlobalRepbasePath,
                      outdir=None, via='lsf'):
    """
    Completes the segment info bed file with the coverage in repeats of each segment.
    For now, works only for mm9, hg19 and dm3.
    """
    if not(isinstance(infile,dict)):
        infile = {"":infile}
    if outdir is None:
        resfile = unique_filename_in()+".bed"
        outf = open(resfile,'w')
    repeatsFile = os.path.join(repeatsPath, genomeName, genomeName+'_rmsk.bed')
    if not(os.path.exists(repeatsFile)):
        print("coverage in repeats not calculated as file "+repeatsFile+" does not exist.")
        if outdir is None:
            outf.close()
            cat([inf[0] for inf in infile.values()],out=resfile)
        else:
            for chrom,inf in infile.iteritems():
                shutil.copy(inf[0], os.path.join(outdir,chrom+".bed"))
            resfile = outdir
        return resfile
    futures = {}
    for chrom,inf in infile.iteritems():
        tmpfile = unique_filename_in()
        futures[chrom] = (tmpfile,coverageBed.nonblocking(ex,repeatsFile,inf[0],via=via,stdout=tmpfile))
    for chrom,fut in futures.iteritems():
        if not(outdir is None):
            resfile = os.path.join(outdir,chrom+".bed")
            outf = open(resfile,'w')
        fut[1].wait()
        coverout = track(fut[0],format='text',fields=['chr','start','end','name','c1','c2','c3','c4'])
        for s in sorted_stream(coverout.read(),[chrom]):
            s_split = s[3].split('|')
            infos = '|'.join(s_split[0:(len(s_split)-4)]+list(s[4:8]))
            outf.write('\t'.join([str(x) for x in s[0:3]+(infos,)])+'\n')
        if not(outdir is None):
            outf.close()
    if outdir is None: outf.close()
    else: resfile = outdir
    return resfile

def getEnzymeSeqId(enzyme_id,byId=False,enzymes_list=[],url=GlobalHtsUrl):
    """
    Returns the restriction site corresponding to a given enzyme id (from existing enzymes)
    """

    if len(enzymes_list) == 0 and not(url is None):
        enzymes_list.extend( json.load(urllib2.urlopen( url+"/enzymes.json" )) )
    if byId:
        src = 'site'
        trg = 'id'
        default = 0
    else:
        src = 'id'
        trg = 'site'
        default = None
    for x in enzymes_list:
        for v in x.values():
            if v[src]==enzyme_id:
                return v[trg]
    return default

def lib_exists( params, libs_list=None, url=GlobalHtsUrl ):
    """
    Return id or filename corresponding to the library described in params.
    """
    if not(isinstance(libs_list,list) or url is None):
#        with open(os.path.join(path,'libraries.json')) as f:
        libs_list = json.load(urllib2.urlopen( url+"/libraries.json" ))
    enzymes_list = []
    for lib in libs_list:
        enz1=getEnzymeSeqId(lib['library']['enzyme1_id'],False,enzymes_list,url)
        enz2=getEnzymeSeqId(lib['library']['enzyme2_id'],False,enzymes_list,url)
        if params['species'] == lib['library']['assembly_name'] \
                and params['primary'] == enz1 \
                and params['secondary'] == enz2 \
                and int(params['length']) == int(lib['library']['segment_length']) \
                and params['type'] == lib['library'].get('type','typeI'):
            return (lib['library'].get('id',0), lib['library'].get('filename'))
    return (0, None)

def createLibrary(ex, assembly_or_fasta, params, url=GlobalHtsUrl, via='local'):
    """
    Main call to create the library
    """
    if len(params['primary'])<2:
        print('Some parameters are missing, cannot create the library')
        print('primary='+params['primary']+" ; "+'secondary='+params['secondary'])
        return [None,None,None,None]

    if isinstance(assembly_or_fasta,genrep.Assembly):
        chrnames = assembly_or_fasta.chrnames
        allfiles = assembly_or_fasta.fasta_by_chrom  #assembly_or_fasta.untar_genome_fasta()
    else:
        allfiles["lib"] = assembly_or_fasta
        chrnames = ["lib"]

    libfiles = dict((c, getRestEnzymeOccAndSeq.nonblocking( ex, f,
                                                            params['primary'], params['secondary'],
                                                            params['length'],  params['type'],
                                                            via=via ))
                    for c, f in allfiles.iteritems())
    resfile = unique_filename_in()
    os.mkdir(resfile)
    bedfiles = {}
    for chrom, future in libfiles.iteritems():
        libfiles[chrom] = future.wait()
        if not os.path.getsize(libfiles[chrom][1])>0:
            time.sleep(60)
            touch(ex,libfiles[chrom][1])
        bedfiles[chrom] = parse_fragFile(libfiles[chrom][1])
    rescov = coverageInRepeats(ex, bedfiles, params['species'], outdir=resfile, via=via)
    bedchrom = [os.path.join(resfile,chrom+".bed") for chrom in chrnames]
    cat(bedchrom,out=resfile+".bed")
    gzipfile(ex,[resfile+".bed"]+bedchrom)
#    resfile_sql = resfile+".sql"
#    track.convert((resfile,'bed'),(resfile_sql,'sql'),assembly=params['species'])
    enz_list = []
    infos_lib = { 'assembly_name':  params['species'],
                  'enzyme1_id':     getEnzymeSeqId(params['primary'], True, enz_list, url),
                  'enzyme2_id':     getEnzymeSeqId(params['secondary'], True, enz_list, url),
                  'segment_length': params['length'],
                  'type':           params['type'],
                  'filename':       resfile }
    return [ libfiles, bedfiles, resfile, infos_lib ]

def get_libForGrp(ex, group, fasta_or_assembly, new_libraries, grpId, url=None, lib_dir=None, via='lsf'):
#wd_archive="/archive/epfl/bbcf/mleleu/pipeline_vMarion/pipeline_3Cseq/vWebServer_Bein/"
    def _libfile(id_lib):
        libs_list = json.load(urllib2.urlopen( url+"/libraries.json" ))
        for lib in libs_list:
            if lib['library']['id']==int(id_lib):
                return lib['library']['filename']
        return None

    def _paramsFile(paramsfile):
        """Returns a dictionary with the parameters required for the creation of a new library"""
        paramslib={'name': 'myLibrary', 'length': '30', 'type': 'typeI'}
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

    if url is None: url = GlobalHtsUrl
    if lib_dir is None: lib_dir = os.path.split(ex.remote_working_directory)[0]
    if not(group.get('library_param_file','null') in ["null",'', None]):
        library_filename = os.path.join(lib_dir,'group_'+group['name']+"_paramsFileLibrary.txt")
        paramslib = _paramsFile(library_filename)
        lib_id, ex_libfile = lib_exists( paramslib, new_libraries, url )
        if lib_id == 0 and ex_libfile == None:
            libfiles = createLibrary(ex, fasta_or_assembly, paramslib, url, via=via)
            reffile = libfiles[2]
            ex.add( libfiles[2]+".bed.gz",
                    description=set_file_descr( group['name']+"_new_library.bed.gz", groupId=grpId,
                                                step="library", type="bed" ))
#            ex.add(reffile,description=set_file_descr("new_library.sql",groupId=grpId,step="library",type="sql",view='admin'))
            new_libraries.append( {'library': libfiles[3]} )
        elif lib_id > 0:
            reffile = _libfile(lib_id)
        else:
            reffile = ex_libfile
    elif 'library_id' in group and group['library_id']> 0 and not str(group['library_id'])=="":
        reffile = _libfile(group['library_id'])
        if reffile is None:
            raise TypeError("No valid parameter passed for the library.")
        if not(os.path.exists(reffile) or os.path.exists(reffile+'.bed.gz')):
            raise TypeError("library file ("+reffile+") is not valid")
        if not os.path.exists(reffile):
            reffile += '.bed.gz'
    elif 'library_file_url' in group and group['library_file_url'] != "":
        reffile=group['library_file_url']
    else:
        reffile=None
    return reffile
