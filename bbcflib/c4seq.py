"""
=======================
Module: bbcflib.c4seq
=======================

This module provides functions to run a 4c-seq analysis
from reads mapped on a reference genome.

"""

import sys, os, json, re, tarfile, time
from bein import program
from bein.util import touch
from bbcflib.createlib import get_libForGrp
from bbcflib.track import track, convert
from bbcflib.mapseq import parallel_density_sql
from bbcflib.common import unique_filename_in, gzipfile, merge_sql, cat, gfminer_run, set_file_descr

from shutil import copyfile

# *** Create a dictionary with infos for each primer (from file primers.fa)
# ex: primers_dict=loadPrimers('/archive/epfl/bbcf/data/DubouleDaan/finalAnalysis/XmNGdlXjqoj6BN8Rj2Tl/primers.fa')
def loadPrimers(primersFile):
    '''Create a dictionary with infos for each primer (from file primers.fa) '''
    primers = {}
    name = ''
    with open(primersFile,'rb') as f:
        for s in f:
            s = re.sub(r'\s','',s)
            if not(s): continue
            if s[0] != ">":
                if name:
                    primers[name]['seq'] = s
                    name = ''
                continue
            infos = s.split('|')
            name = infos[0][1:]
            primerInfos = { 'fullseq': infos[1],
                            'baitcoord': infos[2],
                            'primary': infos[3],
                            'seq': '' }
            if 'Exclude' in infos[-1]:
                primerInfos['regToExclude'] = infos[-1].split('=')[1]
            primerInfos['seqToFilter'] = infos[4:-1]
            if not name in primers: primers[name] = primerInfos
    return primers

def removeNA( fileToClean ):
    ''' remove NA present in the 4th column of a file '''
    ''' mainly used with bedgraph'''
    fileNoNA = unique_filename_in()
    resfile = open(fileNoNA, 'w')
    with open( fileToClean ) as f:
        for s in f:
            if s[0:5] == 'track': resfile.write(s)
            if s[0:5] != 'track' and s.strip().split('\t')[3] != "NA": resfile.write(s)
    resfile.close()
    return fileNoNA

@program
def segToFrag( countsPerFragFile, regToExclude="", script_path='' ):
    '''
    This function calls segToFrag.awk (which transforms the counts per segment to a normalised count per fragment).
    Provide a region to exclude if needed.
    '''
    args = ["awk","-f",os.path.join(script_path,'segToFrag.awk')]
    if regToExclude: args += ["-v","reg2Excl="+regToExclude]
    return {'arguments': args+[countsPerFragFile], 'return_value': None}

@program
def BRICKSToFrag( BRICKSfile, fragFile, resfile, script_path='' ):
    '''
    This function calls BRICKS2Frags.sh (which return for each fragment the -log10(pval) of the associated BRICKS)
    For now it only concerns fragments located on the viewpoint chromosome as the domainogram is not run on 'trans' fragments
    '''
    args = ['BRICKS2Frags.sh', BRICKSfile, fragFile, resfile]
    print(args)
    return {'arguments': args, 'return_value': None}


def _RCMD(path,script):
    return ["R","--vanilla","--slave","-f",os.path.join(path,script),"--args"]

@program
#--args segToFrag_PFBM8_Sp9_all.bedGraph normFile.bedGraph chr2:73114371-73115873 1000000 name chr2:73112371-73118419
def normFrags( inputFile, outputFile, baitCoord, sizeExt, name, regToExclude, script_path=''):
    args = _RCMD(script_path,"normFrags.R")+[inputFile,outputFile,baitCoord,str(sizeExt),name, regToExclude]
    print(args)
    return {'arguments': args, 'return_value': None}

@program
#mergeRep.R --args file1;file2;file3 mergedfile.bedGraph 4 chr2:73911294-73916104
def mergeRep(inputFiles, outputFile, regToExclude, idCol="4", name="", script_path=''):
    #time.sleep(20)
    args = _RCMD(script_path,"mergeRep.R")+[inputFiles,outputFile,idCol,name,regToExclude]
    print(args)
    return {'arguments': args, 'return_value': None}

@program
#combineTracks.R --args $fragsFiles combinedFile.txt $names 4,5 $defVal $all_regToExclude $reffile
def makeTable(inputFiles, outputFile, names, idCols="4", out_chromosomes = 'NA', defVal=0, all_regToExclude='NA', script_path=''):
    #time.sleep(30)
    args = _RCMD(script_path,"combineTracks.R")+[inputFiles,outputFile,names,idCols,str(defVal),out_chromosomes, all_regToExclude]
    print(args)
    return {'arguments': args, 'return_value': None}

@program
def profileCorrection( inputFile, baitCoord, name, outputFile, reportFile, tableFile, script_path='' ):
    #time.sleep(20)
    args = _RCMD(script_path,"profileCorrection.R")+[inputFile,baitCoord,name,outputFile,reportFile,tableFile]
    print(args)
    return {'arguments': args, 'return_value': None}

@program
def smoothFragFile( inputFile, nFragsPerWin, curName, outputFile, regToExclude="", script_path='' ):
    #time.sleep(30)
    args = _RCMD(script_path,"smoothData.R")+[inputFile,str(nFragsPerWin),curName,outputFile,regToExclude]
    print(args)
    return {'arguments': args, 'return_value': None}

@program
def runDomainogram( infile, name, prefix=None, regCoord="",
                    wmaxDomainogram=500, wmax_BRICKS=50, skip=0, script_path='' ):
    #time.sleep(60)
    if prefix is None: prefix = name
    args = _RCMD(script_path,"runDomainogram.R")+[infile,name,prefix,regCoord,str(wmaxDomainogram),str(wmax_BRICKS),str(skip),script_path]
    print(args)
    return {'arguments': args, 'return_value': prefix+".log"}


def density_to_countsPerFrag( ex, file_dict, groups, assembly, regToExclude, script_path, via='lsf' ):
    '''
    Main function to compute normalised counts per fragments from a density file.
    '''
    futures = {}
    results = {}
    for gid, group in groups.iteritems():
        reffile = file_dict['lib'][gid]
        futures[gid] = {}
        results[gid] = {}
        for rid,run in group['runs'].iteritems():
            density_file = file_dict['4cseq']['density_files'][gid][rid]
            gm_futures = []
            for ch in assembly.chrnames:
                chref = os.path.join(reffile,ch+".bed.gz")
                if not(os.path.exists(chref)): chref = reffile
    #            features = track(chref,'bed')
    #            outbed.write(gMiner.stream.mean_score_by_feature(
    #                    scores.read(selection=ch),
    #                    features.read(selection=ch)), mode='append')
                bedfile = unique_filename_in()+".bed"
                gfminer_job = {"operation": "score_by_feature",
                               "output": bedfile,
                               "datatype": "qualitative",
                               "args": "'"+json.dumps({"trackScores":density_file,
                                                       "trackFeatures":chref,
                                                       "chromosome":ch})+"'"}
                gm_futures.append((gfminer_run.nonblocking(ex,gfminer_job,via=via),
                                   bedfile))
            outsql = unique_filename_in()+".sql"
            sqlouttr = track( outsql, chrmeta=assembly.chrmeta,
                              info={'datatype':'quantitative'},
                              fields=['start', 'end', 'score'] )
            outbed_all = []
            for n,f in enumerate(gm_futures):
                f[0].wait()
                fout = f[1]
                if not(os.path.exists(fout)):
                    time.sleep(60)
                    touch(ex,fout)
                outbed_all.append(fout)
                outbed = track(fout, chrmeta=assembly.chrmeta)
                sqlouttr.write( outbed.read(fields=['start', 'end', 'score'],
                                            selection={'score':(0.01,sys.maxint)}),
                                chrom=assembly.chrnames[n] )
            sqlouttr.close()
            countsPerFragFile = unique_filename_in()+".bed"
            countsPerFragFile = cat(outbed_all,out=countsPerFragFile)
            results[gid][rid] = [ countsPerFragFile, outsql ]
            FragFile = unique_filename_in()
            touch(ex,FragFile)
            futures[gid][rid] = (FragFile,
                            segToFrag.nonblocking( ex, countsPerFragFile, regToExclude[gid],
                                                   script_path, via=via, stdout=FragFile ,
                                                   memory=4 ))
    def _parse_select_frag(stream):
        for s in stream:
            sr = s.strip().split('\t')
            if 'IsValid' in sr[2] and not any([w in sr[8] for w in ['_and_','BothRepeats','notValid']]):
                patt = re.search(r'([^:]+):(\d+)-(\d+)',sr[1])
                if patt:
                    coord = patt.groups()
#                    if float(sr[11])>0.0:
                    yield (coord[0], int(coord[1])-1, int(coord[2]), float(sr[11]))

    for gid, dict_gid in futures.iteritems():
        for rid, res in dict_gid.iteritems():
            res[1].wait()
            touch(ex,res[0])
            segOut = open(res[0],"r")
            resBedGraph = unique_filename_in()+".sql"
            sqlTr = track( resBedGraph, fields=['start','end','score'],
                           info={'datatype':'quantitative'}, chrmeta=assembly.chrmeta )
            sqlTr.write(_parse_select_frag(segOut),fields=['chr','start','end','score'])
            sqlTr.close()
            segOut.close()
            results[gid][rid].extend([res[0],resBedGraph])
    return results #[countsPerFrag_allBed, countsPerFrag_selectSql, segToFrag_out, segToFrag_sql]

############################################################################
def c4seq_workflow( ex, job, primers_dict, assembly,
                    c4_url=None, script_path='', logfile=sys.stdout, via='lsf' ):
    '''
    Main
    * open the 4C-seq minilims and create execution
    * 0. get/create the library
    * 1. if necessary, calculate the density file from the bam file (mapseq.parallel_density_sql)
    * 2. calculate the count per fragment for each denstiy file with gfminer:score_by_feature to calculate)
    '''

    mapseq_files = job.files
### outputs
    processed = {'lib': {}, 'density': {}, '4cseq': {}}
    processed['4cseq'] = {'density_files' : {},
                          'countsPerFrag' : {},
                          'countsPerFrag_grp' : {},
                          'norm' : {},
                          'norm_grp' : {},
                          'profileCorrection': {},
                          'profileCorrection_grp' : {},
                          'smooth_grp' : {},
                          'domainogram_grp' : {},
                          'bricks2frags' : {}}
                            # was 'smoothFrag': {}, 'domainogram': {}}
    regToExclude = {}
    new_libs=[]
### options
    run_domainogram = {}
    before_profile_correction = {}
    if not job.options.get('viewpoints_chrs',False):
        out_chromosomes = ','.join([ch for ch in assembly.chrnames])
    else:
        out_chromosomes = ','.join([primers_dict.get(group['name'],{}).get('baitcoord').split(':')[0] for gid,group in job.groups.iteritems()])
    print "out_chromosomes=" + out_chromosomes + "\n"
### do it
    for gid, group in job.groups.iteritems():
        run_domainogram[gid] = group.get('run_domainogram',False)
        if isinstance(run_domainogram[gid],basestring):
            run_domainogram[gid] = (run_domainogram[gid].lower() in ['1','true','on','t'])
        before_profile_correction[gid] = group.get('before_profile_correction',False)
        if isinstance(before_profile_correction[gid],basestring):
            before_profile_correction[gid] = (before_profile_correction[gid].lower() in ['1','true','on','t'])
        processed['lib'][gid] = get_libForGrp(ex, group, assembly,
                                              new_libs, gid, c4_url, via=via)
#reffile='/archive/epfl/bbcf/data/DubouleDaan/library_Nla_30bps/library_Nla_30bps_segmentInfos.bed'
        processed['4cseq']['density_files'][gid] = {}
        regToExclude[gid] = primers_dict.get(group['name'],{}).get('regToExclude',"").replace('\r','')

        # if no regToExclude defined, set it as mid_baitCoord +/-5kb
        if len(regToExclude[gid])==0 :
            baitcoord_mid = int(0.5 * (int(primers_dict.get(group['name'],{}).get('baitcoord').split(':')[1].split('-')[0]) + int(primers_dict.get(group['name'],{}).get('baitcoord').split(':')[1].split('-')[1]) ))
            regToExclude[gid] = primers_dict.get(group['name'],{}).get('baitcoord').split(':')[0] + ':' + str(baitcoord_mid-5000) + '-' + str(baitcoord_mid+5000)

        #print(';'.join([k+"="+v for k,v in primers_dict.get(group['name'],{}).iteritems()]))
        print(primers_dict.get(group['name'],{}))
        print "regToExclude["+str(gid)+"]="+regToExclude[gid]
        for rid,run in group['runs'].iteritems():
            libname = mapseq_files[gid][rid]['libname']
            if job.options.get('merge_strands') != 0 or not('wig' in mapseq_files[gid][rid]):
                density_file=parallel_density_sql( ex, mapseq_files[gid][rid]['bam'],
                                                   assembly.chrmeta,
                                                   nreads=mapseq_files[gid][rid]['stats']["total"],
                                                   merge=0,
                                                   read_extension=mapseq_files[gid][rid]['stats']['read_length'],
                                                   convert=False,
                                                   via=via )
                density_file += "merged.sql"
                ex.add( density_file,
                        description=set_file_descr("density_file_"+libname+".sql",
                                                   groupId=gid,step="density",type="sql",view='admin',gdv="1") )
            else:
                density_file = mapseq_files[gid][rid]['wig']['merged']
            #density_files.append(density_file)
            processed['4cseq']['density_files'][gid][rid]=density_file

        # back to grp level!
        # not anymore:
        # processed['density'][gid] = merge_sql(ex, density_files, via=via)

    processed['4cseq']['countsPerFrag'] = density_to_countsPerFrag( ex, processed, job.groups, assembly, regToExclude, script_path, via )
    ## access per gid+rid

    futures_norm = {}
    countsPerFrags_bedGraph = {}
    futures_merged_raw = {}
    for gid, group in job.groups.iteritems():
        futures_norm[gid] = {}
        countsPerFrags_bedGraph[gid] = {}
        processed['4cseq']['norm'][gid] = {}
        for rid,run in group['runs'].iteritems():
            normfile = unique_filename_in()
            touch(ex, normfile)
            resfile = unique_filename_in()+".bedGraph"
            resfiles = processed['4cseq']['countsPerFrag'][gid][rid] # _all.sql
            convert(resfiles[3],resfile)
            countsPerFrags_bedGraph[gid][rid] = resfile

            print "call normFrags: infiles="+resfile+", normfile="+normfile+"baitCoord="+primers_dict[group['name']]['baitcoord']+", sizeExt=1000000, name="+ group['name']+"rep_"+str(rid) + "regToExclude="+regToExclude[gid]+"\n"
            futures_norm[gid][rid] = normFrags.nonblocking( ex, resfile, normfile, baitCoord=primers_dict[group['name']]['baitcoord'], sizeExt=1000000, name=group['name']+"rep_"+str(rid) ,regToExclude=regToExclude[gid], script_path=script_path, via=via )
            processed['4cseq']['norm'][gid][rid] = normfile

        if len(group) > 1:
            ## merge replicates before normalisation.
            mergefile = unique_filename_in()
            touch(ex, mergefile)
            titleName=group['name']+"_raw_mergedRep"
            print "gid="+group['name']
            print "call mergeRep for replicates before normalisation: infiles="+",".join([res_rid for rid,res_rid in countsPerFrags_bedGraph[gid].iteritems()])+", mergedfile="+mergefile+", regToExclude="+regToExclude[gid]+"\n"
            futures_merged_raw[gid] = mergeRep.nonblocking( ex, ",".join([res_rid for rid,res_rid in countsPerFrags_bedGraph[gid].iteritems()]), mergefile, regToExclude[gid], name=titleName, script_path=script_path, via=via , memory= 8)
            processed['4cseq']['countsPerFrag_grp'][gid] = mergefile
        else:
            futures_merged_raw[gid] = None
            processed['4cseq']['countsPerFrag_grp'][gid] = countsPerFrags_bedGraph[gid][0] #if no replicates, then the file we want is the 1st one

    print "***** profile correction / sample + merge normalised data"
    futures_merged = {} # per gid
    futures_profcor = {} # per gid, per rid
    for gid, group in job.groups.iteritems():
        ## run profile correction per run then merge them
        futures_profcor[gid] = {}
        processed['4cseq']['profileCorrection'][gid] = {}
        for rid, run in group['runs'].iteritems():
            # wait for normalisation of all replicates to be finished
            futures_norm[gid][rid].wait() ## normalised files, per grp, per rep
            normfile = processed['4cseq']['norm'][gid][rid]
            file1 = unique_filename_in() #track file
            touch(ex,file1)
            file2 = unique_filename_in() #report file
            touch(ex,file2)
            file3 = unique_filename_in() #table file
            touch(ex, file3)
            print "call profileCorrection: normfile="+normfile+", baitCoord="+primers_dict[group['name']]['baitcoord']+", name="+group['name']+", file1="+file1+", file2="+file2+", file3= "+file3+"\n"
            futures_profcor[gid][rid] = profileCorrection.nonblocking( ex, normfile,
                                        primers_dict[group['name']]['baitcoord'],
                                        group['name'], file1, file2, file3, script_path,
                                        via=via )
            processed['4cseq']['profileCorrection'][gid][rid] = [file1, file2, file3]

        ## merge replicates before profile correction. Needs all normalisation for the given grp to be finished, this is why it comes after the rid loop.
        if len(group)>1:
            mergefile = unique_filename_in()
            touch(ex, mergefile)
            titleName=group['name']+"_norm_mergedRep"
            print "gid="+group['name']
            print "call mergeRep: infiles="+",".join([res_rid for rid,res_rid in processed['4cseq']['norm'][gid].iteritems()])+", mergedfile="+mergefile+", regToExclude="+regToExclude[gid]+"\n"
            futures_merged[gid] = mergeRep.nonblocking( ex, ",".join([res_rid for rid,res_rid in processed['4cseq']['norm'][gid].iteritems()]), mergefile, regToExclude[gid], name=titleName, script_path=script_path, via=via , memory= 8)
            processed['4cseq']['norm_grp'][gid] = mergefile
        else:
            futures_merged[gid] = None
            processed['4cseq']['norm_grp'][gid] = processed['4cseq']['norm'][gid][0] ##if no replicates, then the file we want is the 1st one

    print "***** merge profile corrected data"
    futures_profcor_merged = {} # per gid
    for gid, group in job.groups.iteritems():
        processed['4cseq']['profileCorrection_grp'][gid] = {}
        for rid, run in group['runs'].iteritems():
            futures_profcor[gid][rid].wait()   ## wait for ProfileCorrection to be finished

        ## merge replicates after profile correction
        if len(group)>1:
            mergefile = unique_filename_in()
            touch(ex, mergefile)
            titleName=group['name']+"_ProfCor_mergedRep"
            pcfiles = [ processed['4cseq']['profileCorrection'][gid][rid][0] for rid,res_rid in processed['4cseq']['profileCorrection'][gid].iteritems()]
            print "call mergeRep (for PC tables): infiles="+",".join(pcfiles)+", mergedfile="+mergefile+", regToExclude="+regToExclude[gid]+"\n"
            futures_profcor_merged[gid] = mergeRep.nonblocking( ex, ",".join(pcfiles), mergefile, regToExclude[gid], name=titleName, script_path=script_path, via=via , memory= 8)
            processed['4cseq']['profileCorrection_grp'][gid] = mergefile
        else:
            futures_profcor_merged[gid] = None
            processed['4cseq']['profileCorrection_grp'][gid] = processed['4cseq']['profileCorrection'][gid][0] ##if no replicates, then the file we want is the 1st one


    print "***** smooth data"
    futures_smoothed = {}
    for gid, group in job.groups.iteritems():
        file1 = unique_filename_in()
        touch(ex,file1)
        file2 = unique_filename_in()
        touch(ex, file2)
        file3 = unique_filename_in()
        touch(ex, file3)
        nFragsPerWin = group['window_size']
        futures_merged_raw[gid].wait() ## wait for merging of raw_grp to be completed
        futures_smoothed[gid] = ( smoothFragFile.nonblocking( ex, processed['4cseq']['countsPerFrag_grp'][gid], nFragsPerWin, group['name'],
                                                    file1, regToExclude[gid], script_path=script_path, via=via, memory=6 ), )
        futures_merged[gid].wait() ## wait for merging of norm_grp to be completed
        futures_smoothed[gid] += ( smoothFragFile.nonblocking( ex, processed['4cseq']['norm_grp'][gid], nFragsPerWin, group['name']+"_norm",
                                                    file2, regToExclude[gid], script_path=script_path, via=via, memory=6 ), )
        futures_profcor_merged[gid].wait() # wait for the merging of profile corrected data to be done
        futures_smoothed[gid] += ( smoothFragFile.nonblocking( ex, processed['4cseq']['profileCorrection_grp'][gid], nFragsPerWin, group['name']+"_fromProfileCorrected",
                                                    file3, regToExclude[gid], script_path=script_path, via=via, memory=6 ), )
        processed['4cseq']['smooth_grp'][gid] = [file1,file2,file3] #[smoothed_file_before_Norm, smoothed file before PC, smoothed file after PC]

    print "***** Domainograms"
    futures_domainograms = {}
    for gid, group in job.groups.iteritems():
        grName = job.groups[gid]['name']
        if run_domainogram[gid]:
            regCoord = regToExclude[gid] or primers_dict[grName]['baitcoord']
            if before_profile_correction[gid]:
               futures_domainograms[gid] = runDomainogram.nonblocking( ex, processed['4cseq']['norm_grp'][gid],
                                                                            grName, regCoord=regCoord, skip=1,
                                                                            script_path="/scratch/cluster/monthly/mleleu/tmp/share/", via=via, memory=15 )
            else:
                futures_domainograms[gid] = runDomainogram.nonblocking( ex, processed['4cseq']['profileCorrection_grp'][gid],
                                                                            grName, regCoord=regCoord.split(':')[0], skip=1,
                                                                            script_path="/scratch/cluster/monthly/mleleu/tmp/share/", via=via, memory=15 )

    ## prepare tar files for domainogram results (if any)
    ## and create "BRICKS to frags" files
    print "***** BRICKS to Frags"
    futures_BRICKS2Frags = {}
    for gid, f in futures_domainograms.iteritems():
        if run_domainogram[gid]: # if domainogram has been run
            resFiles = []
            logFile = f.wait()
            start = False
            tarname = job.groups[gid]['name']+"_domainogram.tar.gz"
            res_tar = tarfile.open(tarname, "w:gz")
            futures_BRICKS2Frags[gid] = []
            processed['4cseq']['bricks2frags'][gid] = []
            if logFile is None: continue
            with open(logFile) as f:
                for s in f:
                    s = s.strip()
                    if '####resfiles####' in s:
                        start = True
                    elif start and "RData" not in s:
                        resFiles.append(s)
                        res_tar.add(s)
                    if start and "foundBRICKS" in s:
                        bricks2fragsfile = unique_filename_in()+".bedGraph"
                        touch(ex, bricks2fragsfile)
                        futures_BRICKS2Frags[gid] += [ BRICKSToFrag.nonblocking(ex, s, processed['4cseq']['norm_grp'][gid], bricks2fragsfile, script_path=script_path, via=via, memory=4 ) ]
                        processed['4cseq']['bricks2frags'][gid] += [ bricks2fragsfile ]
            res_tar.close()
            processed['4cseq']['domainogram_grp'][gid] = resFiles + [tarname]




############### prepare tables for global results
    print "***** combine results into tables "
    allNames=[]
    allFiles=[]
    allRegToExclude=[]
    for gid, group in job.groups.iteritems():
        for rid,run in group['runs'].iteritems():
            allNames += [ group['name']+"_rep"+str(rid)+"_norm", group['name']+"_rep"+str(rid)+"_fit" ]
            allFiles += [ processed['4cseq']['profileCorrection'][gid][rid][2] ]
            allRegToExclude += [ regToExclude[gid] ]
    tablePC=unique_filename_in()+".txt"
    print("***will call makeTable with:")
    print(",".join(allFiles))
    print("resfile="+tablePC)
    print(",".join(allNames))
    touch(ex,tablePC)

    #regToExclude[gid]

    futures_tables = (makeTable.nonblocking(ex, ",".join(allFiles), tablePC, ",".join(allNames), idCols="4,5", all_regToExclude=','.join(allRegToExclude), script_path=script_path, via=via, memory=8 ), )

    # wait for all smoothing to be done
    for gid, fg in futures_smoothed.iteritems():
        for f in fg: f.wait()

    ## make Table raw/smoothed_raw
    print("** make Table raw/smoothed_raw")
    allNames=[]
    allFiles=[]
    allRegToExclude=[]
    for gid, group in job.groups.iteritems():
        futures_merged_raw[gid].wait()
        allNames += [ group['name']+"_raw", group['name']+"_rawSmoothed" ]
        allFiles += [ processed['4cseq']['countsPerFrag_grp'][gid], processed['4cseq']['smooth_grp'][gid][0] ]
        allRegToExclude += [ 'NA', regToExclude[gid] ]

    tableSmoothedRaw_grp=unique_filename_in()+".txt"
    touch(ex,tableSmoothedRaw_grp)
    futures_tables += (makeTable.nonblocking(ex, ",".join(allFiles), tableSmoothedRaw_grp, ",".join(allNames), idCols="4", out_chromosomes = out_chromosomes, all_regToExclude=','.join(allRegToExclude), script_path=script_path, via=via, memory=8 ), )

    ## make Table norm/smoothed_norm before PC
    print("** make Table norm/smoothed_norm befor PC")
    allNames=[]
    allFiles=[]
    allRegToExclude=[]
    for gid, group in job.groups.iteritems():
        allNames += [ group['name']+"_norm", group['name']+"_smoothed" ]
        allFiles += [ processed['4cseq']['norm_grp'][gid], processed['4cseq']['smooth_grp'][gid][1] ]
        allRegToExclude += [ regToExclude[gid], regToExclude[gid] ]

    tableSmoothed_grp=unique_filename_in()+".txt"
    touch(ex,tableSmoothed_grp)
    futures_tables += (makeTable.nonblocking(ex, ",".join(allFiles), tableSmoothed_grp, ",".join(allNames), idCols="4", out_chromosomes = out_chromosomes, all_regToExclude=','.join(allRegToExclude), script_path=script_path, via=via, memory=8 ), )

    ## make Table norm/smoothed_norm after PC
    print("** make Table norm/smoothed_norm after PC")
    allNames=[]
    allFiles=[]
    allRegToExclude=[]
    for gid, group in job.groups.iteritems():
        allNames += [ group['name']+"_normPC", group['name']+"_smoothedPC" ]
        allFiles += [ processed['4cseq']['profileCorrection_grp'][gid], processed['4cseq']['smooth_grp'][gid][2] ]
        allRegToExclude += [ regToExclude[gid], regToExclude[gid] ]

    tableSmoothedPC_grp=unique_filename_in()+".txt"
    touch(ex,tableSmoothedPC_grp)
    futures_tables += (makeTable.nonblocking(ex, ",".join(allFiles), tableSmoothedPC_grp, ",".join(allNames), idCols="4", out_chromosomes = out_chromosomes, all_regToExclude=','.join(allRegToExclude), script_path=script_path, via=via, memory=8 ), )

    ## combine BRICKS2Frags files
    allNames=[]
    allFiles=[]
    for gid, fg in futures_BRICKS2Frags.iteritems():
        for f in fg: f.wait()
        allNames += [ job.groups[gid]['name']+"_BRICKSpval" ]
        cat_bricks2frags = unique_filename_in()+".txt"
        print ','.join(processed['4cseq']['bricks2frags'][gid])
        cat(processed['4cseq']['bricks2frags'][gid],out=cat_bricks2frags)
        allFiles += [ cat_bricks2frags ]

    for gid, fg in futures_smoothed.iteritems():
        for f in fg: f.wait()

    tableBRICKS2Frags = unique_filename_in()+".txt"
    touch(ex,tableBRICKS2Frags)
    futures_tables += (makeTable.nonblocking(ex, ",".join(allFiles), tableBRICKS2Frags, ",".join(allNames), idCols="4", out_chromosomes = out_chromosomes, defVal="NA", script_path=script_path, via=via, memory=8 ), )


    for f in futures_tables: f.wait()


################ Add everything to minilims below!
    step = "density"
    for gid in processed['4cseq']['density_files'].keys():
        for rid, sql in processed['4cseq']['density_files'][gid].iteritems():
            fname = "density_file_"+job.groups[gid]['name']+"_merged_rep"+str(rid)
            ex.add( sql, description=set_file_descr( fname+".sql",
                                                 groupId=gid,step=step,type="sql",gdv="1" ) )
            wig = unique_filename_in()+".bw"
            convert( sql, wig )
            ex.add( wig, description=set_file_descr( fname+".bw",
                                                 groupId=gid,step=step,type="bigWig",ucsc="1") )
    step = "counts_per_frag" #was _norm_counts_per_frags # before normalisation process, per replicate
    for gid in processed['4cseq']['countsPerFrag'].keys():
        for rid, resfiles in processed['4cseq']['countsPerFrag'][gid].iteritems():
            fname = "meanScorePerFeature_"+job.groups[gid]['name']+"_rep"+str(rid)
            ex.add( resfiles[1], description=set_file_descr( fname+".sql",
                                                             groupId=gid,step=step,type="sql",view="admin",gdv='1'))
            #gzipfile(ex,resfiles[0])
            #ex.add( resfiles[0]+".gz", description=set_file_descr( fname+".bed.gz",
            #                                                       groupId=gid,step=step,type="bed",view="admin" ))
            fname = "segToFrag_"+job.groups[gid]['name']+"_rep"+str(rid)
            ex.add( resfiles[3], description=set_file_descr( fname+"_all.sql",
                                                             groupId=gid,step=step,type="sql",
                                                             comment="all informative frags - null included" ))
            trsql = track(resfiles[3])
            bwig = unique_filename_in()+".bw"
            trwig = track(bwig,chrmeta=trsql.chrmeta)
            trwig.write(trsql.read(fields=['chr','start','end','score'],
                                   selection={'score':(0.01,sys.maxint)}))
            trwig.close()
            ex.add( bwig, set_file_descr(fname+".bw",groupId=gid,step=step,type="bigWig",ucsc='1'))
        ## add segToFrags before normalisation
        futures_merged_raw[gid].wait()
        trbedgraph = track(removeNA(processed['4cseq']['countsPerFrag_grp'][gid]),format='bedgraph')
        bwig = unique_filename_in()+".bw"
        trwig = track(bwig,chrmeta=assembly.chrmeta)
        trwig.write(trbedgraph.read(fields=['chr','start','end','score'],
                               selection={'score':(0.01,sys.maxint)}))
        trwig.close()
        fname = "segToFrag_"+job.groups[gid]['name']
        ex.add( bwig, description=set_file_descr( fname+".bw",
                                                             groupId=gid,step=step,type="bigWig",
                                                             comment="segToFrag file before normalisation" ))

    step = "norm_counts_per_frags"  # after new normalisation process, combined replicates
    for gid, resfile in processed['4cseq']['norm_grp'].iteritems():
        fname = "normalised_scorePerFeature_"+job.groups[gid]['name']
        gzipfile(ex,resfile)
        ex.add( resfile+".gz", description=set_file_descr( fname+".bedGraph.gz", groupId=gid,step=step, type="bedGraph",ucsc='1'))
    # norm files, per replicates (might be removed)
    for gid, dict_gid in processed['4cseq']['norm'].iteritems():
        for rid, resfile in dict_gid.iteritems():
            fname = "normalised_scorePerFeature_"+job.groups[gid]['name']+"_rep"+str(rid)
            gzipfile(ex,resfile)
            ex.add(resfile+".gz",
                    description=set_file_descr(fname+".bedGraph.gz",groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))
    step = "profile_correction" # Profile corrected data, combined replicates
    for gid, profileCorrectedFile in processed['4cseq']['profileCorrection_grp'].iteritems():
        fname = "segToFrag_"+job.groups[gid]['name']+"_profileCorrected"
        gzipfile(ex,profileCorrectedFile)
        ex.add( profileCorrectedFile+".gz",
                description=set_file_descr(fname+".bedGraph.gz",groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))
    # Profile corrected, per replicate (might be removed)
    for gid, dict_gid in processed['4cseq']['profileCorrection'].iteritems():
        for rid, resfiles in dict_gid.iteritems():
    #        profileCorrectedFile = resfiles[0]
            reportProfileCorrection = resfiles[1]
            fname = "segToFrag_"+job.groups[gid]['name']+"_profileCorrected_rep"+str(rid)
    #        gzipfile(ex,profileCorrectedFile)
     #       ex.add( profileCorrectedFile+".gz",
      #              description=set_file_descr(fname+".bedGraph.gz",groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))
            ex.add( reportProfileCorrection, description=set_file_descr(fname+".pdf",
                                                                    groupId=gid,step=step,type="pdf"))
    step = "smoothing"
    for gid, resfiles in processed['4cseq']['smooth_grp'].iteritems():
        rawSmoothFile = resfiles[0]
        smoothFile = resfiles[1]
        afterProfileCorrection = resfiles[2]
        nFrags = str(job.groups[gid]['window_size'])
        ## smoothed file before normalisation
        fname = "segToFrag_"+job.groups[gid]['name']+"_smoothed_"+nFrags+"FragsPerWin.bedGraph.gz"
        gzipfile(ex,rawSmoothFile)
        ex.add(rawSmoothFile+".gz",
               description=set_file_descr(fname,groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))
        ## smoothed file after normalisation, before Profile correction
        fname = "segToFrag_"+job.groups[gid]['name']+"_norm_smoothed_"+nFrags+"FragsPerWin.bedGraph.gz"
        gzipfile(ex,smoothFile)
        ex.add(smoothFile+".gz",
               description=set_file_descr(fname,groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))
        ## smoothed file after normalisation, after Profile correction
        fname = "segToFrag_"+job.groups[gid]['name']+"_profileCorrected_smoothed_"+nFrags+"FragsPerWin.bedGraph.gz"
        gzipfile(ex,afterProfileCorrection)
        ex.add(afterProfileCorrection+".gz",
               description=set_file_descr(fname,groupId=gid,step=step,type="bedGraph",ucsc='1',gdv='1'))

    step = "domainograms"
    for gid, resfiles in processed['4cseq']['domainogram_grp'].iteritems():
        tarFile = resfiles.pop()
        fname = job.groups[gid]['name']+"_domainogram.tar.gz"
        ex.add(tarFile, description=set_file_descr(fname,
                                                   groupId=gid,step=step,type="tgz"))
        for s in resfiles:
            if s[-8:] == "bedGraph":
                gzipfile(ex,s)
                s += ".gz"
                ex.add( s, description=set_file_descr( s, groupId=gid,step=step,type="bedGraph",ucsc="1",gdv="1"))

    step = "combined_results"
    gzipfile(ex,tableSmoothedRaw_grp)
    ex.add(tableSmoothedRaw_grp+".gz", description=set_file_descr("table_segToFrags_smoothed_combined_replicates.txt.gz",step=step,type="txt"))

    gzipfile(ex,tableSmoothed_grp)
    ex.add(tableSmoothed_grp+".gz", description=set_file_descr("table_normalised_smoothed_combined_replicates.txt.gz",step=step,type="txt"))

    gzipfile(ex,tableSmoothedPC_grp)
    ex.add(tableSmoothedPC_grp+".gz", description=set_file_descr("table_profileCorrected_smoothed_combined_replicates.txt.gz",step=step,type="txt"))

    gzipfile(ex,tablePC)
    ex.add(tablePC+".gz", description=set_file_descr("table_normalised_fit_per_replicates.txt.gz",step=step,type="txt"))

    gzipfile(ex,tableBRICKS2Frags)
    ex.add(tableBRICKS2Frags+".gz", description=set_file_descr("table_frags_in_BRICKS_combined_replicates.txt.gz",step=step,type="txt"))

    return processed

