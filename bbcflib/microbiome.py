"""
==========================
Module: bbcflib.microbiome
==========================
"""

def bam_to_annot_counts ( bamfiles, annotations_file, pref_name='' ):
    '''
    Scan each bam file of a list and calculate the corrected counts for each annotation key
    present in the "annotations_file".
    '''

    # read annotations
    map = {}
    i=0
    with open(annotations_file,'r') as f:
        for line in f:
            s = line.strip('\n').split()
            if(i==0):
                header=line
                i = 1
                continue
            map[s[0]] = s
    f.close()

    # get norm counts
    counts = {}
    # initialize to 0
    for k in map.keys():
        counts[k]=0

    # read each bam and update counts
    tot=0
    for bamfile in bamfiles:
        infile = pysam.Samfile( bamfile, "rb" )
        for read in infile:
            nh = dict(read.tags).get('NH',1)
            if isinstance(nh,str): nh=1
            if nh < 1:
                continue
            rname = infile.getrname(read.rname)
            counts[rname] += 1.0/nh
            tot += 1.0/nh
        infile.close()

    resfile = unique_filename_in()
    out = open(resfile,'w')
    out.write(header.split('\t')[0]+'\tcounts_'+pref_name+'\t%counts_'+pref_name+'\t'+'\t'.join(header.split()[1:len(header.split())])+'\n')
    for k,v in map.iteritems():
        if tot > 0.0:
            pc = 100*(counts[k]/tot)
        else:
            pc = 0.0
        out.write(k+'\t'+str(round(float(counts[k]),2))+'\t'+str(round(float(pc),3))+'\t'+'\t'.join(map[k][1:len(map[k])])+'\n')

    out.close()

    return resfile

#####################################
def getCountsPerLevel( infile, level=None ):
    print("get counts per "+level)
    i = 0
    counts = {}
    map = {}
    tot = 0
    idColCounts = 1
    name = ''
    with open(infile,'r') as f:
        for line in f:
            s = line.strip('\n').split('\t')
            if(i==0):
                header=s
                try:
                    level_idx = header.index(level)
                except:
                    print("no column corresponds to "+level+" in file "+infile)
                    break #no column corresponds to the level required - nothing to do
                i = 1
                level_top = header.index('Kingdom')
                name = s[idColCounts]
                continue
            if len(s) < len(header): s=s+['']*(len(header)-len(s)) #make sure there is enough items
            tot += float(s[idColCounts])
            counts[s[level_idx]] = counts.get(s[level_idx],0.0) + float(s[idColCounts])
#            if s[level_idx] in counts.keys():
#                counts[s[level_idx]] += float(s[idColCounts]) # !! s[2] might change
#            else:
#                counts[s[level_idx]] = float(s[idColCounts]) # !! idem
            if level_top < level_idx:
                map[s[level_idx]] = (s[level_top:level_idx])[::-1]
                header_out = (header[level_top:level_idx])[::-1]
            elif level_idx == level_top:
                map[s[level_idx]] = s[level_idx]
                header_out = header[level_idx]
            else:
                map[s[level_idx]] = s[level_idx:level_top]
                header_out = header[level_idx:level_top]
    f.close()

    outfile = os.path.split(infile)[1].split('.')[0]+"_"+level+".txt"
    out = open(outfile,'w')

    if level_idx == level_top :
        out.write(level+"\tcounts_"+name+"\t%counts_"+name+"\n")
    else:
        out.write(level+'\t'+'\t'.join(header_out)+"\tcounts_"+name+"\t%counts_"+name+"\n")
    for k,v in map.iteritems():
        if tot > 0.0:
            pc = 100*(counts[k]/tot)
        else:
            pc = 0.0
        curk = k
        if len(k)==0: curk='Unnanotated'
        if level_idx == level_top :
            out.write(curk+'\t'+str(round(float(counts[k]),0))+'\t'+str(round(float(pc),3))+'\n')
        else:
            out.write(curk+'\t'+'\t'.join(v)+'\t'+str(round(float(counts[k]),0))+'\t'+str(round(float(pc),3))+'\n')
    out.close()

    return outfile


##########################################
def combine_counts(counts, idsColsKey, idsColsCounts, idColsInfos = None, resfile = "combined_counts.txt"):
    all_counts = {}
    infos = {}
    idColsInfos = None
    h_key = ''
    h_counts = ''
    h_infos = ''

    if not isinstance(idsColsKey,list) or len(idsColsKey) < 2 :  # 0 or [0]
        idsColsKey = [ idsColsKey, (idsColsKey+1) ]
    else: # [0,1,2,3]
        idsColsKey = [ idsColsKey[0], (idsColsKey[-1]+1) ] # first and last item

    if not isinstance(idsColsCounts,list) or len(idsColsCounts) < 2 : # [0,1,2,3]
        idsColsCounts = [ idsColsCounts, (idsColsCounts+1) ]
    else:
        idsColsCounts = [ idsColsCounts[0], (idsColsCounts[-1]+1) ] # first and last item

    i = 0
    l = 0
    for filename in counts:
        j = 0
        with open(filename,'r') as f:
            print(filename)
            for line in f:
                line = line.strip('\n').replace("[","").replace("]","")
                s = line.split()
                if i == 0: #1st file: initialization of counts and infos
                    if j == 0: # 1st line of 1st file - init idColsInfos and headers
                        l = len(s)
                        idColsInfos=range(l)
                        [idColsInfos.remove(x) for x in range(idsColsKey[0],idsColsKey[1])+range(idsColsCounts[0],idsColsCounts[1])]
                        if len(idColsInfos) > 0 :
                            #h_infos = '\t'.join(s[idColsInfos[0]:idColsInfos[-1]])
                            h_infos = '\t'.join(s[idColsInfos[0]:(idColsInfos[-1]+1)])
                        h_counts = '\t'.join(s[idsColsCounts[0]:idsColsCounts[1]]) #the +1 is already taken into account during the ini    tialisation of idColsCounts
                        h_key='\t'.join(s[idsColsKey[0]:idsColsKey[1]])
                        j = 1
                    else: # init all_counts for each line
                        curKey = '\t'.join(s[idsColsKey[0]:idsColsKey[1]])
                        all_counts[curKey] = ['']*len(counts)
                        all_counts[curKey][i]='\t'.join(s[idsColsCounts[0]:idsColsCounts[1]])
                        if len(idColsInfos) > 0 :
                            curInfo = s[idColsInfos[0]:(idColsInfos[-1]+1)]
                            if(len(curInfo) < len(idColsInfos)): curInfo = curInfo + ['']*(l-len(curInfo))
                            infos['\t'.join(s[idsColsKey[0]:idsColsKey[1]])] = '\t'.join(curInfo)
                elif j == 0:
                    h_counts = h_counts + '\t' + '\t'.join(s[idsColsCounts[0]:idsColsCounts[1]])
                    j = 1
                else:
                    curKey = '\t'.join(s[idsColsKey[0]:idsColsKey[1]])
                    all_counts[curKey][i]='\t'.join(s[idsColsCounts[0]:idsColsCounts[1]])
            f.close()
            i = i + 1

    if resfile in [None,'']: resfile = unique_filename_in()
    out = open(resfile,'w')
    out.write(h_key+'\t'+h_counts+'\t'+h_infos+'\n')
    for k,v in all_counts.iteritems():
        out.write(k+'\t'+'\t'.join(str(s) for s in all_counts[k])+'\t'+infos.get(k,'')+'\n')
    out.close()

    return(resfile)

###############################################################
@program
def microbiome_workflow( ex, job, assembly,
                         microbiome_url=None, script_path='', logfile=sys.stdout, via='lsf' ):
    '''
    Main
    * 0. retrieve bam files from mapseq job
    *   0.a. merge bam files (=> 1 bam file per group)
    * 1. for each group:
    *   1.a get counts per group (=> 1 file per group)
    *   1.b get counts per Level (Kingdom, Phylum, Class, Order, Family, Genus and Species) (=> 1 file per level / per group)
    * 2. combine counts
    *   2.a combine counts for all groups (=> 1 combined file)
    *   2.b combine counts per level for all groups (=> 1 combined file per Level)
    * 3. generate barplots (=> 1 plot per group + per level + per combined files)
    '''
### params
    #annotations_file = "GG_13_5_otu_annotations.txt"
    #/db/genrep/nr_assemblies/annot_txt/cf6689ddbf5ce2ee47dffcf44c0fc00d2869d5e4.txt
    annotations_file = "/db/genrep/nr_assemblies/annot_txt/"+ assembly.md5 +".txt"
    levels = ['Kingdom','Phylum','Class','Order','Family', 'Genus', 'Species']
    infosCols = {'Kingdom':[0,[1,2]],
                'Phylum':[[0,1],[2,3]],
                'Class':[[0,1,2],[3,4]],
                'Order':[[0,1,2,3],[4,5]],
                'Family':[[0,1,2,3,4],[5,6]]}
### outputs
    processed = {'cnts': {}, 'cnts_level': {}, 'plots': {}}

### do it
    mapseq_files = job.files

    # 1.a get counts per group (=> 1 file per group)
    futures = {}
    for gid, group in job.groups.iteritems():
        group_name = group['name']
        futures[gid] = bam_to_annot_counts.nonblocking(ex, mapseq_files[gid], annotations, group_name, via=via)

    # 1.b get counts per Level (Kingdom, Phylum, Class, Order, Family, Genus and Species) (=> 1 file per level / per group)
    for gid, future in futures.iteritems():
        res = future.wait()
        processed['cnts'][gid] = res # group_name + "_counts_annot.txt"
        processed['cnts_level'][gid] = [getCountsPerLevel.nonblocking(res,level) for level in levels]

    # 2.a combine counts for all groups (=> 1 combined file)
    files = [[processed['cnts'][gid],group['name']] for gid, group in job.groups.iteritems()]
    combine_counts(files, idsColsKey = 0, idsColsCounts = [1,2], idColsInfos = None, resfile="combined_counts.txt") # "combined_counts.txt"

    # 2.b combine counts per level for all groups (=> 1 combined file per Level)
    for n,level in enumerate(levels):
        files = [processed['cnts_level'][gid][n].wait() for gid in processed['cnts_level'].keys()]
        #files = [processed['cnts_level'][gid][n].wait() for gid in processed['cnts_level'].keys()]
        #files = [file for i, file in enumerate(all_levelsfiles) if re.search(level, file)]
        combine_counts(files, idsColsKey = infosCols.get(level)[0], idsColsCounts = infosCols.get(level)[1], idColsInfos = None, resfile="combined_counts_"+level+".txt") # "combined_counts_"+level+".txt"

### add results files to lims
    step='counts'
    for gid, resfile in processed['cnts'].iteritems():
        fname = job.groups[gid]['name']+"_counts_annot.txt"
        ex.add( resfile, description=set_file_descr( fname, groupId=gid, step=step, type="txt" ) )

    for gid, resfiles in processed['cnts_levels'].iteritems():
        for i in range(len(levels)): # 1 file per level
            fname = job.groups[gid]['name']+"_counts_annot_"+ levels[i] +".txt"
            ex.add( resfiles[i], description=set_file_descr( fname, groupId=gid, step=step, type="txt" ) )

    step='combined'
    ex.add("combined_counts.txt", description=set_file_descr( "combined_counts.txt", step=step, type="txt") )
    for level in enumerate(levels):
        ex.add("combined_counts"+ level +".txt", description=set_file_descr( "combined_counts"+ level +".txt", step=step, type="txt") )


    return processed

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
