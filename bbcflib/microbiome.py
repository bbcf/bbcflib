"""
==========================
Module: bbcflib.microbiome
==========================
"""

import sys, pysam
from bbcflib.common import set_file_descr, unique_filename_in
from bein import program

@program
def run_microbiome( options=[], output=None ):
    if output is None: output = unique_filename_in()
    options = [",".join([str(x) for x in o]) if isinstance(o,(list,tuple)) else str(o)
               for o in options]
    return {'arguments': ["run_microbiome.py"]+options+[output], 'return_value': output }


def bam_to_annot_counts( bamfiles, annotations_file, pref_name='', output=None ):
    '''
    Scan each bam file of a list and calculate the corrected counts for each annotation key
    present in the "annotations_file".
    '''
    if output is None: output = unique_filename_in()
    map = {}
    counts = {}
    with open(annotations_file) as f:
        header = f.next().strip('\n').split("\t")
        for line in f:
            s = line.strip('\n').split("\t")
            k = s.pop(0)
            map[k] = s
            counts[k] = 0

    tot = 0
    for bamfile in bamfiles:
        infile = pysam.Samfile( bamfile )
        for read in infile:
            nh = dict(read.tags).get('NH',1)
            if isinstance(nh,basestring): nh = 1
            if nh < 1: continue
            inh = 1.0/nh
            rname = infile.getrname(read.rname)
            if rname in counts:
                counts[rname] += inh
## still increment if not in counts?
            tot += inh
        infile.close()

    with open(output,'w') as out:
        out.write('\t'.join([header[0],'counts_'+pref_name,'%counts_'+pref_name]+header[1:])+'\n')
        for k,v in map.iteritems():
            pc = 100*counts[k]/tot
            out.write('\t'.join([k,"%.2f"%counts[k],"%.3f"%pc]+map[k])+'\n')

    return output

#####################################
def getCountsPerLevel( infile, level=None, output=None ):
    if output is None: output = unique_filename_in()
    counts = {}
    map = {}
    tot = 0
    idColCounts = 1
    name = ''
    with open(infile) as f:
        header = f.next().strip('\n').split('\t')
        try:
            level_idx = header.index(level)
        except:
            raise ValueError("No column corresponds to "+level+" in file "+infile)

        level_top = header.index('Kingdom')
        colrange = range(level_idx,level_top,2*int(level_top>level_idx)-1)
        header_out = [header[n] for n in colrange]
        name = header[idColCounts]
        for line in f:
            s = line.strip('\n').split('\t')
            if len(s) < len(header): s.extend(['']*(len(header)-len(s)))
            tot += float(s[idColCounts])
            counts[s[level_idx]] = counts.get(s[level_idx],0.0)+float(s[idColCounts])
            map[s[level_idx]] = [s[n] for n in colrange]

    with open(output,'w') as out:
        header = [level]+header_out+["counts_"+name,"%counts_"+name]
        out.write("\t".join(header)+"\n")
        for k,v in map.iteritems():
            pc = 100*counts[k]/tot
            curk = k or 'Unnanotated'
            out.write("\t".join([curk]+v+["%.2f"%counts[k],"%.3f"%pc])+"\n")

    return output


##########################################
def combine_counts( counts, idsColsKey, idsColsCounts, output="combined_counts.txt" ):
    if output in [None,'']: output = unique_filename_in()
    all_counts = {}
    infos = {}
    leninfos = 0
    if not isinstance(idsColsKey,(list,tuple)): idsColsKey = [idsColsKey]
    if not isinstance(idsColsCounts,(list,tuple)): idsColsCounts = [idsColsCounts]

    for i,filename in enumerate(counts):
        with open(filename) as f:
            s = f.next().strip('\n').replace("[","").replace("]","").split("\t")
            if i == 0: #1st file: initialization of counts and infos
                _colinfos = [ss for n,ss in enumerate(s) if n not in idsColsKey+idsColsCounts]
                leninfos = len(_colinfos)
                h_infos = '\t'.join(_colinfos)
                h_counts = '\t'.join([s[n] for n in idsColsCounts])
                h_key='\t'.join([s[n] for n in idsColsKey])
            else:
                h_counts += '\t'.join(['']+[s[n] for n in idsColsCounts])
            for line in f:
                s = line.strip('\n').replace("[","").replace("]","").split("\t")
                curKey = '\t'.join([s[n] for n in idsColsKey])
                if i == 0: #1st file: initialization of counts and infos
                    all_counts[curKey] = ['']*len(counts)
                    curInfo = [ss for n,ss in enumerate(s) if n not in idsColsKey+idsColsCounts]
                    if len(curInfo) < leninfos: curInfo.extend(['']*(leninfos-len(curInfo)))
                    infos[curKey] = '\t'.join(curInfo)
                all_counts[curKey][i] = '\t'.join([s[n] for n in idsColsCounts])

    with open(output,'w') as out:
        out.write(h_key+'\t'+h_counts+'\t'+h_infos+'\n')
        for k,v in all_counts.iteritems():
            out.write(k+'\t'+'\t'.join(str(s) for s in all_counts[k])+'\t'+infos.get(k,'')+'\n')

    return(output)

###############################################################
def microbiome_workflow( ex, job, assembly, logfile=sys.stdout, via='lsf' ):
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
        bamfiles = [m['bam'] for m in mapseq_files[gid].values()]
        futures[gid] = run_microbiome.nonblocking(ex, ["bam_to_annot_counts", bamfiles,
                                                       assembly.annotations_path, group_name], 
                                                  via=via)

    # 1.b get counts per Level (Kingdom, Phylum, Class, Order, Family, Genus and Species) (=> 1 file per level / per group)
    step = 'counts'
    for gid, future in futures.iteritems():
        res = future.wait()
        processed['cnts'][gid] = res # group_name + "_counts_annot.txt"
        fname = job.groups[gid]['name']+"_counts_annot.txt"
        ex.add( res, description=set_file_descr( fname, groupId=gid, step=step, type="txt" ) )
        processed['cnts_level'][gid] = [run_microbiome.nonblocking(ex, ["getCountsPerLevel", res, level], via=via) 
                                        for level in levels]

    # 2.a combine counts for all groups (=> 1 combined file)
    files = processed['cnts'].values()
    combined_out = [run_microbiome.nonblocking(ex, ["combine_counts", files, 0, [1,2]], via=via)]

    # 2.b combine counts per level for all groups (=> 1 combined file per Level)
    for n,level in enumerate(levels):
        files = dict([(gid,f[n].wait()) for gid,f in processed['cnts_level'].iteritems()])
        combined_out.append(run_microbiome.nonblocking(ex, ["combine_counts", files.values()]+infosCols.get(level,[0,[1,2]]),
                                                       via=via))
        for gid,f in files.iteritems():
            fname = job.groups[gid]['name']+"_counts_annot_"+level+".txt"
            ex.add( f, description=set_file_descr( fname, groupId=gid, step=step, type="txt" ) )


    step = 'combined'
    ex.add(combined_out[0].wait(), description=set_file_descr( "combined_counts.txt", step=step, type="txt") )
    for nl,level in enumerate(levels):
        ex.add(combined_out[nl+1].wait(), 
               description=set_file_descr( "combined_counts"+level+".txt", step=step, type="txt") )
    return 0

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
