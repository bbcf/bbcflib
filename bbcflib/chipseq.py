import sqlite3
from bein import *
from bein.util import *
import gMiner as gm

############ gMiner operations ############
def merge_sql( ex, sqls, ids, description="merged.sql" ):
    """Run gMiner merge_score function on a set of sql files
    """
    out = unique_filename_in()
    req = "[gMiner]\nversion=0.1.3\n"
    req += "\n".join(['track'+str(i+1)+'='+sqls[i]+"\ntrack"+str(i+1)+'_name="'+str(ids[i])+'"' for i in range(len(sqls))])
    req += '''
operation_type=genomic_manip 
manipulation=merge_scores
output_location='''+out
    job = gm.gmRequest(req)
    error, result, type = job.prepare()
    if error != 200: 
        raise ValueError(result)
    error, result, type = job.run()
    if error != 200: 
        raise RuntimeError(result)
    ex.add( out, description="sql:"+description )
    return out

############################################################

 
############ Peaks and annotation ############
@program
def macs( read_length, genome_size, bamfile, ctrlbam=None, shift=80 ):
    outname = unique_filename_in()
    macs_args = ["macs14","-t",bamfile]
    if ctrlbam != None:
        macs_args += ["-c",ctrlbam]
    macs_args += ["-n",outname,"-f","BAM",
                  "-g",str(genome_size),"-s",str(read_length),
                  "--nomodel","-m","5,50","-p",".001","--bw=200"]
    if shift>0:
        macs_args += ["--shiftsize="+str(shift)]
    return {"arguments": macs_args, "return_value": outname}

def add_macs_results( ex, chromosomes, read_length, genome_size, bamfile,
                      ctrlbam=[None], name=None, shift=80, alias=None ):
    """Adds macs result files to the repository.
    """
    if not(isinstance(bamfile,list)):
        bamfile = [bamfile]
    if not(isinstance(ctrlbam,list)):
        ctrlbam = [ctrlbam]
    futures = {}
    i = 0
    for bam in bamfile:
        n = name['tests'][i]
        i += 1
        j = 0
        for cam in ctrlbam:
            m = name['controls'][j]
            if m == None:
                n = (n,)
            else:
                n = (n,m)
            j += 1
            futures[n] = macs.lsf( ex, read_length, genome_size, bam, cam, shift )
    prefixes = dict((n,f.wait()) for n,f in futures.iteritems())
    for n,p in prefixes.iteritems():
        description = "_vs_".join(n)
        touch( ex, p )
        ex.add( p, description="none:"+description, alias=alias )
        ex.add( p+"_peaks.xls", description="macs:"+description+"peaks.xls",
                associate_to_filename=p, template='%s_peaks.xls' )
        ex.add( p+"_peaks.bed", description="macs:"+description+"peaks.bed",
                associate_to_filename=p, template='%s_peaks.bed' )
        ex.add( p+"_summits.bed", description="macs:"+description+"summits.bed",
                associate_to_filename=p, template='%s_summits.bed' )
        if len(n)>1:
            ex.add( p+"_negative_peaks.xls", description="macs:"+description+"negative_peaks.xls",
                    associate_to_filename=p, template='%s_negative_peaks.xls' )
    return prefixes

############################################################ 

def cat(files):
    if len(files) > 1:
        out = unique_filename_in()
        with open(out,"w") as f1:
            for inf in files: 
                with open(inf,"r") as f2:
                    [f1.write(l) for l in f2]
                    f2.close()
            f1.close()
    elif len(files) == 1:
        out = files[0]
    else:
        out = None
    return out

@program
def bam_to_density( bamfile, chromosome_accession, chromosome_name, output,
                    nreads=1, merge=-1, read_length=-1, convert=True, sql=False ):
    """Runs the bam2wig program on a bam file and 
    normalizes for the total number of reads
    provided as argument nreads. 

    Returns the name of the output wig or sql file(s) (if 'sql' is True).

    Use 'convert'=False if the bam already uses chromosome names.
    """
    b2w_args = ["-w",str(nreads),"-s",bamfile,"-o",output]
    if convert:
        b2w_args += ["-a",chromosome_accession,"-n",chromosome_name]
    else:
        b2w_args += ["-a",chromosome_name]
    if merge>=0:
        b2w_args += ["-p",str(merge)]
    if read_length>0:
        b2w_args += ["-q",str(read_length)]
    if sql:
        b2w_args += ["-d"]
        if merge<0:
            files = [output+"fwd",output+"rev"]
        else:
            files = [output+"merged"]
    else:
        if merge<0:
            b2w_args += ["-6"]
        files = output
    return {"arguments": ["bam2wig"]+b2w_args, "return_value": files}


def parallel_density_wig( ex, bamfile, chromosomes, 
                          nreads=1, merge=-1, read_length=-1, convert=True, 
                          description="", alias=None ):
    """Runs 'bam_to_density' in parallel 
    for every chromosome in the 'chromosomes' list.
    Returns produces a single text wig file.
    """
    futures = [bam_to_density.lsf(ex,bamfile,k,v['name'],unique_filename_in(),
                                  nreads,merge,read_length,convert,False)
               for k,v in chromosomes.iteritems()]
    results = []
    for f in futures:
        try: 
            results.append(f.wait())
        except ProgramFailed:
            pass
    output = cat(results)
    ex.add( output, description=description, alias=alias )
    return output


def parallel_density_sql( ex, bamfile, chromosomes, 
                          nreads=1, merge=-1, read_length=-1, convert=True, 
                          description="", alias=None ):
    """Runs 'bam_to_density' for every chromosome in the 'chromosomes' list.
    Returns one or two sqlite files depending if 'merge'>=0 or 'merge'<0.
    """
    output = unique_filename_in()
    touch(ex,output)
    if merge<0:
        suffix = ['fwd','rev']
    else:
        suffix = ['merged']
    for s in suffix:
        conn1 = sqlite3.connect( output+s )
        conn1.execute('create table chrNames (name text, length integer)')
        conn1.execute('create table attributes (key text, value text)')
        conn1.execute('insert into attributes (key,value) values ("datatype","quantitative")')
        vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
        conn1.executemany('insert into chrNames (name,length) values (?,?)',vals)
        [conn1.execute('create table '+v['name']+' (start integer, end integer, score real)') 
         for v in chromosomes.values()]
        conn1.commit()
        conn1.close()
    results = []
    for k,v in chromosomes.iteritems():
        future = bam_to_density.lsf(ex,bamfile,k,v['name'],output,
                                    nreads,merge,read_length,convert,True)
        try: 
            results.append(future.wait())
        except ProgramFailed:
            pass
    ex.add( output, description="none:"+description, alias=alias )
    for outf in results[0]:
        suffix = outf[len(output):len(outf)]
        ex.add( outf, description="sql:"+description+"_"+suffix,
                associate_to_filename=output, template='%s'+suffix )
    return output

def get_chipseq_files( hts_key, minilims, lib_root, hts_url, gdv_url,
                       mapseq_url=None ):
    def path_from_lib( M, id, root ):
        return os.path.relpath(M.path_to_file(id),root)
    files = {}
    bamfiles = {}
    htss = frontend.Frontend( url=hts_url )
    job = htss.job( hts_key )
    try: 
        exid = max(minilims.search_executions(with_text=hts_key))
    except ValueError, v:
        raise ValueError("No execution with key "+hts_key)
    if job.options['select_source'] == 'mapseq_key':
        files['url'] = {mapseq_url+job.options['mapseq_key']+"/get_results": "Mapseq results"}
    gid2name = {}
    for gid,group in job.groups.iteritems():
        name = re.sub(r'\s+','_',group['name'])
        gid2name[gid] = name
        if job.options['select_source'] != 'mapseq_key':
            for rid,run in group['runs'].iteritems():
                both = minilims.search_files(with_text=str(rid)+" filtered bam", 
                                             source=('execution',exid))
                file = "_".join([name,run['machine'],str(run['run']),str(run['lane'])])+".bam"
                bamfiles[path_from_lib(minilims,both.pop(),lib_root)] = file+".bai"
                bamfiles[path_from_lib(minilims,both.pop(),lib_root)] = file
            files["bam"] = bamfiles
            names = "_".join([g['name'] for g in job.groups.values()])
            fid = minilims.search_files(with_text="mapping pdf report",
                                        source=('execution',exid)).pop()
            files["pdf"] = {path_from_lib(minilims,fid,lib_root):
                            names+".pdf"}
        sql_text = str(gid)+" group density"
        files["sql"] = {}
        if job.options['merge_strands']>=0:
            sql_fid = minilims.search_files(with_text=sql_text+" merged", 
                                            source=('execution',exid)).pop()
            files["sql"][path_from_lib(minilims,sql_fid,lib_root)] = name+"_merged.sql"
        else:
            sql_fid = minilims.search_files(with_text=sql_text+" fwd", 
                                            source=('execution',exid)).pop()
            files["sql"][path_from_lib(minilims,sql_fid,lib_root)] = name+"_fwd.sql"
            sql_fid = minilims.search_files(with_text=sql_text+" rev", 
                                            source=('execution',exid)).pop()
            files["sql"][path_from_lib(minilims,sql_fid,lib_root)] = name+"_rev.sql"
    if job.options['peak_calling']:
        files['macs'] = {}
        for fid in minilims.search_files(with_text=" macs", source=('execution',exid)):
            fd = minilims.fetch_file(fid)['description'].split(' ')
            fn = "_vs_".join([gid2name[int(x)] for x in fd[0].split(':')])
            files['macs'][path_from_lib(minilims,fid,lib_root)] = fn+"_".join(fd[1:])
    files['url'] = {gdv_url+"/log/htschipseq/"+hts_key: 'GDV view'}
    return files


