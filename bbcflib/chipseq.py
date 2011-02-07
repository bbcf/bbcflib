import sqlite3
from bein import *
from bein.util import *
import gMiner as gm

############ gMiner operations ############
def merge_sql( ex, sqls, ids, description="merged sql" ):
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
    ex.add( out, description=description )
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
                  "--nomodel","--shiftsize="+str(shift),"-m","5,50",
                  "-p",".001","--bw=200"]
    return {"arguments": macs_args, "return_value": outname}

def add_macs_results( ex, chromosomes, read_length, genome_size, bamfile,
                      ctrlbam=[None], name=None, shift=80, alias=None ):
    """Adds macs result files to the repository.
    """
    genomsize = sum([c['length'] for c in chromosomes.values()])
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
        description = ":".join(n)+" macs"
        touch( ex, p )
        ex.add( p, description=description, alias=alias )
        ex.add( p+"_peaks.xls", description=description+" positive peaks",
                associate_to_filename=p, template='%s_peaks.xls' )
        ex.add( p+"_peaks.bed", description=description+" peaks bed file",
                associate_to_filename=p, template='%s_peaks.bed' )
        ex.add( p+"_summits.bed", description=description+" summits bed file",
                associate_to_filename=p, template='%s_summits.bed' )
        if len(n)>1:
            ex.add( p+"_negative_peaks.xls", description=description+" negative peaks",
                    associate_to_filename=p, template='%s_negative_peaks.xls' )
    return prefixes

############################################################ 

def cat(files):
    out = unique_filename_in()
    with open(out,"w") as f1:
        for inf in files: 
            with open(inf,"r") as f2:
                [f1.write(l) for l in f2]
                f2.close()
        f1.close()
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


def parallel_density( ex, bamfile, chromosomes, 
                      nreads=1, merge=-1, read_length=-1, convert=True, sql=False,
                      description="", alias=None ):
    """Runs 'bam_to_density' in parallel for every chromosome in the 'chromosomes' list.
    """
    output = unique_filename_in()
    futures = [bam_to_density.lsf(ex,bamfile,k,v["name"],output,nreads,merge,read_length,convert,sql)
               for k,v in chromosomes.iteritems()]
    results = []
    for f in futures:
        try: 
            results.append(f.wait())
        except ProgramFailed:
            pass
    if sql:
        touch( ex, output )
        ex.add( output, description=description, alias=alias )
        for outf in results[0]:
            connection = sqlite3.connect( outf )
            connection.execute('create table chrNames (name text, length integer)')
            connection.execute('create table attributes (key text, value text)')
            connection.execute('insert into attributes (key,value) values ("datatype","quantitative")')
            vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
            connection.executemany('insert into chrNames (name,length) values (?,?)',vals)
            connection.commit()
            connection.close()
            suffix = outf[len(output):len(outf)]
            ex.add( outf, description=description+" "+suffix,
                    associate_to_filename=output, template='%s'+suffix )
    else:
        if len(results) > 1:
            output = cat(results)
        elif len(results) == 1:
            output = results[0]
        else:
            output = None
        ex.add( output, description=description, alias=alias )
    return output
