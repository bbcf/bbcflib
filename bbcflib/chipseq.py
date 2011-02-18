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
                  "--nomodel","-m","5,50","-p",".001","--bw=200"]
    if shift>0:
        macs_args += ["--shiftsize="+str(shift)]
    return {"arguments": macs_args, "return_value": outname}

def add_macs_results( ex, read_length, genome_size, bamfile,
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

@program
def sql_prepare_deconv(sql_prefix,macs_bed,chr_name,chr_length,cutoff,read_length):
    out = unique_filename_in()
    args = [sql_prefix,macs_bed,out,chr_name,str(chr_length),
            str(cutoff),str(read_length)]
    return {"arguments": ["sql_prepare_deconv.py"]+args, 
            "return_value": out}

@program
def run_deconv_r( counts_file, read_length, chr_name, script_path ):
    pdf_file = unique_filename_in()
    output_file = unique_filename_in()
    rargs = [counts_file, pdf_file, str(read_length), 
             chr_name, output_file, script_path]
    return {'arguments': ["R","--vanilla","--slave","-f",
                          os.path.join(script_path,"deconv.R"),
                          "--args"] + rargs,
            'return_value': {'pdf': pdf_file, 'rdata': output_file}}

@program
def sql_finish_deconv(sqlout,rdata):
    return {"arguments": ["sql_finish_deconv.py",rdata,sqlout], 
            "return_value": sqlout}

def run_deconv(ex,sql,macs,chromosomes,read_length,script_path):
    prep_futures = dict((chr['name'],
                         sql_prepare_deconv.lsf(ex,sql,macs,chr['name'],
                                                chr['length'],1500,read_length))
                        for chr in chromosomes.values())
    rdeconv_futures = dict((chr,
                            run_deconv_r.lsf(ex,f.wait(), read_length,
                                             chr, script_path))
                           for chr,f in prep_futures.iteritems())
    rdeconv_out = dict((chr, f.wait()) for chr,f in rdeconv_futures.iteritems())
    pdf_future = join_pdf.lsf(ex,[x['pdf'] for x in rdeconv_out.values()])
    sqlout = unique_filename_in()
    outfiles = {}
    conn = sqlite3.connect( sqlout )
    conn.execute('create table chrNames (name text, length integer)')
    conn.execute('create table attributes (key text, value text)')
    conn.execute('insert into attributes (key,value) values ("datatype","quantitative")')
    vals = [(v['name'],str(v['length'])) for v in chromosomes.values()]
    conn.executemany('insert into chrNames (name,length) values (?,?)',vals)
    [conn.execute('create table '+v['name']+' (start integer, end integer, score real)') 
     for v in chromosomes.values()]
    conn.commit()
    conn.close()
    for fout in rdeconv_out.values():
        f = sql_finish_deconv.lsf(ex,sqlout,fout['rdata'])
        f.wait()
    outfiles['sql'] = sqlout
    outfiles['pdf'] = pdf_future.wait()
    return outfiles
    
############################################################ 
@program
def merge_two_bed(file1,file2):
    out = unique_filename_in()
    args = ['intersectBed','-a',file1,'-b',files2]
    def catch_out(p):
        with open(out,'w') as f:
            [f.write(l) for l in p.stdout]
            f.close()
        return out
    return {"arguments": args, "return_value": catch_out}

def merge_many_bed(ex,files):
    out = files[0]
    for f in files[1:]:
        future = merge_two_bed.lsf(ex,out,f)
        out = future.wait()
    return out

@program
def join_pdf(files):
    out = unique_filename_in()
    gs_args = ['gs','-dBATCH','-dNOPAUSE','-q','-sDEVICE=pdfwrite',
               '-sOutputFile=%s'%out]
    gs_args += files
    return {"arguments": gs_args, "return_value": out}

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
    ex.add( output, description='none:'+description+'.sql', alias=alias )
    for outf in results[0]:
        suffix = outf[len(output):len(outf)]
        ex.add( outf, description='sql:'+description+'_'+suffix+'.sql',
                associate_to_filename=output, template='%s_'+suffix+'.sql' )
    return output

