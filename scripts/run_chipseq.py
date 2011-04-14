from bbcflib import daflims, genrep, frontend, email, gdv, common
from bbcflib.mapseq import *
from bbcflib.chipseq import *
import sys, getopt, os, json

opts = dict(getopt.getopt(sys.argv[1:],"k:d:w:",[])[0])
hts_key = opts['-k']
limspath = opts['-d']
working_dir = opts['-w']
ms_minilims = "/srv/mapseq/public/data/mapseq_minilims"
mapseq_url = ''

os.chdir(working_dir)
M = MiniLIMS( limspath )
gl = use_pickle(M, "global variables")
htss = frontend.Frontend( url=gl["hts_url"] )
job = htss.job( hts_key )
job.options['ucsc_bigwig'] = True
g_rep = genrep.GenRep( gl["genrep_url"], gl["bwt_root"] )
#_ = [M.delete_execution(x) for x in M.search_executions(with_text=hts_key)]
with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
    # should remove this option
    if job.options['select_source'] == 'lims':
        g_rep_assembly = g_rep.assembly( job.assembly_id )
        dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'],
                                          password=gl['lims']['passwd'][loc] ))
                    for loc in gl['lims']['passwd'].keys())
        job.options['discard_pcr_duplicates'] = True
        mapseq_files = map_groups( ex, job, dafl, ex.working_directory, g_rep )
        pdf = add_pdf_stats( ex, mapseq_files,
                             dict((k,v['name']) for k,v in job.groups.iteritems()), 
                             gl['script_path'] )
    elif job.options['select_source'] == 'bam_url':
        ms_files = {}
        for gid, group in job.groups:
            my_files[gid] = {}
            group_name = re.sub(r'\s+','_',group['name'])
            for rid,run in group['runs'].iteritems():
                bamfile = unique_filename_in(ex.working_directory)
                with open(bamfile,"w") as out:
                    out.write(urllib2.urlopen(run['url']).read())
                with open(bamfile+".bai","w") as out:
                    out.write(urllib2.urlopen(run['url']+".bai").read())
                s = bamstats.nonblocking( ex, bamfile, via='lsf' )
                ms_files[gid][rid] = {'bam': bamfile, 
                                      'stats': s, 
                                      'libname': name+'_'+str(rid)}
        for gid, group in job.groups:
            for rid,run in group['runs'].iteritems():
                ms_files[gid][rid]['stats'] = ms_files[gid][rid]['stats'].wait()
        pdf = add_pdf_stats( ex, ms_files,
                             dict((k,v['name']) for k,v in job.groups.iteritems()), 
                             gl['script_path'] )
        job.options['compute_densities'] = False
#        raise RuntimeError("bam_url not implemented yet!")
    elif job.options['select_source'] == 'mapseq_key':
        M_ms = MiniLIMS(ms_minilims)
        gl_ms = use_pickle(M_ms, "global variables")
        mapseq_url = gl_ms['hts_url']
        (ms_files, ms_job) = import_mapseq_results( str(job.options['mapseq_key']), M_ms, 
                                                    ex.working_directory, gl_ms["hts_url"] )
        g_rep_assembly = g_rep.assembly( ms_job.assembly_id )
        job.groups = ms_job.groups
        for gid, group in job.groups:
            job.groups[gid]['name'] = re.sub(r'\s+','_',group['name'])
        job.options['merge_strand'] = ms.job.options['merge_strand']
        job.options['compute_densities'] = ms.job.options['compute_densities']
        if 'read_extend' in ms.job.options and ms.job.options['read_extend']>0:
            job.options['b2w_args'] = ["-q",str(options['read_extend'])]
    else:
        raise RuntimeError("Didn't know you could do that: "+job.options['select_source'])
    chipseq_files = workflow_groups( ex, job, ms_files, g_rep_assembly.chromosomes, gl['script_path'] )
allfiles = common.get_files( ex.id, M )
if job.options['select_source'] == 'mapseq_key':
    allfiles['url'].update({mapseq_url+"jobs/"+job.options['mapseq_key']+"/get_results":
                            "Mapseq results"})
if 'gdv_project' in job.options and 'sql' in allfiles:
    allfiles['url'] = {job.options['gdv_project']['public_url']: 'GDV view'}
    download_url = gl['hts_download']
    _ = [gdv.add_gdv_sqlite( gl['gdv']['key'], gl['gdv']['email'],
                             job.options['gdv_project']['project_id'],
                             url=download_url+str(k), 
                             name = re.sub('\.sql','',str(f)),
                             gdv_url=gl['gdv']['url'], datatype="quantitative" ) 
         for k,f in allfiles['sql'].iteritems()]
print json.dumps(allfiles)
r = email.EmailReport( sender=gl['email']['sender'],
                       to=str(job.email),
                       subject="Chipseq job "+str(job.description),
                       smtp_server=gl['email']['smtp'] )
r.appendBody('''
Your chip-seq job is finished,
you can retrieve the results at this url:
'''+gl["hts_url"]+"jobs/"+hts_key+"/get_results")
r.send()
