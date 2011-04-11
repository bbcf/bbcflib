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
with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
    if job.options['select_source'] == 'lims':
        g_rep_assembly = g_rep.assembly( job.assembly_id )
        dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'],
                                          password=gl['lims']['passwd'][loc] ))
                    for loc in gl['lims']['passwd'].keys())
        job.options['discard_pcr_duplicates'] = True
        processed = map_groups( ex, job, dafl, ex.working_directory, g_rep )
        pdf = add_pdf_stats( ex, processed,
                             dict((k,v['name']) for k,v in job.groups.iteritems()), 
                             gl['script_path'] )
    elif job.options['select_source'] == 'bam_url':
#        bamfile = unique_filename_in(gl['fastq_root'])
#        with open(bamfile,"w") as out:
#            out.write(urllib2.urlopen(job.options['bam_url']).read())
#        with open(bamfile+".bai","w") as out:
#            out.write(urllib2.urlopen(job.options['bam_url']+".bai").read())
#    ***run bamstat
#        pdf = add_pdf_stats( ex, {??}, dict(??), gl['script_path'] )
        raise RuntimeError("bam_url not implemented yet!")
    elif job.options['select_source'] == 'mapseq_key':
        M_ms = MiniLIMS(ms_minilims)
        gl_ms = use_pickle(M_ms, "global variables")
        mapseq_url = gl_ms['hts_url']
        (processed, ms_job) = import_mapseq_results( str(job.options['mapseq_key']), M_ms, 
                                                     ex.working_directory, gl_ms["hts_url"] )
        g_rep_assembly = g_rep.assembly( ms_job.assembly_id )
        job.groups = ms_job.groups
    else:
        raise RuntimeError("Didn't know you could do that: "+job.options['select_source'])
    (p,g) = workflow_groups( ex, job, processed, g_rep_assembly.chromosomes, gl['script_path'] )
    gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                          job.description, hts_key, 
                                          g_rep_assembly.nr_assembly_id,
                                          gdv_url=gl['gdv']['url'], public=True )
    add_pickle( ex, gdv_project, description='py:gdv_json' )
allfiles = common.get_files( ex.id, M )
allfiles['url'] = {gdv_project['public_url']: 'GDV view'}
if job.options['select_source'] == 'mapseq_key':
    allfiles['url'].update({mapseq_url+"jobs/"+job.options['mapseq_key']+"/get_results":
                            "Mapseq results"})
print json.dumps(allfiles)
download_url = gl['hts_download']
_ = [gdv.add_gdv_sqlite( gl['gdv']['key'], gl['gdv']['email'],
                         gdv_project['project_id'],
                         url=download_url+str(k), 
                         name = re.sub('\.sql','',str(f)),
                         gdv_url=gl['gdv']['url'], datatype="quantitative" ) 
     for k,f in allfiles['sql'].iteritems()]
r = email.EmailReport( sender=gl['email']['sender'],
                       to=str(job.email),
                       subject="Chipseq job "+str(job.description),
                       smtp_server=gl['email']['smtp'] )
r.appendBody('''
Your chip-seq job is finished,
you can retrieve the results at this url:
'''+gl["hts_url"]+"jobs/"+hts_key+"/get_results")
r.send()
