from bbcflib import daflims, genrep, frontend, email, gdv, common
from bbcflib.mapseq import *
import sys, getopt, os

opts = dict(getopt.getopt(sys.argv[1:],"k:d:w:",[])[0])
hts_key = opts['-k']
limspath = opts['-d']
working_dir = opts['-w']

os.chdir(working_dir)
M = MiniLIMS( limspath )
gl = use_pickle(M, "global variables")
htss = frontend.Frontend( url=gl["hts_url"] )
job = htss.job( hts_key )
g_rep = genrep.GenRep( gl["genrep_url"], gl["bwt_root"] )
assembly = g_rep.assembly( job.assembly_id )
dafl = dict((loc,daflims.DAFLIMS( username=gl['lims']['user'],
                                  password=gl['lims']['passwd'][loc] ))
            for loc in gl['lims']['passwd'].keys())
job.options['ucsc_bigwig'] = True
#_ = [M.delete_execution(x) for x in M.search_executions(with_text=hts_key)]
with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
    files = map_groups( ex, job, dafl, ex.working_directory, assembly )
    pdf = add_pdf_stats( ex, files,
                         dict((k,v['name']) for k,v in job.groups.iteritems()),
                         gl['script_path'] )
    if job.options['compute_densities']:
        if not(job.options.get('read_extend')>0):
            job.options['read_extend'] = files.values()[0].values()[0]['stats']['read_length']
        files = densities_groups( ex, job, files, assembly )
        gdv_project = gdv.create_gdv_project( gl['gdv']['key'], gl['gdv']['email'],
                                              job.description, hts_key, 
                                              g_rep_assembly.nr_assembly_id,
                                              gdv_url=gl['gdv']['url'], public=True )
        add_pickle( ex, gdv_project, description='py:gdv_json' )
allfiles = common.get_files( ex.id, M )
if job.options['compute_densities']:
    allfiles['url'] = {gdv_project['public_url']: 'GDV view'}
    download_url = gl['hts_download']
    _ = [gdv.add_gdv_sqlite( gl['gdv']['key'], gl['gdv']['email'],
                             gdv_project['project_id'],
                             url=download_url+str(k), 
                             name = re.sub('\.sql','',str(f)),
                             gdv_url=gl['gdv']['url'], datatype="quantitative" ) 
         for k,f in allfiles['sql'].iteritems()]
print json.dumps(allfiles)
r = email.EmailReport( sender=gl['email']['sender'],
                       to=str(job.email),
                       subject="Mapseq job "+str(job.description),
                       smtp_server=gl['email']['smtp'] )
r.appendBody('''
Your mapseq job is finished,
you can retrieve the results at this url:
'''+gl["hts_url"]+"jobs/"+hts_key+"/get_results")
r.send()
