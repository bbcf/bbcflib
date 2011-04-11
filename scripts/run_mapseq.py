from bbcflib import daflims, genrep, frontend, email, common
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
with execution( M, description=hts_key, remote_working_directory=working_dir ) as ex:
    processed = map_groups( ex, job, dafl, ex.working_directory, assembly )
    pdf = add_pdf_stats( ex, processed,
                         dict((k,v['name']) for k,v in job.groups.iteritems()),
                         gl['script_path'] )
    if job.options['compute_densities']:
        b2w_args = []
        if 'b2w_args' in job.options:
            b2w_args = job.options['b2w_args']
        if 'read_extend' in job.options and job.options['b2w_args']>0:
            b2w_args += ["-q",str(job.options['b2w_args'])]
        processed = densities_groups( ex, job, processed, assembly, b2w_args )
allfiles = common.get_files( ex.id, M )
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
