"""
=========================
Module: bbcflib.workflows
=========================

This module provides a base class that can be used to construct workflows from the various
modules in this library.


"""

import sys, os, re, json
from bbcflib import genrep, frontend
from bbcflib.common import normalize_url, get_files, set_file_descr
from bein.util import use_pickle, add_pickle
from bein import execution, MiniLIMS

###DEFAULTS VALUES
_usage = "run_htsstation.py module [-h] [-w wdir] [-k job_key] [-u via] [-c config_file] "
_description = "High-throughput sequencing data analysis workflows."
_basepath = "/archive/epfl/bbcf/data/htsstation"
_via = "lsf"

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

class Workflow(object):
    def __init__(self,**kw):
        self.module = kw.get('module','HTSstation')
#### name can be different from module: e.g. c4seq is the python module,
#### but the name 4cseq is used in outputs and in variable names
        self.name = kw.get('name',self.module)
        self.opts = (
            ("-v", "--via", "Run executions locally or remotely (can be 'local' or 'lsf')",
             {'default': _via}),
            ("-k", "--key", "Alphanumeric key of the %s job" %self.name,
             {'default': None}),
            ("-w", "--working-directory", "Directory to run execution in",
             {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Configuration file", {'default': None}),
            ("--basepath","HTS data basepath (%s MiniLIMS)"%self.name, {'default': _basepath}))
        self.opts += kw.get("opts",())
        self.usage = _usage + str(kw.get('usage',''))
        self.desc = kw.get('desc',_description)
#### By default the workflow will execute the call
####     X_workflow(ex,**self.main_args) from bbcflib.X where X is the module name
#### Can be overloaded in derived classes
        __import__('bbcflib.'+self.module)
        self.sysmod = getattr(sys.modules['bbcflib'],self.module)
        self.main_func = getattr(sys.modules["bbcflib."+self.module],self.module+"_workflow")
        self.main_args = {}
        self.logfile = None
        self.debugfile = None


    def __call__(self,opts):
        self.opts = opts
        if os.path.exists(self.opts.wdir):
            os.chdir(self.opts.wdir)
        else:
            raise Usage("Working directory '%s' does not exist." %self.opts.wdir)

##### Connect to Minilims, recover global variables, fetch job info
        self.minilims = os.path.join(self.opts.basepath,self.name+"_minilims")
        M = MiniLIMS(self.minilims)
        if not((self.opts.key != None or (self.opts.config and os.path.exists(self.opts.config)))):
            raise Usage("Need a job key or a configuration file")
        if self.opts.key:
            self.globals = use_pickle(M, "global variables")
            htss = frontend.Frontend( url=self.globals['hts_mapseq']['url'] )
            self.job = htss.job( self.opts.key )
            [M.delete_execution(x) for x in \
                 M.search_executions(with_description=self.opts.key,fails=True)]
            if self.job.options.get("config_file"):
                if os.path.exists(self.job.options["config_file"]):
                    self.opts.config = os.path.abspath(self.job.options["config_file"])
                elif os.path.exists("config.txt"):
                    self.opts.config = os.path.abspath("config.txt")
            if self.opts.config and os.path.exists(self.opts.config):
                (self.job,self.globals) = frontend.parseConfig( self.opts.config, self.job, self.globals )
        elif os.path.exists(self.opts.config):
            (self.job,self.globals) = frontend.parseConfig( self.opts.config )
            self.opts.key = self.job.description
        else:
            raise Usage("Need either a job key (-k) or a configuration file (-c).")
##### Genrep instance
        if 'fasta_file' in self.job.options:
            if os.path.exists(self.job.options['fasta_file']):
                self.job.options['fasta_file'] = os.path.abspath(self.job.options['fasta_path'])
            else:
                for ext in (".fa",".fa.gz",".tar.gz"):
                    if os.path.exists("ref_sequence"+ext):
                        self.job.options['fasta_file'] = os.path.abspath("ref_sequence"+ext)
            if not os.path.exists(self.job.options['fasta_file']):
                raise Usage("Don't know where to find fasta file %s." %self.job.options["fasta_file"])
        g_rep = genrep.GenRep( url=self.globals.get("genrep_url"),
                               root=self.globals.get("bwt_root") )
##### Configure facility LIMS
        if 'lims' in self.globals:
            from bbcflib import daflims
            self.job.dafl = dict((loc,daflims.DAFLIMS( username=self.globals['lims']['user'],
                                                       password=pwd ))
                                 for loc,pwd in self.globals['lims']['passwd'].iteritems())
        else:
            self.job.dafl = None
########################################################################
##########################  EXECUTION  #################################
########################################################################
##### Logging
        self.logfile = open(self.opts.key+".log",'w')
        self.debugfile = open(self.opts.key+".debug",'w')
        self.debug_write(json.dumps(self.globals)+"\n")
        with execution( M, description=self.opts.key,
                        remote_working_directory=self.opts.wdir ) as ex:
            self.log_write("Enter execution. Current working directory: %s" %ex.working_directory)
            self.job.assembly = genrep.Assembly( assembly=self.job.assembly_id,
                                                 genrep=g_rep,
                                                 fasta=self.job.options.get('fasta_file'),
                                                 annot=self.job.options.get('annot_file'),
                                                 intype=self.job.options.get('input_type_id',0),
                                                 ex=ex, via=self.opts.via )
##### Check all the options
            if not self.check_options():
                raise Usage("Problem with options %s" %self.opts)
            self.debug_write(json.dumps(self.job.options))
            self.init_files( ex )
##### Run workflow
            self.log_write("Starting workflow.")
            self.main_func(ex,**self.main_args)
            if self.job.options['create_gdv_project']: self.gdv_create(ex)

########################################################################
########################  POSTPROCESSING  ##############################
########################################################################
        allfiles = get_files( ex.id, M )
        if self.job.options['gdv_project'].get('project',{}).get('id',0)>0:
            allfiles['url'] = self.gdv_upload(allfiles.get('sql',{}))
        self.logfile.close()
        self.debugfile.close()
        print json.dumps(allfiles)
        with open(self.opts.key+".done",'w') as done: json.dump(allfiles,done)
        self.send_email()
        return 0
########################################################################
#############################  DONE!  ##################################
########################################################################

    def log_write(self,mesg):
        self.logfile.write(mesg+"\n")
        self.logfile.flush()

    def debug_write(self,mesg):
        self.debugfile.write(mesg+"\n")
        self.debugfile.flush()

    def check_options(self,more_defs={}):
    ### { option_name: (default_value, additional_conditions, ...) }
        defaults = {'compute_densities':(True,),
                    'ucsc_bigwig':(True,(self.job.assembly.intype == 0)),
                    'create_gdv_project':(False,)}
        defaults.update(more_defs)
        for op,val in defaults.iteritems():
            self.job.options.setdefault(op,val[0])
            if isinstance(self.job.options[op],basestring):
                self.job.options[op] = self.job.options[op].lower() in ['1','true','t']
            self.job.options[op] &= all(val[1:])
        self.job.options['gdv_project'] = {'project':{'id': self.job.options.get('gdv_project_id',0)}}
        self.job.options.setdefault('gdv_key',"")
        if hasattr(self.opts,"mapseq_minilims") and self.opts.mapseq_minilims:
            self.mapseq_minilims = self.opts.mapseq_minilims
        else:
            self.mapseq_minilims = os.path.join(self.opts.basepath,"mapseq_minilims")
        self.suffix = ['fwd','rev']
        return

    def init_files(self,ex):
        """ Default behaviour for most modules: get bam files from mapping."""
        from bbcflib.mapseq import get_bam_wig_files
        msurl = self.globals.get('hts_mapseq',{}).get('url','')
        scpath = self.globals.get('script_path','')
        self.job = get_bam_wig_files( ex, self.job,
                                      minilims=self.mapseq_minilims,
                                      hts_url=msurl, suffix=self.suffix,
                                      script_path=scpath, via=self.opts.via )
        return True

    def gdv_create(self,ex):
        from bbcflib import gdv
        project = gdv.get_project(mail=self.globals['gdv']['email'],
                                  key=self.globals['gdv']['key'],
                                  project_key=self.job.options['gdv_key'])
        if 'error' in project:
            self.log_write("Creating GDV project.")
            project = gdv.new_project( self.globals['gdv']['email'],
                                       self.globals['gdv']['key'],
                                       self.job.description,
                                       self.job.assembly.id,
                                       self.globals['gdv']['url'] )
            self.debug_write("\nGDV project: "+json.dumps(project))
            add_pickle( ex, project, description=set_file_descr("gdv_json",step='gdv',type='py',view='admin') )
        self.job.options['gdv_project'] = project
        return True


    def gdv_upload(self,files):
        glg = self.globals['gdv']
        project = self.job.options['gdv_project']['project']
        download_url = normalize_url(self.globals['hts_'+self.name]['download'])
        urls_names = dict(("%s/%s" %(download_url,k),re.sub('\.sql.*','',str(f)))
                          for k,f in files.iteritems())
        self.log_write("Uploading GDV tracks:\n"+" ".join(urls_names.keys())+"\n"+" ".join(urls_names.values()))
        try:
            tr = gdv.multiple_tracks(mail=glg['email'], key=glg['key'], serv_url=glg['url'],
                                     project_id=project['id'],
                                     extensions=['sql']*len(urls_names),
                                     urls=urls_names.keys(), names=urls_names.values(), force=True )
            self.debug_write("GDV Tracks Status\n"+"\n".join([str(v) for v in tr]))
        except Exception, err:
            self.debug_write("GDV Tracks Failed: %s" %err)
            pass
        gdv_project_url = "%s/public/project?k=%s&id=%s" %(normalize_url(glg['url']),project['download_key'],project['id'])
        return {gdv_project_url: 'GDV view'}

    def send_email(self):
        if not('email' in self.globals and hasattr(self.job,"email")):
            return
        from bbcflib import email
        r = email.EmailReport( sender=self.globals['email']['sender'],
                               to=str(self.job.email).split(','),
                               subject="%s job %s" %(self.name,str(self.job.description)),
                               smtp_server=self.globals['email']['smtp'] )
        r.appendBody('''
Your %s job has finished.

The description was:
%s
and its unique key is %s

You can now retrieve the results at this url:
%s/jobs/%s/get_results
''' %(self.name,self.job.description,self.opts.key,normalize_url(self.globals['hts_'+self.name]['url']),self.opts.key))
        r.send()




#------------------------------------------------------#
# This code was written by                             #
# Bioinformatics and Biostatistics Core Facility, EPFL #
# http://bbcf.epfl.ch/                                 #
# webmaster.bbcf@epfl.ch                               #
#------------------------------------------------------#
