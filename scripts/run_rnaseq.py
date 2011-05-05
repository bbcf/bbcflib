#!/bin/env python
"""
deseq_workflow.py

[FIXME: Description]
"""
import os
import sys
import getopt
from bbcflib.rnaseq import *

usage = """workflow.py [-h] [-u via] [-w wdir] minilims job_key

-h           Print this message and exit
-u via       Run executions using method 'via' (can be "local" or "lsf")
-w wdir      Create execution working directories in wdir
minilims     MiniLIMS where RNASeq executions and files will be stored.
job_key      Alphanumeric key specifying the job
"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv=None):
    via = "local"
    minilims_path = None
    job_key = None

    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hu:w:d:k:", ["help","via",
                                                                   "working-directory",
                                                                   "minilims","key"])
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-u", "--via"):
                if a=="local":
                    via = "local"
                elif a=="lsf":
                    via = "lsf"
                else:
                    raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o in ("-w", "--working-directory"):
                if os.path.exists(a):
                    os.chdir(a)
                else:
                    raise Usage("Working directory '%s' does not exist." % a)
            elif o in ("-d", "--minilims"):
                minilims_path = a
            elif o in ("-k", "--key"):
                job_key = a
            else:
                raise Usage("Unhandled option: " + o)

        if len(args) != 0:
            raise Usage("workflow.py takes no arguments without specifiers.")

        if job_key == None:
            raise Usage("Must specify a job key with -k")
        if minilims_path == None:
            raise Usage("Must specify a MiniLIMS to attach to")

        frontend = Frontend('http://htsstation.vital-it.ch/rnaseq/')
        try:
            job = frontend.job(job_key)
        except TypeError, t:
            raise Usage("No such job with key %s at frontend %s" % \
                            (job_key, 'http://htsstation.vital-it.ch/rnaseq/'))

        json = rnaseq_workflow(job, minilims_path, via=via)
        print >>sys.stdout, json

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
    

if __name__ == '__main__':
    sys.exit(main())
