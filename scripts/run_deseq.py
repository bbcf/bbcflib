#!/archive/epfl/bbcf/bin/bin.x86_64/python
import getopt
import os
import sys
import pickle
import numpy
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc

usage = """run_deseq.py cond1_file cond2_file transcript_names_file cond1_label cond2_label method output

"""

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    narg = 7
    if argv is None:
        argv = sys.argv[1:]
    try:
        if len(argv) != narg:
            raise Usage("run_deseq.py takes exactly %d argument." % narg)
        cond1_file = argv[0]
        with open(cond1_file,'rb') as f:
            cond1 = pickle.load(f)
        cond2_file = argv[1]
        with open(cond2_file,'rb') as f:
            cond2 = pickle.load(f)
        transcript_names_file = argv[2]
        with open(transcript_names_file,'rb') as f:
            transcript_names = pickle.load(f)
        cond1_label = str(argv[3])
        cond2_label = str(argv[4])
        method = str(argv[5])
        output = argv[6]

        data_frame_contents = rlc.OrdDict([(cond1_label+'-'+str(i), robjects.IntVector(c))
                                           for i,c in enumerate(cond1)] +
                                          [(cond2_label+'-'+str(i), robjects.IntVector(c))
                                           for i,c in enumerate(cond2)])
        data_frame = robjects.DataFrame(data_frame_contents)
        data_frame.rownames = transcript_names
        
        conds = robjects.StrVector([cond1_label for x in cond1] + [cond2_label for x in cond2]).factor()
        
        # Import the library
        deseq = rpackages.importr('DESeq')
        cds = deseq.newCountDataSet(data_frame, conds)
        cds = deseq.estimateSizeFactors(cds)
        cds = deseq.estimateVarianceFunctions(cds,method=method)
        res = deseq.nbinomTest(cds, cond1_label, cond2_label)
        res.to_csvfile(output)
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())

