from bbcflib.bFlatMajor.figure import rplots
import random, numpy
from bbcflib import btrack as track
from bbcflib.bFlatMajor.numeric import signal

tfwd = track.track("libv2/temp/rev.sql")
trev = track.track("libv2/temp/fwd.sql")
ccf=signal.correlation([tfwd.read(selection="chr7"),trev.read(selection="chr7")],
                       52960031,56607681)
pdffile=rplots.lineplot(range(-1000,1001),ccf,
                        main='Fwd/rev cross-correlation chr7',
                        ylab='Lag',xlab='Correlation')
#x=numpy.real([random.lognormvariate(0,1) for n in range(1000)])
#y=[numpy.real([random.lognormvariate(0,1) for n in range(1000)]) for k in range(5)]
#pdffile=rplots.scatterplot(x,y,main='Random lognorm',ylab='Y',xlab='X')

import numpy, random

from bbcflib.bFlatMajor.figure import rplots
M=numpy.array([[random.lognormvariate(0,1) for n in range(10)] for m in range(5)])
pdffile=rplots.scatterplot(M[1],M[2],output="scatterpy.pdf")
pdffile=rplots.heatmap(M,rows=["A","B","C","D","E"],orderRows=False,output="hm2py.pdf")
#pdffile=rplots.heatmap(M,output="hm2py.pdf")


