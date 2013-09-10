"""
These functions use `rpy2` to bind *R* plotting functions with data from numpy arrays.
Each function takes the following arguments:

* output: the filename, a random name will be generated if this is None (default None),
* format: the image format, 'pdf' (default) or 'png',
* new: boolean indicating if a new figure must be started (default) or if the plot is added to the current figure,
* last: boolean, if true this is the last plot on this figure, the file will be finalized and closed on return.
* **kwargs: additional parameters for *R*, such as 'xlab', 'ylab', 'main, 'mfrow', 'log', 'legend'.
"""

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from bbcflib.common import unique_filename_in
from itertools import combinations
from numpy import array

def list2r(L):
    """Transform a Python list into a string in R format: [1,2,'C'] -> "c(1,2,'C')" ."""
    if not isinstance(L,(list,tuple)): L = [L] # everything is a vector in R
    LL = ["'%s'"%v if isinstance(v,basestring) else str(v) for v in L] # add quotes if elements are strings
    return "c(%s)" % ','.join(LL)

def _begin(output,format,new,ratio=1.375,**kwargs):
    """Initializes the plot in *R*."""
    if new:
        if output is None:
            output = unique_filename_in()
        if format == 'pdf':
            robjects.r('pdf("%s",paper="a4",height=8*%f,width=8)' %(output,ratio))
        elif format == 'png':
            robjects.r('png("%s",height=800*%f,width=800,type="cairo")' %(output,ratio))
        else:
            raise ValueError("Format not supported: %s" %format)
        pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(3,3,1,1),las=1,pch=20"
        if len(kwargs.get('mfrow',[])) == 2:
            pars += ",mfrow=c(%i,%i)" %tuple(kwargs['mfrow'])
        robjects.r('par(%s)' %pars)
    opts = ''
    if 'log' in kwargs: opts += ',log="%s"' %kwargs['log']
    if 'xlim' in kwargs: opts += ',xlim=c(%f,%f)' %tuple(kwargs['xlim'])
    if 'ylim' in kwargs: opts += ',ylim=c(%f,%f)' %tuple(kwargs['ylim'])
    opts += ',main="%s"' %kwargs.get('main','')
    opts += ',xlab="%s"' %kwargs.get('xlab','')
    opts += ',ylab="%s"' %kwargs.get('ylab','')
    return opts, output

def _end(lopts,last,**kwargs):
    """Adds the legend and closes the figure."""
    if not(last): return
    if 'legend' in kwargs:
        names = kwargs.get('legend')
        if names:
            robjects.r("legend(x='topright', legend=%s, col=1:%i%s)" %(list2r(names),len(names),lopts))
    robjects.r("dev.off()")

###################V#########################################
############################################################
def venn(D,options={},output=None,format='png',new=True,last=True,**kwargs):
    """Creates a Venn diagram of D using the *VennDiagram* R package. Up to 4-way.

    :param D: dict `{group_name: count}`. `group_name` is either one name (e.g. 'A')
        or a combination such as 'A|B' - if the items belong to groups A and B.
        `count` must be an integer. Group 'A' means *everything* that is in group 'A',
        including for instance 'A|B' elements.
    :param opts: VennDiagram options, given as a dict `{'option': value}`.
        Ex. `{'euler.d':'TRUE', 'fill':['red','blue']}`
        See `<http://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf>`_
    """
    plotopt,output = _begin(output=output,format=format,new=new,ratio=1,**kwargs)
    groups = sorted([x for x in D if len(x.split('|'))==1])
    ngroups = len(groups)
    combn = [combinations(groups,k) for k in range(1,ngroups+1)]
    combn = ['|'.join(sorted(y)) for x in combn for y in x]
    combn = tuple([D.get(c,0) for c in combn])
    if   ngroups == 1:
        fun = 'single'
        rargs = "area=%d" % D[groups[0]]
    elif ngroups == 2:
        fun = 'pairwise'
        rargs = "area1=%d, area2=%d, cross.area=%d" %combn
        rargs += ", cat.dist=c(0.05,0.05)" # default distance of labels to contour is bad
    elif ngroups == 3:
        fun = 'triple'
        rargs = "area1=%d, area2=%d, area3=%d, n12=%d, n13=%d, n23=%d, n123=%d" %combn
        A,B,C = groups[:3]
        if D[B] == D[A+'|'+B] + D[B+'|'+C] - D[A+'|'+B+'|'+C]:
            rargs += ", cat.dist=c(0.05,0.05,-0.45)"
        else:
            rargs += ", cat.dist=c(0.05,0.05,0.05)"
    elif ngroups == 4:
        fun = 'quad'
        rargs = """area1=%d, area2=%d, area3=%d, area4=%d,
                n12=%d, n13=%d, n14=%d, n23=%d, n24=%d, n34=%d,
                n123=%d, n124=%d, n134=%d, n234=%d, n1234=%d""" %combn
    else:
        return
    rargs += ", category=%s, cex=3, cat.cex=3, fill=1:%i, margin=0.15"%(list2r(groups),ngroups)
    robjects.r("""
plot.new()
library(VennDiagram)
palette(c('grey','blue','green','orange'))
venn.plot = draw.%s.venn(%s)""" % (fun,rargs))
    _end(", pch=%s" %list2r(groups),last,**kwargs)
    return output

############################################################
############################################################
def scatterplot(X,Y,output=None,format='pdf',new=True,last=True,ratio=1.375,**kwargs):
    """Creates a scatter plot of X vs Y.
     If Y is a list of arrays, a different color will be used for each of them."""
    plotopt,output = _begin(output=output,format=format,new=new,ratio=ratio,**kwargs)
    robjects.r.assign('xdata',numpy2ri.numpy2ri(X))
    if not(isinstance(Y,(list,tuple))): Y = [Y]
    robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[0]))
    robjects.r("plot(xdata,ydata%s)" %plotopt)
    for n in range(1,len(Y)):
        robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[n]))
        robjects.r("points(xdata,ydata,col=%i)" %(n+1))
    _end(",pch=20",last,**kwargs)
    return output

############################################################
############################################################
def lineplot(X,Y,output=None,format='pdf',new=True,last=True,ratio=1.375,**kwargs):
    """
    Creates a line plot of X vs Y.
    If Y is a list of arrays, a different color will be used for each of them.
    """
    plotopt,output = _begin(output=output,format=format,new=new,ratio=ratio,**kwargs)
    robjects.r.assign('xdata',numpy2ri.numpy2ri(X))
    if not(isinstance(Y,(list,tuple))): Y = [Y]
    robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[0]))
    robjects.r("plot(xdata,ydata,t='l'%s)" %plotopt)
    for n in range(1,len(Y)):
        robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[n]))
        robjects.r("lines(xdata,ydata,col=%i)" %(n+1))
    _end(",lty=1",last,**kwargs)
    return output

############################################################
############################################################
def boxplot(values,labels,output=None,format='pdf',new=True,last=True,**kwargs):
    """Creates a box-and-whiskers plot of *values* split by *labels*."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('values',numpy2ri.numpy2ri(values))
    robjects.r.assign('labels',numpy2ri.numpy2ri(labels))
    robjects.r("boxplot(values ~ labels,lty=1,varwidth=T)")
    _end("",last,**kwargs)
    return output

############################################################
############################################################
def heatmap(M,output=None,format='pdf',new=True,last=True,
            rows=None,columns=None,orderRows=True,orderCols=True,
            **kwargs):
    """Creates a heatmap of the matrix *M* using *rows* as row labels and *columns* as column labels.
    If either *orderRows* or *orderCols* is True, will cluster accordingly and display a dendrogram."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    if not(rows is None):
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if not(columns is None):
        robjects.r.assign('labCol',numpy2ri.numpy2ri(columns))
        plotopt += ",labCol=labCol"
    if orderCols and orderRows:
        plotopt += ",dendrogram='both',lhei=c(2,10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(0,2,0,0,3,1,0,4),ncol=2)"
    elif orderCols:
        plotopt += ",Rowv=F,dendrogram='column',lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)"
    elif orderRows:
        plotopt += ",Colv=F,dendrogram='row',lhei=c(10,1,2),lwid=c(1,3),mar=c(2,2),lmat=matrix(c(2,0,3,1,0,4),ncol=2)"
    else:
        plotopt += ",Colv=F,Rowv=F,dendrogram='none',lhei=c(10,1,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(1,2,3,4),ncol=1)"
    if kwargs.get('ymin') is not None:
        robjects.r("ymin=%f" %float(kwargs['ymin']))
    else:
        robjects.r("ymin=floor(min(Mdata,na.rm=T))")
    if kwargs.get('ymax') is not None:
        robjects.r("ymax=%f" %float(kwargs['ymax']))
    else:
        robjects.r("ymax=ceiling(max(Mdata,na.rm=T))")
    robjects.r("""
      library(gplots)
      library(RColorBrewer)
      myBreaks = seq(ymin,ymax,length.out=15)
      myColors = rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(length(myBreaks)-1))
#      myCor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
      par(cex.main=1)
      heatmap.2(as.matrix(Mdata),
                col=myColors, trace="none", breaks=myBreaks, #distfun=myCor,
                na.rm=TRUE, density.info='none'%s)""" %plotopt)
    _end("",last,**kwargs)
    return output

############################################################
############################################################
def pairs(M,X=None,labels=None,highlights=([],[]),
          output=None,format='pdf',new=True,last=True,**kwargs):
    """Pairs plot. Returns the name of the output figure.

    :param M,X: If *X* is a vector of length *m*, *M* must be an *n(n+1)/2 x m* matrix
        and the *(i,j)* plot will show *M[ni+j,] ~ X* as a line plot.
        If *X* is `None` then *M* must be an *n x m* matrix and the *(i,j)* plot will
        show *M[i,] ~ M[j,]* as a density plot.
    :param labels: a vector of *n* strings used to label plots, defaults to *1,...,n*.
    :rtype: str
    """
    plotopt,output = _begin(output=output,format=format,new=new,ratio=1,**kwargs)
    if X is None:
        robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
        robjects.r("n = ncol(Mdata); if (exists('X')) rm(X)")
    else:
        robjects.r.assign('X',numpy2ri.numpy2ri(X))
        robjects.r.assign('Mdata',numpy2ri.numpy2ri(M[0][0]))
        for jrow in M[0][1:]:
            robjects.r.assign('Mdj',numpy2ri.numpy2ri(jrow))
            robjects.r("Mdata = cbind(Mdata,Mdj)")
        for imat in M[1:]:
            for jrow in imat:
                robjects.r.assign('Mdj',numpy2ri.numpy2ri(jrow))
                robjects.r("Mdata = cbind(Mdata,Mdj)")
        robjects.r("""
n = as.integer((-1+sqrt(1+8*ncol(Mdata)))/2)  ### ncol(Mdata) = n*(n+1)/2
rowcol = matrix(0,nrow=n,ncol=n)
rowcol[lower.tri(rowcol,diag=T)] = 1:ncol(Mdata)
rowcol = rbind(1+1:n,rowcol)
""")
    if labels is None:
        robjects.r("labels=as.character(1:n)")
    else:
        robjects.r.assign('labels',numpy2ri.numpy2ri(labels))
    if kwargs.get('col') is not None:
        robjects.r("col = '%s'" %kwargs['col'])
    else:
        robjects.r("col = 'orange'")
    if len(highlights[0]) > 1:
        robjects.r("""
hili = %s
hilabs = %s
pch_list = c(3:6,8,15,19,21,22)
col_list = c("green3","blue","red","cyan","magenta","black","yellow","gray")
""" %(list2r(highlights[0]),list2r(highlights[1])))
    robjects.r("""
library(RColorBrewer)
pline1 = function (y, M, X, col, ...) lines(X,M[,y[y[1]]],col=col,...)
pline2 = function (x, y, M, X, col, ...) lines(X,M[,y[x[1]]],col=col,...)
pcor = function(x, y, M, X, ...) {
    usr = par("usr")
    par(usr=c(0, 1, 0, 1))
    cmax = max(M[,x[y[1]]])
    ctext = paste("corr=",format(cmax, digits=4),sep='')
    cmax = abs(cmax)
    if (cmax>0.5) text(0.5, 0.5, ctext, cex=par('cex')*3*cmax)
    par(usr=usr)
}
ppoints = function (x, y, col, ...) {
    if (par("xlog")) x1 = log(x)
    else x1 = x
    if (par("ylog")) y1 = log(y)
    else y1 = y
    colramp = colorRampPalette(c("lightgrey",col),interpolate="spline")
    points(x,y,col=densCols(x1,y1,colramp=colramp),...)
    if (length(hili)>0) {
        for (n in 1:(length(hili)-1)) {
            I = (hili[n]+1):hili[n+1]
            pchn = pch_list[1+(n-1) %% length(pch_list)]
            coln = col_list[1+(n-1) %% length(col_list)]
            points(x[I],y[I],col=coln,pch=pchn,...)
            text(x[I],y[I],hilabs[I-hili[1]],col='black',pos=3,cex=.6)
        }
    }
    abline(0,1,col='black',lty=2)
}
qpoints = function (x, y, col, ...) {
    qq = qqplot(x,y,plot.it=FALSE)
    points(qq$x,qq$y,...)
}
phist = function (x, col, ...) {
    usr = par("usr")
    ylog = par("ylog")
    xlog = par("xlog")
    if (ylog) {
        par(xlog=FALSE,ylog=FALSE)
        h = hist(log(x), plot=FALSE, br=max(10,length(x)/50))
        hbr = (h$breaks-h$breaks[1])/(h$breaks[length(h$breaks)]-h$breaks[1])
        hbr = usr[1]+(usr[2]-usr[1])*hbr
        rect(hbr[-length(hbr)], usr[3], hbr[-1], usr[3]+(usr[4]-usr[3])*h$density,
             col=col, border=NA, ...)
    } else {
        h = hist(x, plot=FALSE, br=max(10,length(x)/50))
        rect(h$breaks[-length(h$breaks)], usr[3],
             h$breaks[-1], usr[3]+(usr[4]-usr[3])*h$density,
             col=col, border=NA, ...)
    }
    par(usr=usr,ylog=ylog,xlog=xlog)
}

if (exists("X")) {
    pairs(rowcol, labels, M=Mdata, X=X, xlim=range(X), ylim=c(-1,1.5), col=col,
          diag.panel=pline1, lower.panel=pcor, upper.panel=pline2)
} else {
    pairs(Mdata, labels, log='xy', col=col,
          diag.panel=phist, lower.panel=qpoints, upper.panel=ppoints)
}
""")
    _end("",last,**kwargs)
    return output

############################################################
############################################################
def hist(X,options={},output=None,format='pdf',new=True,last=True,**kwargs):
    """Create a histogram of the values in vector *X*."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    rargs = ""
    for opt,val in options.iteritems():
        rargs += ", %s=%s" % (opt,list2r(val))
    robjects.r.assign('X',numpy2ri.numpy2ri(X))
    robjects.r("hist(X %s)" % rargs)
    _end("",last,**kwargs)
    return output


############################################################
############################################################
def genomeGraph(chrlist,SP=[],SM=[],F=[],options={},output=None,format='pdf',new=True,last=True,**kwargs):
    """Create a whole genome overview of signals in *SP* (positive), *SM* (negative), and of features in *F*.
    *SP,SM* must be streams with fields 'chr','start','end','score' only,
    *F* must be streams with fields 'chr','start','end','name' only,
    *chrlist* must be a (ordered) list of pairs [(chrname, length),...]."""
    if not(isinstance(SP,(list,tuple))): SP = [SP]
    if not(isinstance(SM,(list,tuple))): SM = [SM]
    if not(isinstance(F,(list,tuple))): F = [F]
    total_tracks = len(SP)+len(SM)+len(F)
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
## cut largest chromosome in 500 segments:
    binsize = 2000*(1+max([c[1] for c in chrlist])/1000000)
    yscale = dict((c,0) for c in [c[0] for c in chrlist])
    robjects.r("spmbins=list()")
    robjects.r("fbins=list()")
#### + signals
    n = 0
    for n,_sp in enumerate(SP):
        _sd = dict((c,[0.]*(l/binsize+1)) for c,l in chrlist)
        for chrom, start, end, score in _sp:
            yscale[chrom] = max(yscale[chrom],score)
            for pos in range(start/binsize,end/binsize+1):
                _sd[chrom][pos] = max(score,_sd[chrom][pos])
        robjects.r("spmbins[[%i]]=list(%s)" %(n+1,",".join("'"+c[0]+"'=list()" \
                                                               for c in chrlist)))
        for chrom in _sd.keys():
            robjects.r.assign('rowtemp',numpy2ri.numpy2ri(array(_sd[chrom])))
            robjects.r("spmbins[[%i]][['%s']]=rowtemp" %(n+1,chrom))
    n0 = n+2
#### - signals
    for n,_sm in enumerate(SM):
        _sd = dict((c,[0.]*(l/binsize+1)) for c,l in chrlist)
        for chrom, start, end, score in _sm:
            yscale[chrom] = max(yscale[chrom],score)
            for pos in range(start/binsize,end/binsize+1):
                _sd[chrom][pos] = min(-score,_sd[chrom][pos])
        robjects.r("spmbins[[%i]]=list(%s)" %(n0+n,",".join("'"+c[0]+"'=list()" \
                                                               for c in chrlist)))
        for chrom in _sd.keys():
            robjects.r.assign('rowtemp',numpy2ri.numpy2ri(array(_sd[chrom])))
            robjects.r("spmbins[[%i]][['%s']]=rowtemp" %(n0+n,chrom))
#### features
    for n,_f in enumerate(F):
        _fd = dict((c[0],[]) for c in chrlist)
        for chrom, start, end, name in _f:
            _fd[chrom].append((start/binsize,end/binsize+1,name))
        robjects.r("fbins[[%i]]=list(%s)" %(n0+n,",".join("'"+c[0]+"'=list()" \
                                                               for c in chrlist)))
        for chrom in _fd.keys():
            robjects.r.assign('rowtemp1',numpy2ri.numpy2ri(array([x[0] for x in _fd[chrom]])))
            robjects.r.assign('rowtemp2',numpy2ri.numpy2ri(array([x[1] for x in _fd[chrom]])))
            robjects.r.assign('rowtemp3',numpy2ri.numpy2ri(array([x[2] for x in _fd[chrom]])))
            robjects.r("fbins[[%i]][['%s']]=list(rowtemp1,rowtemp2,rowtemp3)" %(n+1,chrom))
#### chrlist
    robjects.r("chrlist=list("+",".join(['"%s"=%i'%(c,l) for c,l in chrlist])+")")
    robjects.r.assign('yscale_chrom',numpy2ri.numpy2ri(array(yscale.values())))
    robjects.r.assign('binsize',numpy2ri.numpy2ri(binsize))
    robjects.r("""
n = length(chrlist)
xscale = max(as.numeric(chrlist))
yscale = median(yscale_chrom)
inMb = (xscale>1e7)
if (inMb) {
    unit = "Mb"
    tstep = ceiling(xscale*1e-8)*10
    ticks = seq(tstep*1e6,xscale,by=tstep*1e6)/binsize
    labs = as.character(seq(tstep,xscale*1e-6,by=tstep))
} else {
    unit = "kb"
    tstep = ceiling(xscale*1e-5)*10
    ticks = seq(tstep*1e3,xscale,by=tstep*1e3)/binsize
    labs = as.character(seq(tstep,xscale*1e-3,by=tstep))
}
par(cex.lab=1.5,las=2)
plot(0,0,t='n',xlim=c(0,500),ylim=c(0,n+1),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
abline(v=ticks,lty=2,col="grey")
for (chrom in names(chrlist)) {
    segments(0,n,ceiling(chrlist[[chrom]]/binsize),n,lwd=3)
    mtext(chrom,side=2,at=n)
    colnb = 2
    for (data in spmbins) {
        if (length(data[[chrom]])>0) {
            ynums = n+sapply(sapply(data[[chrom]]/(3*yscale),max,-1),min,1)
            lines(1:length(ynums),ynums,col=colnb)
        }
        colnb = colnb+1
    }
    for (data in fbins) {
        segs = data[[chrom]]
        if (length(segs) == 3 && length(segs[[1]])*length(segs[[2]]) > 0) {
            segments(segs[[1]],n,segs[[2]],n,col=colnb,lwd=5)
            text((segs[[1]]+segs[[2]])/2,n-.2,labels=segs[[3]],cex=.5)
        }
        colnb = colnb+1
    }
    n=n-1
}
axis(side=3,at=ticks,labels=labs,cex.axis=.8,line=-2,las=1,lty=0)
mtext(unit,side=4,at=length(chrlist)+1.6,line=-1.5,las=1)
colnb = colnb-1
""")
    if 'legend' in kwargs:
        names = kwargs.pop('legend')
        robjects.r("legend(x='bottomright',legend=%s,col=2:colnb,pch=15,bg='white')" %list2r(names))
    _end("",last,**kwargs)
    return output
