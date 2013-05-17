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

def list2r(L):
    """Transform a Python list into a string in R format: [1,2,'C'] -> "c(1,2,'C')" ."""
    if  not isinstance(L,(list,tuple)): L = [L] # everything is a vector in R
    LL = ["'%s'"%v if isinstance(v,str) else str(v) for v in L] # add quotes if elements are strings
    return "c(%s)" % ','.join(LL)

def _begin(output,format,new,**kwargs):
    """Initializes the plot in *R*."""
    if new:
        if output is None:
            output = unique_filename_in()
        if format == 'pdf':
            robjects.r('pdf("%s",paper="a4",height=11,width=8)' %output)
        elif format == 'png':
            robjects.r('png("%s",height=800,width=800,type="cairo")' %output)
        else:
            raise ValueError("Format not supported: %s" %format)
        pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(3,3,1,1),oma=c(3,3,0,3),las=1,pch=20"
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
        names = kwargs['legend']
        robjects.r("legend(x='topright', legend=%s, col=1:%i%s)" %(list2r(names),len(names),lopts))
    robjects.r("dev.off()")

###################V#########################################
############################################################
def venn(D,legend=None,options={},output=None,format='png',new=True,last=True,**kwargs):
    """Creates a Venn diagram of D using the *VennDiagram* R package.

    :param D: dict `{group_name: count}`. `group_name` is either one name (e.g. 'A')
        or a combination such as 'A|B' - if the items belong to groups A and B.
        `count` must be an integer. Group 'A' means *everything* that is in group 'A',
        including for instance 'A|B' elements.
    :param legend: (list of str) associate each group name to a legend that will be
        displayed in a corner.
    :param opts: VennDiagram options, given as a dict `{'option': value}`.
        Ex. `{'euler.d':'TRUE', 'fill':['red','blue']}`
        See `<http://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf>`_
    """
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    groups = sorted([x for x in D if len(x.split('|'))==1])
    ngroups = len(groups)
    combn = [combinations(groups,k) for k in range(1,ngroups+1)]
    combn = ['|'.join(sorted(y)) for x in combn for y in x]
    if   ngroups == 1:
        fun = 'draw.single.venn'
        rargs = "area=%d" % D[groups[0]]
    elif ngroups == 2:
        fun = 'draw.pairwise.venn'
        rargs = "area1=%d, area2=%d, cross.area=%d" % tuple([D.get(c,0) for c in combn])
        rargs += ", cat.dist=c(0.05,0.05)" # default distance of labels to contour is bad
    elif ngroups == 3:
        fun = 'draw.triple.venn'
        rargs = """area1=%d, area2=%d, area3=%d,
                n12=%d, n13=%d, n23=%d, n123=%d""" % tuple([D.get(c,0) for c in combn])
        A=groups[0]; B=groups[1]; C=groups[2]
        if D[B] == D[A+'|'+B] + D[B+'|'+C] - D[A+'|'+B+'|'+C]:
            rargs += ", cat.dist=c(0.05,0.05,-0.45)"
        else:
            rargs += ", cat.dist=c(0.05,0.05,0.05)"
    elif ngroups == 4:
        fun = 'draw.quad.venn'
        rargs = """area1=%d, area2=%d,  area3=%d, area4=%d,
                n12=%d, n13=%d, n14=%d, n23=%d, n24=%d, n34=%d,
                n123=%d, n124=%d, n134=%d, n234=%d, n1234=%d""" % tuple([D.get(c,0) for c in combn])
    else: return
    rargs += ", category=%s" % list2r(groups) # group names
    colors = ['grey','blue','green','orange']
    globalopt = {'cex':3, 'cat.cex':3,
                 'fill':colors[:ngroups]}
    for opt,val in globalopt.iteritems():
        rargs += ", %s=%s" % (opt,list2r(options.get(opt,val)))
    if legend: # not in _end() because it requires plot.new()
        legend_opts = "x='topright', pch=%s, legend=%s" % (list2r(groups), list2r(legend))
        legend_opts += ", cex=1.5"
        robjects.r("plot.new()")
        robjects.r("legend(%s)" % legend_opts)
    rargs += ', margin=0.15, pty="s" ' # pty does not work... pdf format is stretched
    robjects.r("library(VennDiagram)")
    robjects.r("venn.plot <- %s(%s%s)" % (fun,rargs,plotopt))
    _end('',last,**kwargs)
    return output

############################################################
############################################################
def scatterplot(X,Y,output=None,format='pdf',new=True,last=True,**kwargs):
    """Creates a scatter plot of X vs Y.
     If Y is a list of arrays, a different color will be used for each of them."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
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
def lineplot(X,Y,output=None,format='pdf',new=True,last=True,**kwargs):
    """Creates a line plot of X vs Y.
     If Y is a list of arrays, a different color will be used for each of them."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
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
    robjects.r("""
      library(gplots)
      library(RColorBrewer)
      myBreaks=seq(floor(min(Mdata,na.rm=T)),ceiling(max(Mdata,na.rm=T)),length.out=15)
      myColors=rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(length(myBreaks)-1))
#      myCor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
      par(cex.main=1)
      heatmap.2(as.matrix(Mdata),
                col=myColors, trace="none", breaks=myBreaks, #distfun=myCor,
                na.rm=TRUE, density.info='none'%s)""" %plotopt)
    _end("",last,**kwargs)
    return output

############################################################
############################################################
def pairs(M,X=None,labels=None,
          output=None,format='pdf',new=True,last=True,**kwargs):
    """Pairs plot. Returns the name of the output figure.

    :param M,X: If *X* is a vector of length *m*, *M* must be an *n(n+1)/2 x m* matrix
        and the *(i,j)* plot will show *M[ni+j,] ~ X* as a line plot.
        If *X* is `None` then *M* must be an *n x m* matrix and the *(i,j)* plot will
        show *M[i,] ~ M[j,]* as a density plot.
    :param labels: a vector of *n* strings used to label plots, defaults to *1,...,n*.
    :param output: (str) name of the output file
    :param format: (str) image format
    :param new: (bool) ?
    :param last: (bool) ?
    :rtype: str
    """
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
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
    colramp = colorRampPalette(c("lightgrey",col),interpolate="spline")
    points(x,y,col=densCols(x,y,colramp=colramp),...)
    abline(0,1,col='black',lty=2)
}
qpoints = function (x, y, col, ...) {
    colramp = colorRampPalette(c("lightgrey",col),interpolate="spline")
    qq = qqplot(x,y,plot.it=FALSE)
    points(qq$x,qq$y,...)
}
ptext = function(x=0.5, y=0.5, txt, cex, font) {
    ylog = par("ylog")
    xlog = par("xlog")
    if (xlog) par(xlog=FALSE,ylog=FALSE)
    text(x, y, txt, cex=cex, font=font)
    par(ylog=ylog,xlog=xlog)
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

col = 'red'
if (exists("X")) {
    pairs(rowcol, labels, M=Mdata, X=X, xlim=range(X), ylim=c(-1,1.5), col=col,
          diag.panel=pline1, lower.panel=pcor, upper.panel=pline2)
} else {
    pairs(Mdata, labels, log='xy', col=col,
          diag.panel=phist, lower.panel=qpoints, upper.panel=ppoints, text.panel=ptext)
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
def genomeGraph(chrmeta,SP=[],SM=[],F=[],options={},output=None,format='pdf',new=True,last=True,**kwargs):
    """Create a whole genome overview of signals in *SP* (positive), *SM* (negative), and of features in *F*."""
    if not(isinstance(SP,(list,tuple))): SP = [SP]
    if not(isinstance(SM,(list,tuple))): SM = [SM]
    if not(isinstance(F,(list,tuple))): F = [F]
    total_tracks = len(SP)+len(SM)+len(F)
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
## cut largest chromosome in 500 segments:
    binsize = 2000*(1+max([c['length'] for c in chrmeta.values()])/1000000)
    yscale = dict((c,0) for c in chrmeta.keys())
    robjects.r("spmbins=list()")
#### + signals
    n = 0
    for n,_sp in enumerate(SP):
        _sd = dict((c,[0.]*(v['length']/binsize+1)) for c,v in chrmeta.iteritems())
        for chrom, start, end, score in _sp:
            yscale[chrom] = max(yscale[chrom],score)
            for pos in range(start/binsize,end/binsize+1):
                _sd[chrom][pos] = max(score,_sd[chrom][pos])
        robjects.r("spmbins[[%i]]=list(%s)" %(n+1,",".join("'"+c+"'=list()" \
                                                               for c in chrmeta.keys())))
        for chrom in _sd.keys():
            robjects.r.assign('rowtemp',numpy2ri.numpy2ri(_sd[chrom]))
            robjects.r("spmbins[[%i]][['%s']]=rowtemp" %(n+1,chrom))
    n0 = n+2
#### - signals
    for n,_sm in enumerate(SM):
        _sd = dict((c,[0.]*(v['length']/binsize+1)) for c,v in chrmeta.iteritems())
        for chrom, start, end, score in _sm:
            yscale[chrom] = max(yscale[chrom],score)
            for pos in range(start/binsize,end/binsize+1):
                _sd[chrom][pos] = min(-score,_sd[chrom][pos])
        robjects.r("spmbins[[%i]]=list(%s)" %(n0+n,",".join("'"+c+"'=list()" \
                                                                for c in chrmeta.keys())))
        for chrom in _sd.keys():
            robjects.r.assign('rowtemp',numpy2ri.numpy2ri(_sd[chrom]))
            robjects.r("spmbins[[%i]][['%s']]=rowtemp" %(n0+n,chrom))
#### features
    for n,_f in enumerate(F):
        _fd = dict((c,[]) for c in chrmeta.keys())
        for chrom, start, end, name in _f:
            _fd[chrom].append(start/binsize,end/binsize+1,name)
        robjects.r("fbins[[%i]]=list(%s)" %(n0+n,",".join("'"+c+"'=list()" \
                                                              for c in chrmeta.keys())))
        for chrom in _fd.keys():
            robjects.r.assign('rowtemp',numpy2ri.numpy2ri(_fd[chrom]))
            robjects.r("fbins[[%i]][['%s']]=rowtemp" %(n+1,chrom))
#### chrlist
    robjects.r("chrlist=list("+",".join(['"%s"=%i'%(k,v['length']) \
                                             for k,v in chrmeta.iteritems()])+")")
    robjects.r.assign('yscale_chrom',numpy2ri.numpy2ri(yscale.values()))
    robjects.r("""
n = length(chrlist)
yscale = median(yscale_chrom)
par(cex.lab=1.5,las=2)
plot(0,0,t='n',xlim=c(0,500),ylim=c(0,n+1),xlab='',ylab='',bty='n',xaxt='n',yaxt='n')
abline(v=seq(2e7,xscale,by=2e7),lty=2,col="grey")
axis(side=3,at=seq(2e7,xscale,by=2e7),
     labels=as.character(seq(2,xscale*1e-7,by=2)),las=1,lty=0)
mtext(paste("x",expression(10^7)),side=4,at=n+1)
for (chrom in names(chrlist)) {
    segments(0,n,chrlist[[chrom]],n,lwd=3)
    mtext(chrom,side=2,at=n)
    colnb = 1
    for (data in spmbins) {
        ynums = n+sapply(sapply(data[[chrom]]/(3*yscale),max,-1),min,1)
        lines(1:length(ynums),ynums,col=colnb)
        colnb = colnb+1
    }
    n=n-1
}
"""2)
    _end("",last,**kwargs)
    return output
