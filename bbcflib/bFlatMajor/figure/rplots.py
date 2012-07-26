"""
These functions use `rpy2` to bind *R* plotting functions with data from numpy arrays. 
Each function takes the following arguments:

* output: the filename, a random name will be generated if this is None (default None),
* format: the image format, 'pdf' (default) or 'png',
* new: boolean indicating if a new figure must be started (default) or if the plot is added to the current figure,
* last: boolean, if true this is the last plot on this figure, the file will be finalized and close on return.
* **kwargs: additional parameters for *R*, such as 'xlab', 'ylab', 'main, 'mfrow', 'log', 'legend'.
"""

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from bbcflib.common import unique_filename_in

def _begin(output,format,new,**kwargs):
    """Initializes the plot in *R*."""
    if new:
        if output is None:
            output = unique_filename_in()
        if format == 'pdf':
            robjects.r('pdf("%s",paper="a4",height=11,width=8)' %output)
        elif format == 'png':
            robjects.r('png("%s",height=800,width=800)' %output)
        else:
            raise ValueError("Format not supported: %s" %format)
    opts = ''
    if 'log' in kwargs: opts += ',log="%s"' %kwargs['log']
    opts += ',main="%s"' %kwargs.get('main','')
    opts += ',xlab="%s"' %kwargs.get('xlab','')
    opts += ',ylab="%s"' %kwargs.get('ylab','')
    pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(1,1,1,1),oma=c(3,3,0,3),las=1,pch=20"
    if len(kwargs.get('mfrow',[])) == 2:
        pars += ",mfrow=c(%i,%i)" %tuple(kwargs['mfrow'])
    robjects.r('par(%s)' %pars)
    return opts, output

def _end(lopts,last,**kwargs):
    """Adds the legend and closes the figure."""
    if not(last): return
    if 'legend' in kwargs:
        names = kwargs['legend']
        robjects.r("legend('topright', c(%s), col=1:%i%s)" %(",".join(names),len(names),lopts))
    robjects.r("dev.off()")

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
        robjects.r("points(xdata,ydata,col=%i)" %n)
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
        robjects.r("lines(xdata,ydata,col=%i)" %n)
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
    """Creates a heatmap of the matrix `M` using `rows` as row labels and `columns` as column labels.
    If either `orderRows` or `orderCols` is True, will cluster accordingly and display a dendrogram."""
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
      rcor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
      heatmap.2(as.matrix(Mdata),
                col=myColors, trace="none", breaks=myBreaks, distfun=rcor,
                na.rm=TRUE, density.info='none'%s)""" %plotopt)
    _end("",last,**kwargs)
    return output
############################################################
############################################################
def pairs(M,X=None,labels=None,
          output=None,format='pdf',new=True,last=True,**kwargs):
    """Pairs plot. 
    If *X* is a vector of length *m*, *M* must be an *n(n+1)/2 x m* matrix and the *(i,j)* plot will show *M[ni+j,] ~ X* as a line plot. 
    If *X* is `None` then *M* must be an *n x m* matrix and the *(i,j)* plot will show *M[i,] ~ M[j,]* as a density plot.
    *labels* is a vector of *n* strings used to label plots, defaults to *1,...,n*.
    """
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    robjects.r("cdim = ncol(Mdata)")
    if X is None: 
        robjects.r("n=ncol(Mdata)")
    else:
        robjects.r.assign('X',numpy2ri.numpy2ri(X))
        robjects.r("""
n = as.integer((-1+sqrt(1+8*ncol(Mdata)))/2)  ### ncol(Mdata) = n*(n+1)/2
rowcol = matrix(0,nrow=n,ncol=n)
rowcol[lower.tri(rowcol,diag=T)]=1:ncol(Mdata)
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
    text(0.5, 0.5, paste("corr=",format(cmax, digits=4),sep=''),cex=par('cex')*3*cmax)
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
phist = function (x, col, ...) {
    h = hist(x, plot=FALSE, br=max(10,length(x)/50))
    usr = par("usr")
    ylog = par("ylog")
    par(ylog=FALSE)
    par(usr=c(usr[1:2], 0, 1.5*max(h$counts)))
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col=col, border=NA)
    par(usr=usr,ylog=ylog)
}

col = 'red'
if (exists("X")) {
    pairs(rowcol, labels, M=Mdata, X=X, xlim=range(X), col=col,
          diag.panel=pline1, lower.panel=pcor, upper.panel=pline2)
} else {
    pairs(Mdata, labels, log='xy', col=col,
          diag.panel=phist, lower.panel=qpoints, upper.panel=ppoints)
}
""")
    _end("",last,**kwargs)
    return output

