import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from bbcflib.common import unique_filename_in

def _begin(output,format,new,**kwargs):
    """Initialize the figure and axes."""
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
    if 'log' in kwargs:
        opts += ',log="%s"' %kwargs['log']
    opts += ',main="%s"' %kwargs.get('main','')
    opts += ',xlab="%s"' %kwargs.get('xlab','')
    opts += ',ylab="%s"' %kwargs.get('ylab','')
    pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(1,1,1,1),oma=c(3,3,0,3),las=1"
    if len(kwargs.get('mfrow',[])) == 2:
        pars += ",mfrow=c(%i,%i)" %tuple(kwargs['mfrow'])
    robjects.r('par(%s)' %pars)
    return opts, output

def _end(lopts,last,**kwargs):
    """Add the legend and close the figure."""
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
    if not(isinstance(Y,(list,tuple))):
        Y = [Y]
    robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[0]))
    robjects.r("plot(xdata,ydata,pch=19%s)" %plotopt)
    for n in range(1,len(Y)):
        robjects.r.assign('ydata',numpy2ri.numpy2ri(Y[n]))
        robjects.r("points(xdata,ydata,pch=19,col=%i)" %n)
    _end(",pch=19",last,**kwargs)
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
    """Creates a box-and-whiskers plot of `values` split by `labels`."""
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('values',numpy2ri.numpy2ri(values))
    robjects.r.assign('labels',numpy2ri.numpy2ri(labels))
    robjects.r("boxplot(values ~ labels,lty=1,pch=19,varwidth=T)")
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

