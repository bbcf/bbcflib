import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from bbcflib.common import unique_filename_in

def _begin(output,format,new,**kwargs):
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
    pars = "lwd=2,cex=1.1,cex.main=1.5,cex.lab=1.3,cex.axis=1.1,mar=c(1,1,1,.5),oma=c(3,3,0,0),las=1"
    if len(kwargs.get('mfrow',[])) == 2:
        pars += ",mfrow=c(%i,%i)" %tuple(kwargs['mfrow'])
    robjects.r('par(%s)' %pars)
    return opts, output

def _end(lopts,last,**kwargs):
    if not(last): return
    if 'legend' in kwargs:
        names = kwargs['legend']
        robjects.r("legend('topright', c(%s), col=1:%i%s)" %(",".join(names),len(names),lopts))
    robjects.r("dev.off()")
    
############################################################
############################################################
def scatterplot(X,Y,output=None,format='pdf',new=True,last=True,**kwargs):
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
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('xdata',numpy2ri.numpy2ri(X))
    if not(isinstance(Y,(list,tuple))):
        Y = [Y]
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
    plotopt,output = _begin(output=output,format=format,new=new,**kwargs)
    robjects.r.assign('Mdata',numpy2ri.numpy2ri(M))
    if rows:
        robjects.r.assign('labRow',numpy2ri.numpy2ri(rows))
        plotopt += ",labRow=labRow"
    if columns:
        robjects.r.assign('labCol',numpy2ri.numpy2ri(columns))
        plotopt += ",labCol=labCol"
    if not(orderRows):
        plotopt += ",Rowv=F"
    if not(orderCols):
        plotopt += ",Colv=F"
    robjects.r("""
      library(gplots)
      library(RColorBrewer)
      myBreaks=seq(floor(min(M,na.rm=T)),ceiling(max(M,na.rm=T)),length.out=15)
      myColors=rev(colorRampPalette(brewer.pal(10,"RdYlBu"))(length(myBreaks)-1))
      rcor = function(x) {as.dist(1-cor(t(x),use="pairwise.complete.ob"))}
      heatmap.2(as.matrix(Mdata), 
                col=myColors, trace="none", breaks=myBreaks, distfun=rcor,
                dendrogram="column", na.rm=TRUE, density.info='none',
                lhei=c(2,10,1,2),lwid=c(1),mar=c(2,2),lmat=matrix(c(3,1,2,4),ncol=1)%s)
    """ %plotopt)
    _end("",last,**kwargs)
    return output

