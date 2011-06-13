from numpy import *
from scipy import stats
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from interactive_plot import *
import rpy2.robjects as ro
import rpy2.robjects.packages as rp
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc

def MAplot(data, mode="interactive"):
    """
    data: rpy DataFrame object; each line is of the form
          [gene name, expression under condition 2, epression under condition 2]
    bins: number of intervals for quantile splines
    mode: if "interactive", click on points to display its name
          if "normal", name of genes over the 99% quantile re displayed
    """
    ###
    ### TO IMPLEMENT: GROUPS OF GENES
    ###
    names = list(data[0])
    ratios = array(data[1])/array(data[2])
    means = sqrt(array(data[1])*array(data[2]))
    points = zip(names,ratios,means)
    bins = len(ratios)/10
    intervals = linspace(min(means),max(means),bins)
    annotes = []
    spline_annotes = []

    fig = plt.figure(figsize=[14,10])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)

    #Points
    ax.plot(log10(means), log2(ratios), ".", color="black")

    #Lines (best fit of percentiles)
    for k in [1,5,25,50,75,95,99]:
        h = ones(bins)*0.0001
        for b in range(bins):
            points_in_b = [p for p in points if p[2]>=intervals[b] and p[2]<intervals[b]+1./bins]
            perc = [p[1] for p in points_in_b]
            if points_in_b != []:
                h[b] = stats.scoreatpercentile(perc, k)
            else: h[b] = h[b-1]
            if k==1:
                for p in points_in_b:
                    if p[1]<h[b]: annotes.append(p)
            if k==99:
                for p in points_in_b:
                    if p[1]>h[b]: annotes.append(p)
        x = log10(intervals)
        y = log2(h)
        #coeffs = polyfit(x,y,4)
        #fit = polyval(coeffs,x)
        ax.plot(x,y, color="green", linestyle="--", alpha=0.3)
        spline = UnivariateSpline(x, y, k=3)
        xs = linspace(min(x), max(x), 10*len(x)) #to increase spline smoothness
        ys = spline(xs)
        #ax.plot(x, fit, color="blue")
        #spline_annotes.append((k,x[0],fit[0]))
        ax.plot(xs, ys, color="blue")
        spline_annotes.append((k,xs[0],ys[0]))

    #Decoration
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    if mode == "interactive":
        af = AnnoteFinder( log10(means), log2(ratios), names )
        plt.connect('button_press_event', af)
    if mode == "normal":
        for p in annotes:
            ax.annotate( p[0], xy=(log10(p[2]),log2(p[1])) )
    for l in spline_annotes:
        ax.annotate( str(l[0])+"%", xy=(l[1],l[2]), xytext=(-27,-5), textcoords='offset points' )

    #Display
    plt.draw()
    plt.show()


def rs(len):
    import string
    import random
    return "".join([random.choice(string.letters+string.digits) for x in range(len)])
    

## N=500
## data = random.rand(N,3)
## names = array([rs(6) for k in range(len(data))])
## od = rlc.OrdDict([('names', names),
##                   ('expr1', data[:,1]),
##                   ('expr2', data[:,2]) ])
## D = ro.DataFrame(od)

## MAplot(D, mode="interactive")




## #*************
## # MA-plot d14 vs. d0
## #*************
## nbin=500
## nlag=1
## ratio=log2(refseq.data$ratio14_0)
## ampl=log2(refseq.data$rpkm.d14*refseq.data$rpkm.d0)/2
## retro.ratio=log2(retro.data$ratio14_0)
## retro.ampl=log2(retro.data$rpkm.d14*retro.data$rpkm.d0)/2
## repbase.ratio=log2(repbase.data$ratio14_0)
## repbase.ampl=log2(repbase.data$rpkm.d14*repbase.data$rpkm.d0)/2
## xx=quantile(ampl,0:nbin/nbin)
## probs=matrix(0,ncol=8,nrow=nbin)
## qtls=c(.01,.05,.25,.5,.75,.95,.99)
## for (n in seq(along=xx[1:(nbin+1-nlag)])) {
##    I=which(ampl >= xx[n] & ampl <= xx[n+nlag])
##    probs[n,]=c(mean(xx[n+0:nlag]),quantile(ratio[I],qtls))
## }
## spl=list(low=smooth.spline(probs[,1],probs[,2],df=3)$fit,
##  high=smooth.spline(probs[,1],probs[,8],df=3)$fit)

## Iselect=c(which(repbase.ratio > predict(spl$high,repbase.ampl)$y),
##  which(repbase.ratio < predict(spl$low,repbase.ampl)$y))
## rnames=repbase.data
## Iselect_retro_d14_d0=c(which(retro.ratio > predict(spl$high,retro.ampl)$y),
##  which(retro.ratio < predict(spl$low,retro.ampl)$y))
## Iselect_refseq=c(which(ratio > predict(spl$high,ampl)$y), which(ratio < predict(spl$low,ampl)$y))

## #write selected refseq
## toWrite <- cbind(round(refseq.data[Iselect_refseq,2:7],3),round(log2(refseq.data[Iselect_refseq,5:7]),3))
## rownames(toWrite) <- refseq.data[Iselect_refseq,1]
## colnames(toWrite)[7:9] <- c("log2FC_d2_d0","log2FC_d14_d0","log2FC_d14_d2")
## write.table(toWrite,file="d14_overepresented_refseq.txt",quote=FALSE,sep="\t")
## refseq.de.d14_d0 <- toWrite


## png("MAplot_d14vsd0.png",width=1200,height=1000)
## plot(ampl,ratio,pch='.',main="d14 vs. d0")
## points(repbase.ampl,repbase.ratio,col='red',pch='.',cex=3)
## points(retro.ampl,retro.ratio,col='red',pch=17)

## for (n in 2:ncol(probs)) {
##    spl=smooth.spline(probs[,1],probs[,n],df=3)
##    I=which(spl$x>-8)
##    lines(spl$x[I],spl$y[I],col='blue',lwd=1)
##    text(spl$x[I[1]],spl$y[I[1]],paste(qtls[n-1]*100,"%",sep=''),col='blue',pos=2)
## }
## labels=gsub("_M*","",unlist(strsplit(x=as.character(rnames[Iselect]),"|",fixed=T))[3*(1:length(Iselect))-1])
## xlab=repbase.ampl[Iselect]
## ylab=repbase.ratio[Iselect]
## pos=4
## for (i in order(labels)) {
##    text(x=xlab[i],y=ylab[i],lab=labels[i],col='red',pos=pos,offset=.2,cex=1)
##    pos=pos+1
##    if (pos>4) pos=1
## }
## dev.off()

## #*********************
## ntot=length(ratio)
## prob=sapply(retro.ratio,FUN=function(x) length(which(ratio>x))/ntot)

## png("d14_vs_d0_expression.png",width=800,height=500)
## par(lwd=1.5,cex=1.5)
## hist(ratio,breaks=100,main="d14 overexpression",col='blue',xlab="log2(d14/d0)")
## ypos=seq(3000,100,length.out=length(retro.ratio))
## O=order(retro.ratio)
## segments(x0=retro.ratio[O],y0=0,x1=retro.ratio[O],y1=ypos,col='red',lty=2)
## segments(x0=retro.ratio[O],y0=ypos,x1=5.5,y1=ypos,col='red',lty=2)
## text(x=5.5,y=ypos,lab=paste(retro.data$name,", p=",round(prob,3),sep='')[O],col='red',cex=.6,pos=4,offset=0)
## dev.off()
