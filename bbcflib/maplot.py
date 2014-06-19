"""
======================
Module: bbcflib.mapseq
======================

Creates an `MA-plot` to compare transcription levels of a set of genes
(or other features) in two different conditions, from a CSV file.
One can enter several datasets (CSV files) in the same format, each of which
will be plotted in a different color, and be annotated if requested.

The class AnnoteFinder is used to create interactive - clickable - plots.
"""

import sys, os, json, math, csv, random, string
import numpy
from scipy import stats
from numpy import asarray,log,log10,log2,exp,sqrt,mean,median,float_,round,nonzero


def unique_filename(len=20, path=None):
    """Return a *len*-characters random filename, unique in the given *path*."""
    if not path: path = os.getcwd()
    def random_string():
        return "".join([random.choice(string.letters + string.digits) for x in range(len)])
    while True:
        filename = random_string()
        files = [f for f in os.listdir(path) if f.startswith(filename)]
        if files == []: break
    return filename

def guess_file_format(f, sep=None):
    """Guess format of a *sep*-delimited file object *f*. Returns a ``csv.Dialect`` object and the
    header of the file."""
    assert isinstance(f,file), "f must be a file object"
    header = 'None'
    f.readline()
    dialect = csv.Sniffer().sniff(f.readline()); f.seek(0)
    if sep: dialect.delimiter = sep
    if csv.Sniffer().has_header(f.readline()):
        f.seek(0); header=f.readline()
    return dialect, header

def _name_or_index(cols, dialect, header):
    """Given an array *cols*, detect if elements are indices of *header* or elements of *header*.
    Returns a dictionary of Python indices of the form {1:[2,3],2:[4,5]} if group 1 is made of
    runs (columns) 2 and 3, and group 2 of runs 4 and 5."""
    if cols[0] is None: return None
    cols = eval(str(cols))
    if isinstance(cols, dict):
        for k,v in cols.iteritems(): cols[k] = [x-1 for x in v]
        return cols
    elif all([c in header.split(dialect.delimiter) for c in cols]):
        cols = [[header.split(dialect.delimiter).index(c)] for c in cols]
    else:
        try: cols = [[int(c)-1] for c in cols]
        except ValueError:
            raise ValueError("Argument must contain column names or indices (got %s). \
                              Detected header: \n%s" % (cols,header))
    cols = dict(zip(range(1,len(cols)+1),cols))
    return cols

def read_data(filename,sep,cols,labels,lower):
    names=[]; counts=[]
    with open(filename,'r') as f:
        # Guess file format & extract columns indexes
        dialect,header = guess_file_format(f,sep)
        try: csvreader = csv.reader(f, dialect=dialect, quoting=csv.QUOTE_NONE)
        except TypeError: csvreader = csv.reader(f, dialect='excel-tab', quoting=csv.QUOTE_NONE)
        pycols = _name_or_index(cols, dialect, header)
        pylabels = _name_or_index(labels, dialect, header)
        pylabels = [x[0] for x in pylabels.values()]
        # Read the file
        for row in csvreader:
            try:
                c1 = mean([float(row[x]) for x in pycols[1]])
                c2 = mean([float(row[x]) for x in pycols[2]])
            except ValueError: continue # Skip line if contains NA, nan, etc.
            except IndexError:
                raise IndexError("list index out of range. Probably wrong separator given with -s (got \"%s\")" % sep)
            if (c1*c2 > lower):
                label = '|'.join([row[x] for x in pylabels])
                counts.append((c1,c2))
                names.append(label)
    counts = asarray(counts)
    return counts, names

def _normalize(counts,mode):
    if mode == 'tags':
        colsum = counts.T.sum(1)
        res = (counts / colsum) * (10**round(log10(max(colsum))))
        sf = None
    elif mode == 'sf':
        loggeomeans = mean(log(counts),1)
        diff = log(counts).T - loggeomeans
        sf = exp(median(diff, 1))
        res = counts / sf
    return res, sf

def MAplot(dataset, cols=['2','3'], labels=['1'], annotate=None, mode="normal", data_format="counts",
           sep=None, limits=[None,None,None,None], slimits=[None,None], deg=3, bins=30,
           assembly_id=None, normalize="sf", quantiles=True, title="MA-plot", extremes=False, prefix=None):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    :param dataset: (list or string) names of up to six text files.
    :param cols: (list or dict) 1-based indices of the two columns containing the numeric data to compare (e.g. `[1,2]`),
        or column names (e.g. `['g1','g2']`), or a dict (e.g. `{1:[`indices of group1`], 2:[`indices of group 2`]}`).
    :param labels: (int or list) 1-based indices of the columns used as labels, or column names.
    :param annotate: (list) in ``normal`` mode, choose which for which datasets you want the
        points to be labeled. Enter 1 to annotate, 0 not to annotate, in the
        same order as datasets were entered. E.g. [0,0,1] to annotate only the third of 3
        datasets.
    :param mode: (str) display mode:
        If ``normal``, name of genes over 99%/under 1% quantile are displayed.
        If ``interactive``, click on a point to display its name.
        If ``json``, a .json file is produced that allows to reproduce the graph.
        in a web interface using Javascript.
    :param data_format: (str) ``counts`` or ``rpkm``.
    :param sep: (str) character delimiting the columns.
    :param limits: (list) bounds of the region displayed on the output graph: ``[min_x, max_x, min_y, max_y]``.
    :param slimits: (list) left and right bounds of the section of the splines to be displayed: ``[smin, smax]``
    :param deg: (int) the degree of the interpolating polynomial splines.
    :param bins: (int) the number of divisions of the x axis for quantiles estimation.
    :param assembly_id: (str or int) if an assembly ID is given,
        the json output will provide links to information on genes, using GenRep.
    :param normalize: (str) ``tags`` (total tag count) or ``sf`` (size factors).
    :param quantiles: (bool) if False, no quantile splines are drawn.
    :param title: (str) title to be written on top of the graph.
    :param extremes: (int) create an output file containing features for which ratios are outside the specified
        percentile (two-sided). The file is named *extreme_ratios_xxxxx* .
    """
    # Constants:
    if data_format == "counts":
        lower = 1 #lower bound on scores
        #bounds of what is taken into account for the computation of the splines
        spline_xmin = math.log10(math.sqrt(10)) #counts of 2,3 -> log(sqrt(2*3))
        spline_xmax = None
        #bounds of the section of the splines that is displayed
        slimits[0] = slimits[0] or 0.8
    elif data_format in ["rpkm","rpk"]:
        lower = 0
        spline_xmin = math.log10(0.1)
        spline_xmax = None
        slimits[0] = slimits[0] or -0.5
    min_pts_per_bin = 20
    prefix = prefix or unique_filename()

    # Extract data from CSV
    if isinstance(dataset,str): dataset = [dataset]
    if isinstance(labels,int): labels = [labels]
    allnames=[]; allmeans=[]; allratios=[]; points=[]
    groups={}; counts={}
    for data in dataset:
        cnts,names = read_data(data,sep,cols,labels,lower)
        cnts = asarray(cnts, dtype=float_)
        cnts = cnts[nonzero(cnts[:,0]*cnts[:,1])] # none of the counts is zero
        if normalize: cnts,sf = _normalize(cnts,mode=normalize)
        means = log10(sqrt(cnts[:,0]*cnts[:,1]))
        ratios = log2(cnts[:,0]/cnts[:,1])
        groups[data] = zip(names, means, ratios)
        counts[data] = dict(zip(names,cnts))
        allnames.extend(names); allmeans.extend(means); allratios.extend(ratios)

    points = zip(allnames, allmeans, allratios)
    try: xmin = min(allmeans); xmax = max(allmeans); ymin = min(allratios); ymax = max(allratios)
    except ValueError:
        raise ValueError("\nError: non-numeric columns selected. Try to specify more suitable columns with [-c].\n")

    # Figure initialization
    colors = iter(["black","red","green","cyan","magenta","yellow"])
    datacolors = dict((data,colors.next()) for data in dataset)
    ax = None
    if not mode == "json":
        if not mode == "interactive":
            import matplotlib
            matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=[14,9])
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.95)

        # Points
    for data in dataset:
        pts = zip(*groups[data])
        if ax: ax.plot(pts[1], pts[2], ".", color=datacolors[data])

    # Lines (best fit of percentiles)
    if quantiles:
        # Create bins; consider only a subset [spline_xmin, spline_xmax] of meaningful points.
        spline_xmin = spline_xmin or xmin
        spline_xmax = spline_xmax or xmax
        intervals = numpy.linspace(spline_xmin, spline_xmax, bins+1)

        points_in = []
        for i in range(len(intervals)-2,-1,-1): #from l-1 to 0, decreasing
            p_in_i = [p for p in points if p[1]>=intervals[i] and p[1]<intervals[i+1]]
            if len(p_in_i) < min_pts_per_bin:
                intervals = numpy.delete(intervals,i)
                bins-=1
            else:
                points_in.append(p_in_i)
        points_in.reverse()
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2. #the middle point of each bin

        # Compute percentiles in each bin
        spline_annotes=[]; spline_coords={}; extreme_ratios=[]
        percentiles = [1,5,25,50,75,95,99]
        if len(x) > 1: #if there are too few points there will be only one bin
            for k in percentiles:
                h=[]
                for b in range(bins):
                    score = stats.scoreatpercentile([p[2] for p in points_in[b]], k)
                    h.append(score)
                    if k == extremes: extreme_ratios.extend([p for p in points_in[b] if p[2] < score])
                    if k == 100-extremes: extreme_ratios.extend([p for p in points_in[b] if p[2] > score])
                coeffs = numpy.polyfit(x, h, deg)
                smin = slimits[0] or x[0]
                smax =  slimits[1] or 0.85*x[-1]
                x_spline = numpy.array(numpy.linspace(smin, smax, 10*bins))
                y_spline = numpy.polyval(coeffs, x_spline)

                if ax: ax.plot(x_spline, y_spline, "-", color="blue")
                #ax.plot(x, h, "o", color="blue") # testing
                spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
                spline_coords[k] = zip(x_spline,y_spline)

        # Print extreme ratios to file
        if extremes:
            extremes_filename = prefix + "_extreme_ratios.txt"
            with open(extremes_filename,"w") as f:
                c = csv.writer(f,delimiter="\t")
                c.writerow(["Name","countsC1","countsC2","log10Mean","log2Fold"])
                for p in extreme_ratios:
                    for data in dataset:
                        if counts[data].get(p[0]) is not None:
                            c.writerow([p[0], # label
                                        str(counts[data][p[0]][0]), str(counts[data][p[0]][1]), # raw counts
                                        str(p[1]), str(p[2])]) # mean & ratio
            print "Ratios r < %s %s < r : %s" % (extremes,100-extremes,extremes_filename)

        # Annotation of splines (percentage)
        if ax:
            for sa in spline_annotes:
                ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                            bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))

    if ax:
        # Decoration
        ax.set_xlabel("Log10 of sqrt(x1*x2)")
        ax.set_ylabel("Log2 of x1/x2")
        ax.set_title(title)
        if limits[0] is not None: xmin = limits[0]
        if limits[1] is not None: xmax = limits[1]
        if limits[2] is not None: ymin = limits[2]
        if limits[3] is not None: ymax = limits[3]
        xlen = abs(xmax-xmin); ylen = abs(ymax-ymin)
        plt.xlim([xmin-0.1*xlen,xmax+0.1*xlen])
        plt.ylim([ymin-0.1*ylen,ymax+0.1*ylen])

    # Annotation of points, draw
    annotes={}
    if not annotate:
        if mode == "normal": annotate = [0 for data in dataset]
        elif mode == "json": annotate = [1 for data in dataset]
    if mode == "interactive":
        af = AnnoteFinder(allmeans,allratios,allnames)
        plt.connect('button_press_event', af)
        plt.draw()
        plt.show()
    elif mode == "normal" or mode == "json":
        for data,annote in zip(dataset,annotate):
            if annote==1:
                annotes[data] = []
                for p in groups[data]:
                    if ax: ax.annotate(p[0], xy=(p[1],p[2]))
                    annotes[data].append(p[0])

    if mode == "normal":
        figname = prefix +".png"
        fig.savefig(figname)
        print "Figure:", figname
        return figname

    # Output for Javascript
    elif mode == "json":
        def rgb_to_hex(rgb):
            return '#%02x%02x%02x' % rgb
        jsdata = []
        # - data points
        for data in dataset:
            gdata = zip(*groups[data])
            datapts = zip(*gdata[1:3])
            jsdata.append({"label": "Data points",
                           "data": datapts,
                           "labels": annotes.get(data),
                           "points": {"symbol":"circle", "show":True},
                           "color": datacolors[data]})
        # - median (bold line)
        jsdata.append({"label": "Mean", "data": spline_coords[50],
                       "lines": {"show":True}, "color": rgb_to_hex((255,0,255)) })

        # - splines
        splinelabels = {"id": "Spline labels",
                        "data": [spline_coords[k][0] for k in percentiles],
                        "points": {"show":False},
                        "lines": {"show":False},
                        "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
        transparencies = iter([0.15, 0.2, 0.35, 0.35, 0.2, 0.15])
        jsdata.append({"id": "1%",
                       "data": spline_coords[1],
                       "lines": {"show":True, "lineWidth":0, "fill": transparencies.next()},
                       "color": rgb_to_hex((255,0,255))})
        percentiles.pop(percentiles.index(50))
        for k in range(1,len(percentiles)):
            jsdata.append({"id": str(percentiles[k])+"%",
                           "data": spline_coords[percentiles[k]],
                           "lines": {"show":True, "lineWidth":0, "fill": transparencies.next(),
                                     "fillBetween": str(percentiles[k-1])+"%"},
                           "color": rgb_to_hex((255,0,255))})

        jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
                 + "var splinelabels = " + json.dumps(splinelabels) + ";\n"

        print >> sys.stdout, jsdata
        return jsdata


#----------------------------------------------------------#

class AnnoteFinder:
  """
  Callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.

  Typical usage::

      plot(xdata, ydata)
      af = AnnoteFinder(xdata, ydata, annotes)
      connect('button_press_event', af)
  """

  def __init__(self, xdata, ydata, annotes, axis=None, xtol=None, ytol=None):
    self.data = zip(xdata, ydata, annotes)
    self.xrange = max(xdata) - min(xdata)
    self.yrange = max(ydata) - min(ydata)
    if xtol is None:
      xtol = len(xdata)*(self.xrange/float(len(xdata)))/10
    if ytol is None:
      ytol = len(ydata)*(self.yrange/float(len(ydata)))/10
    self.xtol = xtol
    self.ytol = ytol
    if axis is None:
      self.axis = plt.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """Distance between two points (x1,y1) and (x2,y2)"""
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

  def __call__(self, event):
    if event.inaxes:
      clickX = event.xdata
      clickY = event.ydata
      if self.axis is None or self.axis==event.inaxes:
        annotes = []
        for x,y,a in self.data:
          if  clickX-self.xtol < x < clickX+self.xtol and  clickY-self.ytol < y < clickY+self.ytol :
            annotes.append((self.distance(x/self.xrange,clickX/self.xrange,y/self.yrange,clickY/self.yrange),x,y, a) )
        if annotes:
          annotes.sort()
          distance, x, y, annote = annotes[0]
          self.drawAnnote(event.inaxes, x, y, annote)
          for l in self.links:
            l.drawSpecificAnnote(annote)

  def drawAnnote(self, axis, x, y, annote):
    """Draw on the plot the label of the point that was clicked"""
    if (x,y) in self.drawnAnnotations:
      markers = self.drawnAnnotations[(x,y)]
      for m in markers:
        m.set_visible(not m.get_visible())
      self.axis.figure.canvas.draw()
    else:
      t = axis.text(x,y, annote)
      m = axis.scatter([x],[y], marker='d', c='r', zorder=100)
      self.drawnAnnotations[(x,y)] =(t,m)
      self.axis.figure.canvas.draw()

  def drawSpecificAnnote(self, annote):
    annotesToDraw = [(x,y,a) for x,y,a in self.data if a==annote]
    for x,y,a in annotesToDraw:
      self.drawAnnote(self.axis, x, y, a)

