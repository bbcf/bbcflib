"""
======================
Module: bbcflib.rnaseq
======================

No documentation
"""

# Built-in modules #
import pickle, json, pysam, urllib, math, time
from itertools import combinations

# Internal modules #
from .mapseq import map_groups
from .genrep import GenRep

#started to take a census of used numpy libs
import numpy
#from numpy import *
from numpy.linalg import pinv 

# Other modules #
from bein.util import *
from scipy import stats
from scipy.interpolate import UnivariateSpline
from rpy2 import robjects
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
import rpy2.rlike.container as rlc
import cogent.db.ensembl as ensembl

################################################################################

def fetch_mappings(ex, path_or_assembly_id):
    """Given an assembly ID, return a tuple
    (gene_ids, gene_names, trans_to_gene_mapping, exons_to_trans_mapping).
    gene_ids is a list of gene IDs; gene_names is a dictionary {gene ID: gene name};
    trans_to_gene_mapping is a dictionary {gene ID: [transcripts IDs]};
    exons_to_trans_mapping is a dictionary {transcript ID: [exons IDs]}.

    *path_or_assembly_id* can be a numeric or nominal ID for GenRep
    (e.g. 76 or 'hg19' for H.Sapiens), or a path to a file containing a
    JSON object which is read to get the mapping.

    Takes about 115 seconds to load.
    """
    assembly_id = path = path_or_assembly_id
    nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
    nr_assemblies = json.loads(nr_assemblies)
    for a in nr_assemblies:
        if a['nr_assembly']['id'] == assembly_id or a['nr_assembly']['name'] == assembly_id:
            mdfive = a['nr_assembly']['md5']; break

    if os.path.exists(str(path)):
        with open(str(path), 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found in", os.path.dirname(str(path))
        return mapping
    elif os.path.exists("/db/genrep/nr_assemblies/exons_pickle/"+str(mdfive)+".pickle"):
        with open("/db/genrep/nr_assemblies/exons_pickle/"+mdfive+".pickle", 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found in /db/"
        return mapping
    else:
        genrep = GenRep('http://bbcftools.vital-it.ch/genrep/','/db/genrep/nr_assemblies/exons_pickle')
        assembly = genrep.assembly(mdfive)
        with open(assembly.index_path + '.pickle', 'rb') as pickle_file:
            mapping = pickle.load(pickle_file)
        print "Mapping found on GenRep"
        return mapping

def map_runs(fun, runs):
    """Parallelization of fun(run) executions"""
    futures = {}
    for group_id, run_list in runs.iteritems():
        futures[group_id] = [fun(run) for run in run_list]
    results = {}
    for group_id, future_list in futures.iteritems():
        results[group_id] = [f.wait() for f in future_list]
    return results

def exons_labels(bamfile):
    """List of the exons labels (SN: 'reference sequence name'; LN: 'sequence length') in *bamfile*."""
    sam = pysam.Samfile(bamfile, 'rb')
    labels = [(t['SN'],t['LN']) for t in sam.header['SQ']]
    sam.close()
    return labels

def pileup_file(bamfile, exons):
    """Return a numpy array of the pileup of *bamfile* in order *exons*."""
    counts = {}
    sam = pysam.Samfile(bamfile, 'rb')

    class Counter(object):
        def __init__(self):
            self.n = 0
        def __call__(self, alignment):
            self.n += 1

    c = Counter()
    for i,exon in enumerate(exons):
        sam.fetch(exon[0], 0, exon[1], callback=c) #(exon_name,0,exon_length,Counter())
        #The callback (c.n += 1) is executed for each alignment in a region
        counts[exon[0]] = c.n
        c.n = 0
    sam.close()
    return counts

def pairs_to_test(controls):
    """*controls* is a dictionary of group_ids to True/False.

    If all the values are True or all the values are False, then it
    returns all unique pairs of group IDs.  Otherwise it returns all
    combinations of one group ID with value True and one with value
    False.
    """
    if all(controls.values()) or not(any(controls.values())):
        return list(combinations(controls.keys(), 2))
    else:
        return [(x,y) for x in controls.keys() for y in controls.keys()
                if controls[x] and not(controls[y])]

@program
def external_deseq(cond1_label, cond1, cond2_label, cond2, method="normal", output=None, gene_names=None):
    if output:
        result_filename = output
    else:
        result_filename = unique_filename_in()
    c1 = unique_filename_in()
    with open(c1,'wb') as f:
        pickle.dump(cond1,f,pickle.HIGHEST_PROTOCOL)
    c2 = unique_filename_in()
    with open(c2,'wb') as f:
        pickle.dump(cond2,f,pickle.HIGHEST_PROTOCOL)
    gn = unique_filename_in()
    with open(gn,'wb') as f:
        pickle.dump(gene_names,f,pickle.HIGHEST_PROTOCOL)
    call = ["run_deseq.py", c1, c2, cond1_label, cond2_label, method, result_filename, gn]
    return {"arguments": call, "return_value": result_filename}

def inference(cond1_label, cond1, cond2_label, cond2, method="normal", output=None, gene_names=None):
    """Runs DESeq comparing the counts in *cond1* and *cond2*.

    Arguments:
    *cond1* and *cond2* are lists of numpy arrays. Each
    array lists the number of reads mapping to a particular
    transcript.
    *cond1_label* and *cond2_label* are string which will be used
    to identify the two conditions in R.
    *gene_names* is a dictionary {gene ID: gene name}.

    ``deseq_inference`` writes a tab delimited text file of the
    conditions to R, the first column being the transcript name,
    followed by one column for each numpy array in *cond1*, then one
    column for each numpy array in *cond2*.

    Then it calls DESeq in R, writes out the results in a new,
    randomly named file, and returns that filename.
    """

    # Pass the data into R as a data frame
    data_frame_contents = rlc.OrdDict([(cond1_label+'-'+str(i), robjects.IntVector(c.values()))
                                       for i,c in enumerate(cond1)] +
                                      [(cond2_label+'-'+str(i), robjects.IntVector(c.values()))
                                       for i,c in enumerate(cond2)])
    data_frame = robjects.DataFrame(data_frame_contents)
    data_frame.rownames = cond1[0].keys()
    conds = robjects.StrVector([cond1_label for x in cond1] + [cond2_label for x in cond2]).factor()
    if cond1[0].keys() != cond2[0].keys(): raise

    ## DESeq full analysis
    deseq = rpackages.importr('DESeq')
    cds = deseq.newCountDataSet(data_frame, conds)
    cds = deseq.estimateSizeFactors(cds)
    try: cds = deseq.estimateVarianceFunctions(cds,method=method)
    except : raise rpy2.rinterface.RRuntimeError("Too few reads to estimate variances with DESeq")
    res = deseq.nbinomTest(cds, cond1_label, cond2_label)

    ## Replace (unique) gene IDs by (not unique) gene names
    if gene_names:
        names = []
        res_ids = list(res[0])
        if res_ids[0].find("ENSG") != -1:
            for s in res_ids:
                start = s.find("ENSG")
                end = s.split("ENSG")[1].find("|")
                gene_id = "ENSG" + s.split("ENSG")[1].split("|")[0]
                names.append(s.replace(gene_id,gene_names.get(gene_id,gene_id)))
        data_frame = robjects.DataFrame({"Name":robjects.StrVector(names)})
        for i in range(2,len(res)+1):
            data_frame = data_frame.cbind(res.rx(i))
        res = data_frame

    if output: result_filename = output
    else: result_filename = unique_filename_in()

    res.to_csvfile(result_filename)
    return result_filename

def rnaseq_workflow(ex, job, assembly, target=["genes"], via="lsf", output=None, maplot="normal"):
    """Run RNASeq inference according to *job_info*.

    *output*: alternative name for output file. Otherwise it is random.
    *maplot*: MA-plot of data.
    - If 'interactive', one can click on a point (gene or exon) to display its name;
    - if 'normal', name of genes over 99.9%/under 0.1% quantiles are displayed;
    - if None, no figure is produced.
    *target*: a string or array of strings indicating the targets to run inference (DESeq) on.
    Targets can be 'genes', 'transcripts', 'exons'. E.g. ['genes','transcripts'], or 'exons'.
    (This part of the workflow may fail if there are too few reads for DESeq to estimate variances.)
    *job* is as Job object (or a dictionary of the same form) as returned from
    HTSStation's frontend.

    Whatever script calls this workflow needs to pass back the JSON string
    it returns in some sensible way.  For the usual HTSStation
    frontend, this just means printing it to stdout.
    """

    names = {}; runs = {}; controls = {}; paths = {}; gids = {}
    groups = job.groups
    assembly_id = job.assembly_id
    for i,group in groups.iteritems():
        names[i] = str(group['name'])
        runs[i] = group['runs'].values()
        controls[i] = group['control']
    if isinstance(target,str): target=[target]

    fastq_root = os.path.abspath(ex.working_directory)
    print "Alignment..."
    bam_files = map_groups(ex, job, fastq_root, assembly_or_dict = assembly)
    print "Reads aligned."

    # All the bam_files were created against the same index, so
    # they all have the same header in the same order.  I can take
    # the list of exons from just the first one and use it for all of them.
    # Format: ('exonID|geneID|start|end|strand', length)
    exons = exons_labels(bam_files.values()[0][1]['bam'])

    print "Loading mappings..."
    # - gene_ids is a list of gene ID's
    # - gene_names is a dict {gene ID: gene name}
    # - transcript_mapping is a dictionary {transcript ID: gene ID}
    # - exon_mapping is a dictionary {exon ID: (transcript ID, gene ID)}
    # - trans_in_gene is a dict {transcript ID: ID of the gene it belongs to}
    # - exons_in_trans is a dict {exon ID: ID of the transcript it belongs to}
    mappings = fetch_mappings(ex, assembly_id)
    (gene_ids, gene_names, transcript_mapping, exon_mapping, trans_in_gene, exons_in_trans) = mappings
    print "Loaded."
    
    exon_pileups = {}; gene_pileups = {}; tran_pileups = {}
    for condition,files in bam_files.iteritems():
        gene_pileups[condition] = []
        exon_pileups[condition] = []
        tran_pileups[condition] = []
        unknown_exons = []
        for f in files.values():
            exon_pileup = pileup_file(f['bam'], exons) #{exon_id: count}
            exon_pileups[condition].append(exon_pileup) #{cond1: [{pileup bam1},{pileup bam2},...], cond2:...}
            gene_pileup = dict(zip(gene_ids,numpy.zeros(len(gene_ids))))
            tran_pileup = dict(zip(transcript_mapping.keys(),numpy.zeros(len(transcript_mapping.keys()))) )
            for e,c in exon_pileup.iteritems():
                e = e.split('|')[0]
                if exon_mapping.get(e):
                    gene_pileup[exon_mapping[e][1]] += c
            gene_pileups[condition].append(gene_pileup)            
            ### Get counts for the transcripts using pseudo-inverse
            for g in gene_ids:
                tg = trans_in_gene[g]
                if len(tg) == 1:
                    t = tg[0]
                    eg = exons_in_trans[t]
                    for e in eg:
                        if exon_pileup.get(e):
                            tran_pileup[t] += exon_pileup[e]
                else:   
                    eg = []
                    for t in tg:
                        ebt = exons_in_trans[t]
                        eg.extend(ebt)
                    M = numpy.zeros((len(eg),len(tg)))
                    exons_expression = numpy.zeros(len(eg))
                    for i,e in enumerate(eg):
                        for j,t in enumerate(tg):
                            ebt = exons_in_trans[t]
                            if e in ebt:
                                M[i,j] = 1
                        if exon_pileup.get(e):
                            exons_expression[i] += exon_pileup[e]
                    transcripts_expression_g = numpy.ceil(numpy.dot(pinv(M),exons_expression)) #pseudo-inverse
                    # DESeq accepts integers only
                    #transcripts_expression_g = numpy.dot(pinv(M),exons_expression)
                    for k,t in enumerate(tg):
                        tran_pileup[t] = transcripts_expression_g[k]
            tran_pileups[condition].append(tran_pileup)

    futures_genes = {}; futures_trans = {}; futures_exons = {}
    for (c1,c2) in pairs_to_test(controls):
        if len(runs[c1]) + len(runs[c2]) > 2:
            method = "normal"
        else: method = "blind"

        print "Inference..."
        if "genes" in target:
            futures_genes[(c1,c2)] = external_deseq.nonblocking(ex,
                                      names[c1], gene_pileups[c1],
                                      names[c2], gene_pileups[c2],
                                      method, output, gene_names=gene_names, via=via)
        if "transcripts" in target:
            futures_trans[(c1,c2)] = external_deseq.nonblocking(ex,
                                      names[c1], tran_pileups[c1],
                                      names[c2], tran_pileups[c2],
                                      method, output, gene_names=gene_names, via=via)
        if "exons" in target:
            futures_exons[(c1,c2)] = external_deseq.nonblocking(ex,
                                      names[c1], exon_pileups[c1],
                                      names[c2], exon_pileups[c2],
                                      method, output, gene_names=gene_names, via=via)
        print "....Wait for results....."
        for c,f in futures_genes.iteritems():
            gene_file = f.wait()
            if gene_file:
                ex.add(gene_file,
                       description="Comparison of GENES in conditions '%s' and '%s' (CSV)" % (names[c[0]], names[c[1]]))
                if maplot:
                    res = robjects.DataFrame.from_csvfile(f.wait())
                    figname, jsname = MAplot(res, mode=maplot, deg=4, bins=100, alpha=0.005, assembly_id=assembly_id)
                    ex.add(figname, description="MA-plot (genes)")
                    ex.add(jsname, description="Data for javascript (genes)")
                print "Genes: Done successfully."
            else: print "Genes: Failed during inference, probably because of too few reads for DESeq stats."

        for c,f in futures_trans.iteritems():
            trans_file = f.wait()
            if trans_file:
                ex.add(trans_file,
                       description="Comparison of TRANSCRIPTS in conditions '%s' and '%s' (CSV)" % (names[c[0]], names[c[1]]))
                if maplot:
                    res = robjects.DataFrame.from_csvfile(f.wait())
                    figname, jsname = MAplot(res, mode=maplot, deg=4, bins=100, alpha=0.005, assembly_id=assembly_id)
                    ex.add(figname, description="MA-plot (transcripts)")
                    ex.add(jsname, description="Data for javascript (transcripts)")
                print "Transcripts: Done successfully."
            else: print "Transcripts: Failed during inference, probably because of too few reads for DESeq stats."

        for c,f in futures_exons.iteritems():
            exons_file = f.wait()
            if exons_file:
                ex.add(exons_file,
                       description="Comparison of EXONS in conditions '%s' and '%s' (CSV)" % (names[c[0]], names[c[1]]))
                if maplot:
                    res = robjects.DataFrame.from_csvfile(f.wait())
                    figname, jsname = MAplot(res, mode=maplot, deg=4, bins=100, alpha=0.005, assembly_id=assembly_id)
                    ex.add(figname, description="MA-plot (exons)")
                    ex.add(jsname, description="Data for javascript (exons)")
                print "Exons: Done successfully."
            else: print "Exons: Failed during inference, probably because of too few reads for DESeq stats."


def MAplot(data, mode="interactive", deg=4, bins=30, alpha=0.005, assembly_id=None):
    """
    Creates an "MA-plot" to compare transcription levels of a set of genes
    in two different conditions. It returns the name of the .png file produced,
    and the name of a json containing enough information to reconstruct the plot using Javascript.

    data:  rpy DataFrame object; two columns, each for a different condition.
    mode:  display mode;
    - if "interactive", click on a point to display its name
    - if "normal", name of genes over 99%/under 1% quantile are displayed
    alpha: a threshold for p-values: points beyond it are emphasized.
    assembly_id  : if an assembly ID is given, the json output will provide links to information on genes.
    """

    import matplotlib.pyplot as plt

    ## Extract data from DataFrame
    pvals = numpy.asarray(data.rx("pval")[0])
    means = numpy.asarray(data.rx("baseMean")[0])
    ratios = numpy.asarray(data.rx("log2FoldChange")[0])
    names = data.rx("Name")[0]
    a = numpy.asarray(names)
    b = numpy.asarray(names.levels)
    names = list(b[a-1])

    points=[]; blackpoints=[]; redpoints=[]; annotes_red=[]; annotes_black=[]; pvals_red=[]; pvals_black=[]
    for i in range(len(ratios)):
        if not math.isnan(ratios[i]) and not math.isinf(ratios[i]):
            p = (names[i],numpy.log10(means[i]),ratios[i],pvals[i])
            if pvals[i] < alpha:
                redpoints.append(p)
                annotes_red.append(p[0])
                pvals_red.append(p[3])
            else:
                blackpoints.append(p)
                annotes_black.append(p[0])
                pvals_black.append(p[3])
            points.append(p)
    names,means,ratios,pvals = zip(*points)
    xmin = min(means); xmax = max(means); ymin = min(ratios); ymax = max(ratios)

    ## Create bins
    N = len(points); dN = N/bins #points per bin
    rmeans = numpy.sort(means)
    intervals = []
    for i in range(bins):
        intervals.append(rmeans[i*dN])
    intervals.append(xmax)
    intervals = numpy.array(intervals)

    points_in = {}; perc = {}
    for b in range(bins):
        points_in[b] = [p for p in points if p[1]>=intervals[b] and p[1]<intervals[b+1]]
        perc[b] = [p[2] for p in points_in[b]]

    ## Figure
    fig = plt.figure(figsize=[14,9])
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.98)
    figname = None

    ### Points
    blackpts = zip(*blackpoints)[1:3]
    redpts = zip(*redpoints)[1:3]
    ax.plot(blackpts[0], blackpts[1], ".", color="black")
    if len(redpts) != 0:
        ax.plot(redpts[0], redpts[1], ".", color="red")

    ### Lines (best fit of percentiles)
    spline_annotes=[]; spline_coords={}
    percentiles = [1,5,25,50,75,95,99]
    for k in percentiles:
        h = numpy.ones(bins)
        for b in range(bins):
            if points_in[b] != []:
                h[b] = stats.scoreatpercentile(perc[b], k)
            else: h[b] = h[b-1]
        x = intervals[:-1]+(intervals[1:]-intervals[:-1])/2.
        spline = UnivariateSpline(x, h, k=deg)
        xs = numpy.array(numpy.linspace(xmin, xmax, 10*bins)) #to increase spline smoothness
        ys = numpy.array(spline(xs))
        l = len(xs)
        xi = numpy.arange(l)[numpy.ceil(l/6):numpy.floor(8*l/9)]
        x_spline = xs[xi]
        y_spline = ys[xi]
        ax.plot(x_spline, y_spline, "-", color="blue") #ax.plot(x, h, "o", color="blue")
        spline_annotes.append((k,x_spline[0],y_spline[0])) #quantile percentages
        spline_coords[k] = zip(x_spline,y_spline)

    ### Decoration
    ax.set_xlabel("Log10 of sqrt(x1*x2)")
    ax.set_ylabel("Log2 of x1/x2")
    for sa in spline_annotes:
        ax.annotate(str(sa[0])+"%", xy=(sa[1],sa[2]), xytext=(-33,-5), textcoords='offset points',
                    bbox=dict(facecolor="white",edgecolor=None,boxstyle="square,pad=.4"))
    if mode == "interactive":
        af = AnnoteFinder( means, ratios, names )
        plt.connect('button_press_event', af)
        plt.draw()
        plt.show()
    else:
        for p in redpoints:
            ax.annotate(p[0], xy=(p[1],p[2]) )
    figname = unique_filename_in()+".png"
    fig.savefig(figname)

    ## Output for Javascript
    def rgb_to_hex(rgb):
        return '#%02x%02x%02x' % rgb
    redpoints = [(p[1],p[2]) for p in redpoints]
    blackpoints = [(p[1],p[2]) for p in blackpoints]
    jsdata = [{"label": "Genes with p-value < " + str(alpha),
               "data": redpoints,
               "labels": annotes_red,
               "pvals": pvals_red,
               "points": {"symbol":"circle", "show":True},
               "color": "red"},
              {"label": "Genes with p-value > " + str(alpha),
               "data": blackpoints,
               "labels": annotes_black,
               "pvals": pvals_black,
               "points": {"symbol":"circle", "show":True},
               "color": "black"},
              {"label": "Mean", "data": spline_coords[50],
               "lines": {"show":True}, "color": rgb_to_hex((255,0,255)) },
              {"id": "1% Quantile", "data": spline_coords[1], "lines": {"show":True, "lineWidth":0, "fill": 0.12},
               "color": rgb_to_hex((255,0,255))},
              {"id": "5% Quantile", "data": spline_coords[5], "lines": {"show":True, "lineWidth":0, "fill": 0.18},
               "color": rgb_to_hex((255,0,255))},
              {"id": "25% Quantile", "data": spline_coords[25], "lines": {"show":True, "lineWidth":0, "fill": 0.3},
               "color": rgb_to_hex((255,0,255))},
              {"id": "75% Quantile", "data": spline_coords[75], "lines": {"show":True, "lineWidth":0, "fill": 0.3},
               "color": rgb_to_hex((255,0,255))},
              {"id": "95% Quantile", "data": spline_coords[95], "lines": {"show":True, "lineWidth":0, "fill": 0.18},
               "color": rgb_to_hex((255,0,255))},
              {"id": "99% Quantile", "data": spline_coords[99], "lines": {"show":True, "lineWidth":0, "fill": 0.12},
               "color": rgb_to_hex((255,0,255))}
             ]
    splinelabels = {"id": "Spline labels",
                    "data": [spline_coords[k][0] for k in percentiles],
                    "points": {"show":False}, "lines": {"show":False},
                    "labels": ["1%","5%","25%","50%","75%","95%","99%"]}
    jsdata = "var data = " + json.dumps(jsdata) + ";\n" \
             + "var splinelabels = " + json.dumps(splinelabels) + ";\n"
    if assembly_id:
        nr_assemblies = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies.json").read()
        nr_assemblies = json.loads(nr_assemblies)
        if isinstance(assembly_id,str):
            for a in nr_assemblies:
                if a['nr_assembly']['name'] == assembly_id:
                    assembly_id = a['nr_assembly']['id']; break
        url_template = urllib.urlopen("http://bbcftools.vital-it.ch/genrep/nr_assemblies/" \
                                      + str(assembly_id) + "/get_links.json?gene_name=%3CName%3E")
        jsdata = jsdata + "var url_template = " + url_template.read() + ";"
    jsname = unique_filename_in()+".js"
    jsname = "data.js"
    with open(jsname,"w") as js:
        js.write(jsdata)

    return figname, jsname

class AnnoteFinder:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.

  Use this function like this:

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
      self.axis = pylab.gca()
    else:
      self.axis= axis
    self.drawnAnnotations = {}
    self.links = []

  def distance(self, x1, x2, y1, y2):
    """distance between two points"""
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
    """
    Draw the annotation on the plot
    """
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
