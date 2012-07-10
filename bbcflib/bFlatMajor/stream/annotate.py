import sys
from bbcflib.bFlatMajor import common
from bbcflib import btrack as track

def getNearestFeature(features, annotations,
                      thresholdPromot=2000, thresholdInter=100000, thresholdUTR=10):
    """
    For each element of *features*, searches the nearest element of *annotations* and returns
    a stream similar to *features*, with additional annotation fields, e.g.::

        ('chr5',12,14) -> ('chr5',12,14,'geneId|geneName','location_type','distance').

    If there are several genes, they are separated by '_': geneId1|geneName1_geneId2|geneName2.
    For each gene, `location_type` is one of:
    * 'Intergenic' if there are no genes within a distance `thresholdInter`,
    * 'Included' if the feature is included in the gene,
    * 'Promot' if the feature is upstream and within `thresholdInter` of the gene start,
    * 'Upstream' if the feature is upstream and beyond the promoter of the gene,
    * '3UTR' if the feature is downstream and within `thresholdUTR`% of the distance to the next downstream gene,
    * 'Downstream' otherwise.
    These annotations can be concatenated with '_' as well.
    The distance to each gene is negative if the feature is included, positive otherwise.

    :param features: (bbcflib.btrack.FeatureStream) features track.
    :param annotations: (bbcflib.btrack.FeatureStream) gene annotation track
        (e.g. as obtained with assembly.gene_track()).
    :param thresholdPromot: (int) associates the promoter of each gene which promoter is within
        this distance of the feature. Above the threshold, associates only the closest. [2000]
    :param thresholdInter: (int) no gene beyond this distance will be considered. [100000]
    :param thresholdUTR: (int) in case the feature is surrounded by two eligible genes on the
        same strand: if distance to gene1's 3'UTR upstream is less than *thresholdUTR*% of the distance
        between gene1 and gene2, associated to 3'UTR of gene1, else to promoter of gene2. [10]
    :rtype: bbcflib.btrack.FeatureStream (..., str, str, str).

    ::

                 <--                   feat                    -->
             ______| thresholdPromot  ++++++   thresholdPromot |______
       -----|______|-------------------------------------------|______|----------
             gene 1                       gene 2

                                         feat
             ______  thresholdInter     ++++++        thresholdInter   ______
       -----|______|----------...------------------...----------------|______|---
             gene 1                                         gene 2

                          feat
            -->          ++++++               -->
            |______  10%             90%     |______
       -----|______|------|------------------|______|-----  (attributed to gene1)
             gene 1      thresholdUTR     gene 2
    """
    def _get_feature(_t,_a):
        F = []
        _a = common.sentinelize(_a, [sys.maxint]*len(_a.fields))
        for peak in _t:
            distMinBefore = distMinAfter = thresholdInter+1
            gene = ""; dist = ""; typeLoc = ""
            geneBefore = None; geneAfter = None
            strandBefore = None; strandAfter = None
            included = 0
            # keep only genes which don't start too far
            for annot in _a:
                F.append(annot)
                if annot[0] > peak[1]+thresholdInter: break
            # remove genes that end too far
            fpop = -1 # always keep one gene before
            for annot in F:
                if annot[1] > peak[0]-thresholdInter: break
                fpop += 1
            if fpop>0: F = F[fpop:]
            for annot in F:
                # if the peak is totally included in the gene
                if (peak[0]>=annot[0]) and (annot[1]>=peak[1]):
                    includedGene = annot[2]
                    includedDist = (annot[3] == -1) and annot[1]-peak[1] or peak[0]-annot[0]
                    included = 1
                # if the gene is totally included in the peak
                elif (annot[0]>peak[0]) and (peak[1]>annot[1]):
                    includedGene = annot[2]
                    includedDist = 0
                    included = 1
                else:
                    # if annot is not too far 3' and no intersection
                    if (peak[0]-annot[1]<distMinBefore) and (peak[0]-annot[1]>0):
                        distMinBefore = peak[0]-annot[1]
                        geneBefore = annot[2]
                        strandBefore = annot[3]
                    # if intersection (annot is before)
                    elif (peak[0]-annot[1]<0) and (peak[0]-annot[0]>0):
                        distMinBefore = 0
                        geneBefore = annot[2]
                        strandBefore = annot[3]
                        #print "gene %s overlaps begin of peak %s" % (geneBefore,peakName)
                    # if annot is not too far 5' and no intersection
                    if (annot[0]-peak[1]<distMinAfter) and (annot[0]-peak[1]>0):
                        distMinAfter = annot[0]-peak[1]
                        geneAfter = annot[2]
                        strandAfter = annot[3]
                    # if intersection (annot is after)
                    elif (annot[0]-peak[1]<0) and (annot[1]-peak[1]>0):
                        distMinAfter = 0
                        geneAfter = annot[2]
                        strandAfter = annot[3]
                        #print "gene %s overlaps end of peak %s" % (geneAfter,peakName)
            # detect intergenic peak
            if not(included) and distMinBefore>thresholdInter and distMinAfter>thresholdInter:
                yield peak+('','Intergenic',thresholdInter)
                continue
            # detect peak before the first chromosome gene or after the last chromosome gene
            if geneBefore == None:
                if distMinAfter<=thresholdInter:
                    gene = geneAfter
                    dist = distMinAfter
                    typeLoc = (strandAfter == 1) and  "Upstream" or "Downstream"
            elif geneAfter == None:
                if distMinBefore<=thresholdInter:
                    gene = geneBefore
                    dist = distMinBefore
                    typeLoc = (strandBefore == -1) and  "Upstream" or "Downstream"
            # detect peak between two genes on the same strand
            elif strandBefore == strandAfter:
                if strandBefore == 1:
                    if thresholdUTR*distMinAfter > 100*distMinBefore:
                        gene = geneBefore
                        dist = distMinBefore
                        typeLoc = "3UTR"
                    else:
                        gene = geneAfter
                        dist = distMinAfter
                        typeLoc = "Promot"
                else:
                    if thresholdUTR*distMinBefore > 100*distMinAfter:
                        gene = geneAfter
                        dist = distMinAfter
                        typeLoc = "3UTR"
                    else:
                        gene = geneBefore
                        dist = distMinBefore
                        typeLoc = "Promot"
            # detect peak between two genes on different strands
            else:
                # detect peak between 2 promoters
                if strandBefore == -1:
                    typeLoc = "Promot"
                    # detect peak close enough to both promoters
                    if distMinBefore < thresholdPromot and distMinAfter < thresholdPromot:
                        typeLoc = typeLoc+"_Promot"
                        gene = geneBefore+"_"+geneAfter
                        dist = str(distMinBefore)+"_"+str(distMinAfter)
                    # detect peak closer to geneBefore
                    elif distMinBefore < distMinAfter:
                        gene = geneBefore
                        dist = distMinBefore
                    # detect peak closer to geneAfter
                    else:
                        gene = geneAfter
                        dist = distMinAfter
                # detect peak between 2 3UTR
                else:
                    typeLoc = "3UTR"
                    # detect peak overlapping the 2 3UTR
                    if distMinBefore == distMinAfter:
                        typeLoc = typeLoc+"_3UTR"
                        gene = geneBefore+"_"+geneAfter
                        dist = str(distMinBefore)+"_"+str(distMinAfter)
                    elif distMinBefore < distMinAfter:
                        gene = geneBefore
                        dist = distMinBefore
                    else:
                        gene = geneAfter
                        dist = distMinAfter
            if included == 1:
                gene = gene and gene+"_"+includedGene or includedGene
                dist = dist and str(dist)+"_"+str(includedDist) or str(includedDist)
                typeLoc = typeLoc and typeLoc+"_Included" or "Included"
            yield peak+(gene,typeLoc,dist)
    if isinstance(features,(tuple,list)): features = features[0]
    if isinstance(annotations,(tuple,list)): annotations = annotations[0]
    features = common.reorder(features,['start','end'])
    annot = common.reorder(annotations,['start','end','name','strand'])
    _fields = features.fields+['gene','location_type','distance']
    return track.FeatureStream(_get_feature(features,annot),fields=_fields)


