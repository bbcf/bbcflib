import sys
from bbcflib.bFlatMajor import common
from bbcflib import btrack as track

def getNearestFeature(features, annotations,
                      thresholdPromot=2000, thresholdInter=100000, thresholdUTR=10):
    def _get_feature(_t,_a):
        F = []
        _a = common.sentinelize(_a, [sys.maxint]*len(_a.fields))
        for peak in _t:
            distMinBefore = distMinAfter = thresholdInter+1
            gene = dist = typeLoc = ""
            geneBefore = geneAfter = None
            included=0
            for x in _a:
                F.append(x)
                if x[0] > peak[1]+thresholdInter: break
            fpop = 0
            for x in F:
                if x[1] < peak[0]-thresholdInter: break
                fpop += 1
            F = F[fpop:]
            for annot in F:
                # if the peak is totaly included in the gene
                if (peak[0]>annot[0]) and (annot[1]>peak[1]):
                    includedGene = annot[2]
                    includedDist = peak[0]-annot[1]
                    included = 1
                # if the gene is totaly included in the peak
                elif (annot[0]>peak[0]) and (peak[1]>annot[1]): 
                    includedGene = annot[2]
                    includedDist = 0
                    included = 1
                else:
                    if (peak[0]-annot[1]<distMinBefore) and (peak[0]-annot[1]>0):
                        distMinBefore = peak[0]-annot[1]
                        geneBefore = annot[2]
                        strandBefore = annot[3]
                    elif (peak[0]-annot[1]<0) and (peak[0]-annot[0]>0):
                        distMinBefore = 0
                        geneBefore = annot[2]
                        strandBefore = annot[3]
                    #print "gene %s overlap begin of peak %s" % (geneBefore,peakName)
                    if (annot[0]-peak[1]<distMinAfter) and (annot[0]-peak[1]>0):
                        distMinAfter = annot[0]-peak[1]
                        geneAfter = annot[2]
                        strandAfter = annot[3]
                    elif (annot[0]-peak[1]<0) and (annot[1]-peak[1]>0):
                        distMinAfter = 0
                        geneAfter = annot[2]
                        strandAfter = annot[3]
                    #print "gene %s overlap end of peak %s" % (geneAfter,peakName)
                    # detect intergenic peak
            if distMinBefore>thresholdInter and distMinAfter>thresholdInter:
                yield peak+('','Intergenic',thresholdInter)
    
            else:
                # detect peak before the first chromosome gene or after the last chromosome gene
                if geneBefore == None:
                    gene = geneAfter
                    dist = distMinAfter
                elif geneAfter == None:
                    gene = geneBefore
                    dist = distMinBefore
                # detect peak between two genes on the same strand
                elif strandBefore == strandAfter:
                    if strandBefore==1:
                        if ((thresholdUTR*distMinAfter)/100)>distMinBefore :
                            gene = geneBefore
                            dist = distMinBefore
                            typeLoc = "3UTR"
                        else:
                            typeLoc = "Promot"
                            gene = geneAfter
                            dist = distMinAfter
                    else:
                        if ((thresholdUTR*distMinBefore)/100)>distMinAfter :
                            gene = geneAfter
                            dist = distMinAfter
                            typeLoc = "3UTR"
                        else:
                            typeLoc = "Promot"
                            gene = geneBefore
                            dist = distMinBefore

                else:
                    # detect peak between 2 promoters
                    if strandBefore == -1:
                        typeLoc = "Promot"
                        # detect peak near the 2 promoters
                        if distMinBefore < thresholdPromot and distMinAfter < thresholdPromot:
                            typeLoc = typeLoc+"_Promot"
                            gene = geneBefore+"_"+geneAfter
                            dist = str(distMinBefore)+"_"+str(distMinAfter)
                        # detect peak closest than geneBefore
                        elif distMinBefore < distMinAfter:
                            gene = geneBefore
                            dist = distMinBefore
                        # detect peak closest than geneAfter
                        else:
                            gene = geneAfter
                            dist = distMinAfter
                    # detect peak between 2 3UTR
                    else:
                        typeLoc = "3UTR"
                        # detect peak overlaping the 2 3UTR
                        if distMinBefore == distMinAfter:
                            typeLoc=typeLoc+"_3UTR"
                            gene = geneBefore+"_"+geneAfter
                            dist = str(distMinBefore)+"_"+str(distMinAfter)
                        elif distMinBefore < distMinAfter:
                            gene = geneBefore
                            dist = distMinBefore
                        else:
                            gene = geneAfter
                            dist = distMinAfter
                    
                if included==1:
                    included = 0
                    gene = gene+"_"+includedGene
                    dist = str(dist)+"_"+str(includedDist)
                    typeLoc = typeLoc+"_Included"
                yield peak+(gene,typeLoc,dist)
    
    if isinstance(features,(tuple,list)): features = features[0]
    if isinstance(annotations,(tuple,list)): annotations = annotations[0]
    features = common.reorder(features,['start','end'])
    annot = common.reorder(annotations,['start','end','gene_name','strand'])
    _fields = features.fields+['gene','location_type','distance']
    return track.FeatureStream(_get_feature(features,annot),fields=_fields)


