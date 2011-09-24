# Built-in modules #
import datetime, ConfigParser, cStringIO

# Internal modules #
from bbcflib.rnaseq import *

# Unitesting module #
import unittest2

# Nosetest flag #
__test__ = True


class Test_Rnaseq(unittest2.TestCase):
    def setUp(self):
	self.counts = {
            "e1":numpy.array([12,3]),
	    "e2":numpy.array([12,3]),
	    "e3":numpy.array([12,3]),
  	    "e4":numpy.array([12,3]),
	    "e5":numpy.array([12,3]),
	    "e6":numpy.array([12,3])}
	self.gene_ids = ["g1","g2"]
	self.exon_mapping = {
	    "e1":(["t1","t2","t3"],"g1"), "e2":(["t1","t2","t3"],"g1"),
	    "e3":(["t2","t3","t4"],"g1"), "e4":(["t3","t4"],"g1"),
	    "e5":(["t4"],"g1"), "e6":(["t"],"g2")}
	self.transcript_mapping = {
	    "t1":"g1", "t2":"g1", "t3":"g1", "t4":"g1", "t":"g2"}
	self.trans_in_gene = {
	    "g1":["t1","t2","t3","t4"], "g2":["t"]}
	self.exons_in_trans = {
	    "t1":["e1","e2"], "t2":["e1","e2","e3"], 
	    "t3":["e1","e2","e3","e4"], "t4":["e3","e4","e5"], 
            "t":["e6"]}
	"""
	|==================g1==================| |===g2===|
	|...e1...| |.e2.| |..e3..| |.e4.| |.e5.| |...e6...|
	|-------t1------|
	|------------t2----------|
	|---------------t3--------------|
	                  |---------t4---------| |----t---|
	"""
    def test_transcripts_expression(self):
	texp = transcripts_expression(self.gene_ids, self.transcript_mapping,
		self.trans_in_gene, self.exons_in_trans, self.counts)
    def test_genes_expression(self):
	gexp = genes_expression(self.gene_ids, self.exon_mapping, self.counts)
	self.assertListEqual(list(gexp["g1"]), [5*12,5*3])
	self.assertListEqual(list(gexp["g2"]), [12,3])
    def test_estimate_size_factors(self):
	counts = numpy.array([[12,12,12],[3,3,3]])
	res, size_factors = estimate_size_factors(counts)
	self.assertListEqual(list(size_factors),[2,0.5])
	self.assertListEqual(list(res[0]),[6,6,6])
	self.assertListEqual(list(res[1]),[6,6,6])

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
