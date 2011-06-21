#source $HOME/.bashrc

import re
import sys
import os
import random
from ..rnaseq import rnaseq_workflow
from unittest2 import TestCase, TestSuite, main, TestLoader, skipIf


class TestWorkflow(TestCase):
    def test_workflow(self):
        root = "../../extras/reads/"
        if not os.exists(os.path.join(root,"s_5_10k.txt")) or not os.exists(os.path.join(root,"s_6_10k.txt")):
            self.skipTest("Files required for RNASeq unittest are missing. Test skipped.")
        data = {"assembly_id": 76,
                "groups": {1: {'name': "TestA", "runs": {'1':os.path.join(root,"s_5_10k.txt")}, 'control': True},
                           2: {'name': "TestB", "runs": {'1':os.path.join(root,"s_6_10k.txt")}, 'control': False}
                          },
                "options": {}
               }
        output = rnaseq_workflow(data, job_or_dict="dict", via="lsf", maplot="normal")

main()


### >>> OUTPUT <<<
## 0|141|4345||0||Warning: Exhausted best-first chunk memory for read C3PO_0038:5:1:16904:1821#0/1 (patid 11655); skipping read
## Warning: Exhausted best-first chunk memory for read C3PO_0038:5:1:6861:2462#0/1 (patid 21325); skipping read
## Warning: Exhausted best-first chunk memory for read C3PO_0038:5:1:10058:2579#0/1 (patid 23118); skipping read
## # reads processed: 25000
## # reads with at least one reported alignment: 15321 (61.28%)
## # reads that failed to align: 9679 (38.72%)
## Reported 44393 alignments to 1 output stream(s)

## 1|141|4347||0||Warning: Exhausted best-first chunk memory for read C3PO_0038:6:1:8405:1045#0/1 (patid 140); skipping read
## Warning: Exhausted best-first chunk memory for read C3PO_0038:6:1:4676:2034#0/1 (patid 13977); skipping read
## Warning: Exhausted best-first chunk memory for read C3PO_0038:6:1:4355:2457#0/1 (patid 20078); skipping read
## # reads processed: 25000
## # reads with at least one reported alignment: 15144 (60.58%)
## # reads that failed to align: 9853 (39.41%)
## # reads with alignments suppressed due to -m: 3 (0.01%)
## Reported 41511 alignments to 1 output stream(s)

## 2|141|4596||0||PyTables not found.  Skipping.

## 3|141|4598||0||PyTables not found.  Skipping.

## 4|141|7300||0||
## 5|141|7302||0||
## 6|141|7564||0||
## 7|141|7566||0||
## 8|141|8398||0|locfit 1.5-6 	 2010-01-20 
## |PyTables not found.  Skipping.

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
