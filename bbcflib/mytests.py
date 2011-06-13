#source $HOME/.bashrc

import re
import sys
import os
import random
import rnaseq
from unittest2 import TestCase, TestSuite, main, TestLoader, skipIf


class TestWorkflow(TestCase):
    def test_workflow(self):
        #root = "/scratch/cluster/monthly/jrougemo/rnaseq-test/fastq_files/"
        root = "/scratch/cluster/monthly/jdelafon-el/"
        #root = "."
        data = {"assembly_id": 11,
                "groups": {1: {'name': "TestA", "runs": {'1':os.path.join(root,"s_5_100k.txt")}, 'control': True},
                           2: {'name': "TestB", "runs": {'1':os.path.join(root,"s_6_100k.txt")}, 'control': False}
                          },
                "options": {}
               }
        output = rnaseq.rnaseq_workflow(data, job_or_dict="dict", via="lsf")

main()
