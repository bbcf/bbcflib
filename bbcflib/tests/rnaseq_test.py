import re
import sys
import os
import random
import rnaseq
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
        output = rnaseq.rnaseq_workflow(data, job_or_dict="dict", via="lsf", maplot="normal")

main()
