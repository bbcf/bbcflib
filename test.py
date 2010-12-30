from bbcflib import *
from unittest import TestCase, TestSuite, main

class TestGenRep(TestCase):
    pass

class TestEmailReport(TestCase):
    pass

class TestDAFLIMS(TestCase):
    pass

class TestFrontend(TestCase):
    pass

def test_all():
    main()

if __name__ == '__main__':
    test_all()
