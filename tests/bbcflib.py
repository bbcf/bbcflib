import unittest

def suites():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite("bbcflib.email")
    return suite

if __name__ == "__main__":
    suites()
