#!/usr/bin/env python
from test_zsi import *

def makeTestSuite():
    return unittest.TestSuite(
        map(lambda t: globals()[t].makeTestSuite(),
            filter(lambda g: g.startswith('test_') and True, globals()))
    )

def main():
    unittest.main(defaultTest="makeTestSuite")
    suite = unittest.TestSuite()

if __name__ == "__main__" : main()
