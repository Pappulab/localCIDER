#!/usr/bin/python

import unittest
import __init__ as test

ALL=True
SP_test=True

if SP_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.test_sequenceParameters.TestSequenceParametersFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
