#!/usr/bin/python

import unittest
import __init__ as test

ALL=True
sequenceParameters_test=True
plots_test=True

if plots_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.test_plots.TestPlotsFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)


if sequenceParameters_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(test.test_sequenceParameters.TestSequenceParametersFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)

    
