"""
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.9                                                         !
   !                                                                          !
   !    Copyright (C) 2014 - 2016                                             !
   !    The localCIDER development team (current and former contributors)     !
   !    Alex Holehouse, James Ahad, Rahul K. Das.                             !
   !                                                                          !
   !    localCIDER was developed in the lab of Rohit Pappu at Washington      !
   !    University in St. Louis. Please see the website for citation          !
   !    information:                                                          !
   !                                                                          !
   !    http://pappulab.github.io/localCIDER/                                 !
   !                                                                          !
   !    For more information please see the Pappu lab website:                !
   !                                                                          !
   !    http://pappulab.wustl.edu/                                            !
   !                                                                          !
   !    localCIDER is free software: you can redistribute it and/or modify    !
   !    it under the terms of the GNU General Public License as published by  !
   !    the Free Software Foundation, either version 3 of the License, or     !
   !    (at your option) any later version.                                   !
   !                                                                          !
   !    localCIDER is distributed in the hope that it will be useful,         !
   !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
   !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
   !    GNU General Public License for more details.                          !
   !                                                                          !
   !    You should have received a copy of the GNU General Public License     !
   !    along with localCIDER.  If not, see <http://www.gnu.org/licenses/>.   !
   !--------------------------------------------------------------------------!
   ! AUTHORSHIP INFO:                                                         !
   !--------------------------------------------------------------------------!
   !                                                                          !
   ! MAIN AUTHOR:   Alex Holehouse                                            !
   !                                                                          !
   !--------------------------------------------------------------------------!


   Runner code for localcider tests. The various flags define which tests are run
   setting ALL to True means that all the tests will be run regardless of the 
   status of the other flags.

   This is set up to make it easy to quickly test specific components, rather
   than running the entire testing suite at once.

"""

import unittest
import __init__ as test


## TEST flags
ALL = True
sequenceParameters_test = True
plots_test = False
sequence_test = False
complexity_test = False

## plotting tests
if plots_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(
        test.test_plots.TestPlotsFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)

## sequence parameter tests
if sequenceParameters_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(
        test.test_sequenceParameters.TestSequenceParametersFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)

## sequence tests
if sequence_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(
        test.test_sequence.TestSequenceFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)

## complexity tests
if complexity_test or ALL:
    suite = unittest.TestLoader().loadTestsFromTestCase(
        test.test_complexity.TestComplexityFunctions)

    unittest.TextTestRunner(verbosity=2).run(suite)


