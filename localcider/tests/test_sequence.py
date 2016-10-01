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


   Unit tests for functions from backened/sequence.py


"""

import os.path
import time
import unittest
import random
from localcider.backend import sequence
import testTools


def roundflt(num):
    return float("%1.2f" % num)

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        self.SO_60 = sequence.Sequence('MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP')
        self.SO_10 = sequence.Sequence('MEEPQSDPSV')


    def test_FCR(self):
        self.assertEqual(self.SO_60.FCR(), 0.25)
        self.assertEqual(self.SO_10.FCR(), 0.3)


    def test_NCPR(self):
        self.assertEqual(roundflt(self.SO_60.NCPR()), -0.22)
        self.assertEqual(roundflt(self.SO_10.NCPR()), -0.3)


    def test_fplus_fminus(self):
        self.assertEqual(roundflt(self.SO_60.Fminus()), 0.23)
        self.assertEqual(roundflt(self.SO_60.Fplus()), 0.02)

        self.assertEqual(roundflt(self.SO_10.Fminus()), 0.3)
        self.assertEqual(roundflt(self.SO_10.Fplus()),  0.0)

    def test_aminoAcidColorMap(self):
        
        # just code coverage
        CM = self.SO_10.aminoAcidColorMap
        self.assertEqual(len(CM), 20)

    def test_amino_acid_fraction(self):
        
        vals = {'A': 0.016666666666666666,
                'C': 0.0,
                'D': 0.11666666666666667,
                'E': 0.11666666666666667,
                'F': 0.03333333333333333,
                'G': 0.016666666666666666,
                'H': 0.0,
                'I': 0.016666666666666666,
                'K': 0.016666666666666666,
                'L': 0.13333333333333333,
                'M': 0.05,
                'N': 0.03333333333333333,
                'P': 0.16666666666666666,
                'Q': 0.06666666666666667,
                'R': 0.0,
                'S': 0.11666666666666667,
                'T': 0.03333333333333333,
                'V': 0.03333333333333333,
                'W': 0.03333333333333333,
                'Y': 0.0}

        fractions = self.SO_60.amino_acid_fraction()
        for i in vals:
            self.assertEqual(roundflt(fractions[i]), roundflt(vals[i]))


    def test_linear_vectors(self):
        
        # charge pattern vector
        self.assertEqual(len(self.SO_60.chargePattern), 60)

        self.assertEqual(len(self.SO_60.linearDistOfFCR(5)[0]), 60)
        self.assertEqual(len(self.SO_60.linearDistOfFCR(6)[0]), 60)

        self.assertEqual(len(self.SO_60.linearDistOfNCPR(5)[0]), 60)
        self.assertEqual(len(self.SO_60.linearDistOfNCPR(6)[0]), 60)

        self.assertEqual(len(self.SO_60.linearDistOfHydropathy(5)[0]), 60)
        self.assertEqual(len(self.SO_60.linearDistOfHydropathy(6)[0]), 60)

        self.assertEqual(len(self.SO_60.linearDistOfHydropathy_2(5)[0]), 60)
        self.assertEqual(len(self.SO_60.linearDistOfHydropathy_2(6)[0]), 60)

        self.assertEqual(len(self.SO_60.linearDistOfSigma(5)[0]), 60)
        self.assertEqual(len(self.SO_60.linearDistOfSigma(6)[0]), 60)
            
        



