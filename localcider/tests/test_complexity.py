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
import numpy as np
import os.path
import time
import unittest
import random
from localcider.sequenceParameters import SequenceParameters
import testTools


class TestComplexityFunctions(unittest.TestCase):
    


    def setUp(self):
        self.SP_60 = SequenceParameters('QWERTYIPASDFGHKLCVNMQWERTYIPASDFGHKLCVNMQWERTYIPASDFGHKLCVNM')
        self.SP_10 = SequenceParameters('KDNIKHVPGG')

            
    def test_get_reduced_alphabet_sequence_predefined_alphabets(self):
        for i in [2,3,4,5,6,8,10,11, 12,15,18,20]:
            
            # checks that a sequence with all 20 amino acids returns a reduced-alphabet sequence 
            # made up of i-residues which contains exactly i residues!
            self.assertEqual(len(set(self.SP_60.get_reduced_alphabet_sequence(i)[0])),i)

            random_seqs = testTools.generate_random_sequence(10)

            for j in random_seqs:

                # build obj
                SP = SequenceParameters(j)

                # check reduced alphabet sequence length matches
                self.assertEqual(len(SP.get_reduced_alphabet_sequence(i)[0]), len(j))


    def test_get_linear_WF_complexity(self):
        
        # general test
        random_seqs = testTools.generate_random_sequence_list(10, minLen=15, maxLen=500)        
        for j in random_seqs:
            SequenceParameters(j).get_linear_complexity('WF',alphabetSize=5)
                    
        # test all alphabets
        for i in [2,3,4,5,6,8,10,11, 12,15,18,20]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('WF',alphabetSize=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('WF',blobLen=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('WF',blobLen=i)

        # test a range of step-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('WF',stepSize=i)

        # test a range of word-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('WF',wordSize=i)

    def test_get_linear_LC_complexity(self):
        
        # general test
        random_seqs = testTools.generate_random_sequence_list(10, minLen=15, maxLen=500)        
        for j in random_seqs:
            SequenceParameters(j).get_linear_complexity('LC',alphabetSize=5)
                    
        # test all alphabets
        for i in [2,3,4,5,6,8,10,11, 12,15,18,20]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LC',alphabetSize=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LC',blobLen=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LC',blobLen=i)

        # test a range of step-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LC',stepSize=i)

        # test a range of word-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LC',wordSize=i)

    def test_get_linear_LZW_complexity(self):
        
        # general test
        random_seqs = testTools.generate_random_sequence_list(10, minLen=15, maxLen=500)        
        for j in random_seqs:
            SequenceParameters(j).get_linear_complexity('LZW',alphabetSize=5)
                    
        # test all alphabets
        for i in [2,3,4,5,6,8,10,11, 12,15,18,20]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LZW',alphabetSize=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LZW',blobLen=i)

        # test a range of window sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LZW',blobLen=i)

        # test a range of step-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LZW',stepSize=i)

        # test a range of word-sizes
        for i in [1,2,3,4,5,10]:
            for j in random_seqs:
                SequenceParameters(j).get_linear_complexity('LZW',wordSize=i)


    def test_complexity_values(self):

        out = {}

        # WF alphatbetsize 4
        out['WF'] = np.array([ 0.92321967,  0.92321967,  0.92321967,  0.88048202,  0.88048202,
                               0.92321967,  0.92321967,  0.92321967,  0.94773092,  0.92321967,
                               0.86096405,  0.86096405,  0.86096405,  0.68048202,  0.68048202,
                               0.86096405,  0.92321967,  0.92321967,  0.96096405,  0.9854753 ,
                               0.92321967,  0.92321967,  0.92321967,  0.88048202,  0.88048202,
                               0.92321967,  0.92321967,  0.92321967,  0.94773092,  0.92321967,
                               0.86096405,  0.86096405,  0.86096405,  0.68048202,  0.68048202,
                               0.86096405,  0.92321967,  0.92321967,  0.96096405,  0.9854753 ,
                               0.92321967,  0.92321967,  0.92321967,  0.88048202,  0.88048202,
                               0.92321967,  0.92321967,  0.92321967,  0.94773092,  0.92321967,
                               0.86096405])
        
        # LZW alphabate size=4
        out['LZW'] = np.array([ 0.9,  0.9,  0.9,  0.9,  0.9,  0.9,  0.9,  0.9,  0.8,  0.8,  0.8,
                                0.8,  0.8,  0.8,  0.8,  0.9,  0.9,  0.9,  1. ,  1. ,  0.9,  0.9,
                                0.9,  0.9,  0.9,  0.9,  0.9,  0.9,  0.8,  0.8,  0.8,  0.8,  0.8,
                                0.8,  0.8,  0.9,  0.9,  0.9,  1. ,  1. ,  0.9,  0.9,  0.9,  0.9,
                                0.9,  0.9,  0.9,  0.9,  0.8,  0.8,  0.8])

        # LC alphabet size = 2
        out['LC'] = np.array([ 0.625,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.625,  0.75 ,
                            0.625,  0.75 ,  0.625,  0.75 ,  0.875,  0.75 ,  0.625,  0.625,
                            0.625,  0.625,  0.75 ,  0.75 ,  0.625,  0.5  ,  0.5  ,  0.5  ,
                            0.5  ,  0.5  ,  0.625,  0.75 ,  0.625,  0.75 ,  0.625,  0.75 ,
                            0.875,  0.75 ,  0.625,  0.625,  0.625,  0.625,  0.75 ,  0.75 ,
                            0.625,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.5  ,  0.625,  0.75 ,
                            0.625,  0.75 ,  0.625])
        
        self.assertEquals((self.SP_60.get_linear_complexity('WF',  alphabetSize=4)[1]  - out['WF'] < 0.00001).all(), True)
        self.assertEquals((self.SP_60.get_linear_complexity('LZW', alphabetSize=4)[1] - out['LZW'] < 0.00001).all(), True)
        self.assertEquals((self.SP_60.get_linear_complexity('LC',  alphabetSize=2)[1]  - out['LC'] < 0.00001).all(), True)


                                                              

