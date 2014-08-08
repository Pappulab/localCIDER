""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.0                                                         !
   !                                                                          !
   !    Copyright (C) 2014, The localCIDER development team (current and      !
   !                        former contributors): Alex Holehouse, James       !
   !                        Ahad, Rahul K. Das.                               !
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


   Unit tests for SequenceParameter functions


"""



import os.path, time
import unittest
from localcider import sequenceParameters

class TestSequenceParametersFunctions(unittest.TestCase):

    def setUp(self):

        
        self.testObj = sequenceParameters.SequenceParameters("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA")



    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    # SECTION 1: Testing SequenceParameter object functions
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    def test_clear_phosphosites(self):
        # set phosphosites
        #self.testObj.
        pass

    def test_get_sequence(self):
        seq=self.testObj.get_sequence()
        self.assertEqual(seq, "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA")
        
    def test_get_mean_hydropathy(self):
        MH = self.testObj.get_mean_hydropathy()
        self.assertEqual(MH,3.242142857142857)

    def test_get_uversky_hydrophobicity(self):
        MUH = self.testObj.get_uversky_hydrophobicity()
        self.assertEqual(MUH, 0.4552380952380952)
        
    def test_get_fraction_disorder_promoting(self):
        FD = self.testObj.get_fraction_disorder_promoting()
        self.assertEqual(FD, 0.7285714285714285)
        
    def test_get_amino_acid_fractions(self):
        AADict=self.testObj.get_amino_acid_fractions()

        PC={'A': 0.1357142857142857,
            'C': 0.0,
            'D': 0.04285714285714286,
            'E': 0.12857142857142856,
            'F': 0.014285714285714285,
            'G': 0.12857142857142856,
            'H': 0.007142857142857143,
            'I': 0.014285714285714285,
            'K': 0.10714285714285714,
            'L': 0.02857142857142857,
            'M': 0.02857142857142857,
            'N': 0.02142857142857143,
            'P': 0.03571428571428571,
            'Q': 0.04285714285714286,
            'R': 0.0,
            'S': 0.02857142857142857,
            'T': 0.07142857142857142,
            'V': 0.1357142857142857,
            'W': 0.0,
            'Y': 0.02857142857142857}

        for AA in AADict:
            self.assertEqual(AADict[AA], PC[AA])
            
        
    def test_get_kappa(self):
        KP=self.testObj.get_kappa()        
        self.assertEqual(KP, 0.17167490749129125)
       
    def test_get_countPos(self):
        CP = self.testObj.get_countPos()
        self.assertEqual(CP,15)

    def test_get_countNeg(self):
        CN = self.testObj.get_countNeg()
        self.assertEqual(CN, 24)

    def test_get_countNeut(self):
        CNeu = self.testObj.get_countNeut()
        self.assertEqual(CNeu, 101)

    def test_get_fraction_positive(self):
        FPos = self.testObj.get_fraction_positive()
        self.assertEqual(FPos, 0.10714285714285714)
        
    def test_get_fraction_negative(self):
        FNeg = self.testObj.get_fraction_negative()
        self.assertEqual(FNeg, 0.17142857142857143)

    def test_get_FCR(self):
        FCR = self.testObj.get_FCR()
        self.assertEqual(FCR, 0.2785714285714286)
        
    def test_get_NCPR(self):
        NCPR = self.testObj.get_NCPR()
        self.assertEqual(NCPR, -0.06428571428571428)
        
    def test_get_mean_net_charge(self):
        MNC = self.testObj.get_mean_net_charge()
        self.assertEqual(MNC, 0.06428571428571428)
        
    def test_get_phasePlotRegion(self):
        PPR = self.testObj.get_phasePlotRegion()
        self.assertEqual(PPR, 2)
        
    def test_get_phosphosites(self):
        PS_pre = self.testObj.get_phosphosites()
        self.assertEqual(len(PS_pre),0)

        
    def test_get_kappa_after_phosphorylation(self):
        
        # check we're starting with the phosphosites set empty
        self.assertEqual(len(self.testObj.get_phosphosites()),0)

        # run and test phosphorylation
        self.testObj.set_phosphosites([33,39,42])
        self.assertEqual(len(self.testObj.get_phosphosites()),3)

        PostPhosK = self.testObj.get_kappa_after_phosphorylation()
        self.assertEqual(PostPhosK, 0.15101360943849682)

        pass
    def test_get_all_phosphorylatable_sites(self):
        
        all_sites = self.testObj.get_all_phosphorylatable_sites()
        
        self.assertEqual(all_sites, [9, 22, 33, 39, 42, 44, 54, 59, 64, 72, 75, 81, 87, 92, 125, 129, 133, 136])
    
    def test_get_full_phosphostatus_kappa_distribution(self):
        
        # just to be sure
        self.testObj.clear_phosphosites()
        self.testObj.set_phosphosites([33,39,42])

        FPKD = [(0.17167490749129125,
          0.10714285714285714,
          0.17142857142857143,
          0.2785714285714286,
          -0.06428571428571428,
          3.242142857142857),
         (0.16068627291329185,
          0.10714285714285714,
          0.17857142857142858,
          0.2857142857142857,
          -0.07142857142857142,
          3.222857142857143),
         (0.16863038053131135,
          0.10714285714285714,
          0.17857142857142858,
          0.2857142857142857,
          -0.07142857142857142,
          3.226428571428571),
         (0.1609064850879989,
          0.10714285714285714,
          0.18571428571428572,
          0.29285714285714287,
          -0.07857142857142857,
          3.207142857142857),
         (0.161014373883184,
          0.10714285714285714,
          0.17857142857142858,
          0.2857142857142857,
          -0.07142857142857142,
          3.2221428571428574),
         (0.15064319064948062,
          0.10714285714285714,
          0.18571428571428572,
          0.29285714285714287,
          -0.07857142857142857,
          3.2028571428571424),
         (0.15825014591368247,
          0.10714285714285714,
          0.18571428571428572,
          0.29285714285714287,
          -0.07857142857142857,
          3.2064285714285714),
         (0.15101360943849682,
          0.10714285714285714,
          0.19285714285714287,
          0.3,
          -0.08571428571428572,
          3.1871428571428573)]

        FPKD_calc = self.testObj.get_full_phosphostatus_kappa_distribution()

        self.assertEqual(FPKD_calc, FPKD)
        

    def test_set_phosphosites(self):
        
        # note we also check that incorrect phosphosites are not
        # assigned
                
        self.testObj.clear_phosphosites()

        print ""
        print "###### WE EXPECT A WARNING ON THE LINE BELOW ######"
        self.testObj.set_phosphosites([33,39,42,43])
        print "###################################################"
        self.assertEqual(self.testObj.get_phosphosites(),[33,39,42])
    

    def test_save_phaseDiagramPlot(self):
        # Should test the creation date but that's for a later update...
        # note that if we can create and save the file we can also display it
        self.testObj.save_phaseDiagramPlot("tmpfiles/phaseplot_test")
        

    def save_uverskyPlot(self):
        # 
        # 
        self.testObj.save_uverskyPlot("tmpfiles/uversky_test")

