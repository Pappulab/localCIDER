""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.0                                                         !
   !                                                                          !
   !    Copyright (C) 2014, The LocalCIDER development team (current and      !
   !                        former contributors): Alex Holehouse, James       !
   !                        Ahad, Rahul K. Das.                               !
   !                                                                          !
   !    localCIDER was developed in the lab of Rohit Pappu at Washington      !
   !    University in St. Louis. Please see the website for citation          !
   !    information:                                                          !
   !                                                                          !
   !    http://pappulab.github.io/LocalCIDER/                                 !
   !                                                                          !
   !    For more information please see the Pappu lab website:                !
   !                                                                          !
   !    http://pappulab.wustl.edu/                                            !
   !                                                                          !
   !    LocalCIDER is free software: you can redistribute it and/or modify    !
   !    it under the terms of the GNU General Public License as published by  !
   !    the Free Software Foundation, either version 3 of the License, or     !
   !    (at your option) any later version.                                   !
   !                                                                          !
   !    LocalCIDER is distributed in the hope that it will be useful,         !
   !    but WITHOUT ANY WARRANTY; without even the implied warranty of        !
   !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
   !    GNU General Public License for more details.                          !
   !                                                                          !
   !    You should have received a copy of the GNU General Public License     !
   !    along with LocalCIDER.  If not, see <http://www.gnu.org/licenses/>.   !
   !--------------------------------------------------------------------------!
   ! AUTHORSHIP INFO:                                                         !
   !--------------------------------------------------------------------------!
   !                                                                          !
   ! MAIN AUTHOR:   Alex Holehouse                                            !
   !                                                                          !
   !--------------------------------------------------------------------------!

   
   File Description:
   ================
   
   This file contains functionality for generating sequence permutants via
   the Wang-Landau algrorithm. There are two seperate operations which can
   be carried out here;


   1) Estimate the Kappa density of states over some region

   2) Create permutants of your sequence with various kappa values

  
   The general approach taken is that you create a sequencePermutant
   object with some sequence, and then call an operation on that object.


"""

class SequencePermutantException(Exception):
    pass


import numpy as np

from backend.sequence import Sequence
from backend.seqfileparser import SequenceFileParser 
from backend.backendtools import status_message

from backend import wang_landau

class SequencePermutants:

    def __init__(self, sequence="", sequenceFile=""):
        """
        sequencePermutant objects, like sequenceParameter objects, can be initialized by either an amino acid sequence as a string, or a sequence file. Such sequence files can be FASTA files, or simply an amino acid sequence in a text file.

        """

        # provide the flexibility to submit either a sequence
        # file or an actual sequence as a string
        if sequence=="" and sequenceFile=="":
            return None

        # if the sequence isn't empty constuct a local 
        # sequence object using the sequence
        if not sequence == "":
            self.SeqObj = Sequence(sequence)
        else:
            parserMachine = SequenceFileParser()
            self.SeqObj = Sequence(parserMachine.parseSeqFile(sequenceFile))
            
        self.WL_ready=False
        

    def initializeWangLandauParameters(self, OUTDIR, frozen=set([]),  nbins=10, binmin=0, binmax=1, flatchck=10000, flatcrit=0.7, convergence=np.exp(.000001), WL_type="NORMAL"):
        self.WLM = wang_landau.WangLandauMachine(self.SeqObj.seq,  OUTDIR, frozen, nbins, binmin,  binmax,  flatchck, flatcrit, convergence, WL_type)

        self.WL_ready=True
                                                 

    def run_wang_landau_DOS(self, binmin, minmax, numberofBins):
        self.__readyCheck()


    def run_histogramZoom_DOS_estimation(self, nbins):
        pass

    
           
    def run_generatePermutants(self, numSequencesPerBin):
        pass
        
    def __readyCheck(self):
        if self.WL_ready == False:
            raise SequencePermutantException("Wang Landau not yet parameterized\nPLEASE RUN <initializeWangLandauParameters> BEFORE CALLING a run function")
        
