"""
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.9                                                         !
   !--------------------------------------------------------------------------!

   File Description:
   ================

   ***********************************************************************
   *                                                                     *
   * WARNING: SequencePermutants is still undergoing testing and will be *
   * released in version 0.2.0                                           *
   *                                                                     *
   ************************************************************************

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

from sequenceParameters import SequenceParameters
from backend.sequence import Sequence
from backend.seqfileparser import SequenceFileParser
from backend.backendtools import status_message, warn_notReadyYet

from backend import wang_landau


class SequencePermutants:

    #...................................................................................#
    def __init__(self, sequence="", sequenceFile=""):
        """
        sequencePermutant objects, like sequenceParameter objects, can be initialized by either an amino acid 
        sequence as a string, or a sequence file. Such sequence files can be FASTA files, or simply an amino 
        acid sequence in a text file.

        """

        warn_notReadyYet()

        # provide the flexibility to submit either a sequence
        # file or an actual sequence as a string
        if sequence == "" and sequenceFile == "":
            return None

        # if the sequence isn't empty constuct a local
        # sequence object using the sequence
        if not sequence == "":
            self.SeqObj = Sequence(sequence)
        else:
            parserMachine = SequenceFileParser()
            self.SeqObj = Sequence(parserMachine.parseSeqFile(sequenceFile))

        self.WL_ready = False

    #...................................................................................#
    def initializeWangLandauParameters(
            self,
            OUTDIR,
            frozen=set(
                []),
            nbins=10,
            binmin=0,
            binmax=1,
            flatchck=10000,
            flatcrit=0.7,
            convergence=np.exp(.000001),
            WL_type="NORMAL"):

        # initialize the Wang-Landau Machine
        self.WLM = wang_landau.WangLandauMachine(
            self.SeqObj.seq,
            OUTDIR,
            frozen,
            nbins,
            binmin,
            binmax,
            flatchck,
            flatcrit,
            convergence,
            WL_type)

        # set the ready flag
        self.WL_ready = True

    #...................................................................................#
    def run_wang_landau_DOS(self, binmin, minmax, numberofBins):
        # self.__readyCheck()
        warn_notReadyYet()

    #...................................................................................#
    def run_histogramZoom_DOS_estimation(self, nbins):
        warn_notReadyYet()

    #...................................................................................#
    def run_generatePermutants(self, numSequencesPerBin):
        warn_notReadyYet()

    #...................................................................................#
    def __readyCheck(self):
        if not self.WL_ready:
            raise SequencePermutantException(
                "Wang Landau not yet parameterized\nPLEASE RUN <initializeWangLandauParameters> BEFORE CALLING a run function")

    #...................................................................................#
    def get_permutant(self):

        # create a permutant Sequence object
        SO = self.SeqObj.full_shuffle([])

        # now create an 'empty' SequenceParameter object
        SP = SequenceParameters("A")

        # manually set the underlying Sequence object
        # and return
        SP.SeqObj = SO
        return SP


