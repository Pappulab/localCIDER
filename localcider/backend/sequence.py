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
   ! MAIN AUTHOR:   James Ahad and Alex Holehouse                             !
   !                                                                          !
   !--------------------------------------------------------------------------!


   File Description:
   ================

   This is the main file which deals with sequence objects, and is where all
   sequence analysis operations occur. These functions should not be called
   by normal users - if there is functionality or logical carried out in here
   you would like to see in localCIDER please submit a feature request and we
   can integrate these changes into the sequenceparameter API.


REFERENCES
R.K. Das, R.V. Pappu. (2013).
"Conformations of intrinsically disordered proteins are  influenced by
linear sequence distributions of oppositely charged residues."

Proceedings of the National Academy of Sciences USA. 110: 13392-13397

"""

import numpy as np
import random as rng
import time
import copy as cp
import os
import itertools
from backendtools import return_absolute_datafile_path, warning_message, verifyType, status_message, warning_message
from restable import ResTable
from data import aminoacids
import data
from localciderExceptions import SequenceException
from sequenceComplexity import SequenceComplexity


######################

# lookuptable is a global table of AA properties
#
# Hydropathy of amino acids is determined by the KD table
#
lkupTab = ResTable()


class Sequence:
    """
    The sequence class is the main class which holds sequence properties
    and carries out sequence parameter logic operations
    """

    #...................................................................................#
    def __init__(self, seq, dmax=-1, chargePattern=[], validateSeq=False):
        """
        seq = amino acid sequence as a string

        """

        if not verifyType(seq, str):
            raise SequenceException(
                "Must pass a string to a new Sequence object")

        # by default don't validate
        if validateSeq:
            seq = seq.upper()
            seq = self.validateSequence(seq)

        self.seq = seq.upper()
        self.len = len(seq)
        self.chargePattern = chargePattern

        if(chargePattern == []):
            for i in np.arange(0, self.len):
                if(lkupTab.lookUpCharge(self.seq[i]) > 0):
                    chargePattern = np.append(chargePattern, 1)
                elif(lkupTab.lookUpCharge(self.seq[i]) < 0):
                    chargePattern = np.append(chargePattern, -1)
                else:
                    chargePattern = np.append(chargePattern, 0)
            self.chargePattern = chargePattern
        # initializing to prevent extra computational time
        self.dmax = dmax

        # set phosphosites as empty
        self.phosphosites = []

        # set color pallete
        self.set_HTMLColorResiduePalette(aminoacids.DEFAULT_COLOR_PALETTE)

        # build a complexity object (just code encapsulation - no state)
        self.ComplexityObject = SequenceComplexity()

    #...................................................................................#
    def __unicode__(self):
        """ Returns the sequences """
        return self.seq

    #...................................................................................#
    def __str__(self):
        """ Returns the sequences """
        return self.seq

    #...................................................................................#
    def validateSequence(self, seq):

        processed = ""

        AAs = data.aminoacids.ONE_TO_THREE.keys()
        pos = 0
        messageWarned = False

        # for each residue in your protein sequence
        for i in seq:
            pos = pos + 1
            if i not in AAs:

                # if we find whitespace
                if i.isspace():
                    if not messageWarned:
                        # only warn once...
                        status_message("Removing whitespace from sequence")
                        messageWarned = True
                    pass
                # if unexpected residue/character bail
                else:
                    raise SequenceException(
                        "Invalid amino acid [" + str(i) + "] found at position " + str(pos))

            # else append sequence to the processed sequence
            else:
                processed = processed + i

        # determine proline content and warn if over 15%
        prolineContent = float(processed.count("P")) / float(len(processed))
        if prolineContent > 0.15:
            warning_message(
                "This sequence has a proline content of greater than 15%.\nThis may render some analyses [notably kappa and phase diagram predictions] incorrect")

        return processed

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                           GENERAL SEQUENCE PROPERTIES
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #...................................................................................#
    def fraction_disorder_promoting(self):
        """
        Get the fraction of 'disorder promoting residues' in a
        sequence
        """

        D = ['T', 'A', 'G', 'R', 'D', 'H', 'Q', 'K', 'S', 'E', 'P']
        # the O list is redundant but I've included it for good measure!
        O = ['W', 'F', 'Y', 'I', 'M', 'L', 'V', 'N', 'C']

        D_count = 0
        O_count = 0
        for i in self.seq:
            if i in D:
                D_count += 1
            else:
                O_count += 1
        return float(D_count) / len(self.seq)

    def amino_acid_fraction(self):
        AADICT = {'A': 0,
                  'C': 0,
                  'D': 0,
                  'E': 0,
                  'F': 0,
                  'G': 0,
                  'H': 0,
                  'I': 0,
                  'K': 0,
                  'L': 0,
                  'M': 0,
                  'N': 0,
                  'P': 0,
                  'Q': 0,
                  'R': 0,
                  'S': 0,
                  'T': 0,
                  'V': 0,
                  'W': 0,
                  'Y': 0}
        for i in self.seq:
            AADICT[i] += 1

        for i in AADICT:
            AADICT[i] = float(AADICT[i]) / float(len(self.seq))

        return AADICT

    #...................................................................................#
    def countPos(self):
        """ Get the number of positive residues in the sequence """
        return len(np.where(self.chargePattern > 0)[0])

    #...................................................................................#
    def countNeg(self):
        """ Get the number of negative residues in the sequence """
        return len(np.where(self.chargePattern < 0)[0])

    #...................................................................................#
    def countNeut(self):
        """ Get the number of neutral residues in the sequence """
        return len(np.where(self.chargePattern == 0)[0])

    #...................................................................................#
    def Fplus(self):
        """ Get the fraction of positive residues in the sequence """
        return self.countPos() / (self.len + 0.0)

    #...................................................................................#
    def Fminus(self):
        """ Get the fraction of negative residues in the sequence """
        return self.countNeg() / (self.len + 0.0)

    #...................................................................................#
    def FCR(self):
        """ Get the fraction of charged residues in the sequence """
        return (self.countPos() + self.countNeg()) / (self.len + 0.0)

    #...................................................................................#
    def FER(self):
        """ Get the fraction of charged residues in the sequence """
        return (self.countPos() + self.countNeg() + self.seq.count('P')) / (self.len + 0.0)

    #...................................................................................#
    def NCPR(self):
        """ Get the net charge per residue of the sequence """
        return (self.countPos() - self.countNeg()) / (self.len + 0.0)

    #...................................................................................#
    def mean_net_charge(self):
        """ Get the absolute magnitude of the mean net charge """
        return abs(self.NCPR())

    #...................................................................................#
    def kappa(self):
        """
        Return the kappa value, as defined in REF 1 \
        """

        if self.deltaMax() == 0:
            warning_message(
                "The sequence has no charged residues - kappa is not a valid/relevant parameter")
            return -1
        else:

            kappaVal = self.delta() / self.deltaMax()

            # so the heuristics for kappa are good BUT may under estimate
            # deltaMax is some cases. If this is a small deviation then we
            # just set it to 0 because the sequence with the highest delta is probably
            # an sequence-Isomer of the sequence we have. If this deviation is larger,
            # however, it may be indicative of a bug in the code which we
            # should address 
            #
            # It's probably worth pointing out that after > 10 Million sequences 
            # analyzed we're yet to find a bug - May 2016
            if kappaVal > 1.0 and kappaVal < 1.1:
                return 1.0
            else:
                return kappaVal


    #...................................................................................#
    def kappa_proline(self):
        """
        Return the kappa-proline value, as defined in Martin et al. 

        """

        # firstly convert the sequence into a 2 letter alphabet of
        # E/D/R/K/P or 'other'
        newseq=''
        for res in self.seq:
            if res == 'P' or res =='E' or res =='D' or res =='K' or res =='R': 
                newseq=newseq+'E'
            else:
                newseq=newseq+'K'

        # then calculate the new kappa value in this new sequence space and
        # return the value
        augmented_seq = Sequence(newseq)

        return augmented_seq.kappa()


    #...................................................................................#
    def kappa_proline_seq(self):
        """
        Return the kappa-proline sequence, where R/K/D/E/P are X and all other residues
        are O

        """

        # firstly convert the sequence into a 2 letter alphabet of
        # E/D/R/K/P or 'other'
        newseq=''
        for res in self.seq:
            if res == 'P' or res =='E' or res =='D' or res =='K' or res =='R': 
                newseq=newseq+'X'
            else:
                newseq=newseq+'O'

        return newseq


    #...................................................................................#
    def kappa_X(self, grp1, grp2=None):
        """
        Return the user defined kappa value, where residues are subdivided into one
        or two arbitrary groups. If only grp1 is defined then the 20 amino acids are
        divided into grp1 and not grp1, and the patterning is determined for the perfect
        binary system. If grp2 is defined the 20 amino acids are divided into grp1, grp2,
        not grp1 or grp2, and the binary patterning of grp1 and grp2 residues is determined
        in the context of a ternary system.

        As an example, to get the original kappa definition, grouping would be
        grp1 = ['E','D']
        grp2 = ['K','D']

        While for the kappa-proline definition the grouping would be
        grp1 = ['P','E','D','K','R']

        """                
        
        grp1 = self.__parse_group(grp1)

        if grp2:
            grp2 = self.__parse_group(grp2)

        if grp2:
            # firstly convert the sequence into a 3 letter alphabet of
            # grp1, grp2, or 'other'
            newseq=''
            for res in self.seq:
                if res in grp1:
                    newseq=newseq+'E'                
                elif res in grp2:
                    newseq=newseq+'K'
                else:
                    newseq=newseq+'G'

        else:
            # firstly convert the sequence into a 2 letter alphabet of
            # grp1 residue or 'other'
            newseq=''
            for res in self.seq:
                if res in grp1:
                    newseq=newseq+'E'                
                else:
                    newseq=newseq+'K'
            
        # then calculate the new kappa value in this new sequence space and
        # return the value
        augmented_seq = Sequence(newseq)

        return augmented_seq.kappa()


    #...................................................................................#
    def phasePlotRegion(self):
        """ Returns the IDP diagram of states [REF 1] region based on the FCR/NCPR
            rules. Possible return values are;


            1 +-----------------------+---------------+
              |                     +                 |
              |                   +                   |
              |      4          +                     |
        F     |               +                       |
        R     |             +                         |
        A     |           +                           +
        C     |         +                           + |
        T     |       +                           +   |
        .     |     +            3              +     |
          0.5 |   +                           +       |
        N     | +                           +         |
        E     +                           +           |
        G     | +                       +             |
              + 2 +                   +               |
        R     | +   +               +       4         |
        E     |   + 2 +           +                   |
        S     |     +   +       +                     |
              | 1     + 2 +   +                       |
           0  +---------+---+-------------------------+

              0                  0.5                  1
                    FRACTION OF POSITIVE RESIDUES


        1) Weak polyampholytes and polyelectrolyte
        2) Janus sequences
        3) Strong polyampholytes
        4) Strong polyelectrolytes

        """
        fcr = self.FCR()
        ncpr = self.NCPR()

        # if we're in region 1
        if(fcr < .25):
            return 1

        # if we're in region 2
        elif(fcr >= .25 and fcr <= .35):
            return 2

        # if we're in region 3
        elif(fcr > .35 and abs(ncpr) < 0.35):
            return 3

        # if we're in region 4 or 5
        elif(self.Fplus() > 0.35):
            if self.Fminus() > 0.35:
                raise SequenceException(
                    "Algorithm bug when coping with phase plot regions")
            return 5

        elif(self.Fminus() > 0.35):
            return 4

        else:  # This case is impossible but here for completeness\
            raise SequenceException(
                "Found inaccessible region of phase diagram. Numerical error")

    #...................................................................................#
    def phasePlotAnnotation(self):
        """ Return the string annotation for the diagram of states region
            which this sequence falls into
        """
        region = self.phasePlotRegion()
        if(region == 1):
            return 'Globule/Tadpole'
        elif(region == 2):
            return 'Boundary Region'
        elif(region == 3):
            return 'Coils,Hairpins and Chimeras'
        elif(region == 4):
            return 'Negatively Charged Swollen Coils'
        elif(region == 5):
            return 'Positively Charged Swollen Coils'
        else:
            return 'ERROR, NOT A REAL REGION'

    #...................................................................................#
    def meanHydropathy(self):
        """ Return the mean hydropathy value for the sequence
        """
        ans = 0
        for i in xrange(0, self.len):
            ans += lkupTab.lookUpHydropathy(self.seq[i]) / self.len
        return ans

    #...................................................................................#
    def uverskyHydropathy(self):
        """
        Return the mean hydropathy as calculated by Uversky (window size of
        five, with normalized Kyte-Doolitle scale.

        This means we can use different hydropathy scales if we choose, but
        it won't break the Uversky plots
        """

        normalizedKD = data.aminoacids.get_KD_uversky()
        translate = data.aminoacids.ONE_TO_THREE

        ans = 0
        for idx in xrange(0, self.len):
            ans += normalizedKD[translate[self.seq[idx]]] / self.len

        return ans

    #...................................................................................#
    def cumMeanHydropathy(self):
        """ Return a vector of the average cumulative hydropathy. i.e.

        [ vector of cumulative hydropathy ] / [ vector of sequence position ]

        This allows you to examine how hydropathy changes through the sequence
        """
        ans = [lkupTab.lookUpHydropathy(self.seq[0])]
        for i in xrange(1, self.len):
            ans.append(ans[i - 1] + lkupTab.lookUpHydropathy(self.seq[i]))
        ans /= (np.arange(0, self.len) + 1)
        return ans
        
    #...................................................................................#
    def FPPII_chain(self):
        """ 
        Returns the overal chain's average PPII propensity as defined by Elam et al[1].

        [1] - Elam WA, Schrank TP, Campagnolo AJ, Hilser VJ. Evolutionary 
        conservation of the polyproline II conformation surrounding intrinsically 
        disordered phosphorylation sites. 
        Protein Sci. 2013; 22: 405- 417. doi: 10.1002/pro.2217 PMID: 23341186
        
        """
        
        # calculate the total PPII sum for the sequence
        total = 0
        for i in xrange(0, self.len):
            total = total + lkupTab.lookUpPPII(self.seq[i])

        # normalize by the sequence length
        return (total/float(self.len))
        
    #...................................................................................#
    def linearDistOfNCPR(self, bloblen):
        """
        Returns a np vertical stack object showing the NCPR over
        blob-sized regions along the sequence
        """

        # detemrine the number of blobs of length blpblen in the sequence
        nblobs = self.len - bloblen + 1

        # determine the flanking positions over which we don't calculate 
        # NCPR values BUT still want to include to build a profile which is
        # the correct length. 
        flank = int(bloblen/2)
        
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank
        else:
            flank_start = flank - 1
            flank_end   = flank 
            
        # initialize an empty list of blobs
        blobncpr = [0] * nblobs

        # for each overlapping blob in the sequence calculate the NCPR
        for i in np.arange(0, nblobs):
            blob = self.chargePattern[i:(i + bloblen)]
            bpos = len(np.where(blob > 0)[0])
            bneg = len(np.where(blob < 0)[0])
            blobncpr[i] = (bpos - bneg) / (bloblen + 0.0)

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blobncpr + [0]*flank_end))

    

    #...................................................................................#
    def linearDistOfFCR(self, bloblen):
        """
        Returns an np vertical stack object showing the FCR over
        blob-sized regions along the sequence
        """

        self.__check_window_to_length(bloblen)

        nblobs = self.len - bloblen + 1

        # determine the flanking positions over which we don't calculate 
        # FCR values BUT still want to include to build a profile which is
        # the correct length. 
        flank = int(bloblen/2)

        # if bloblen is odd
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank 
        else:
            flank_start = flank - 1
            flank_end   = flank 

        blobfcr = [0] * nblobs
        
        # for each overlapping blob in the sequence calculate the FCR
        for i in np.arange(0, nblobs):
            blob = self.chargePattern[i:(i + bloblen)]
            bpos = len(np.where(blob > 0)[0])
            bneg = len(np.where(blob < 0)[0])
            blobfcr[i] = (bpos + bneg) / (bloblen + 0.0)

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blobfcr + [0]*flank_end))

    #...................................................................................#
    def linearDistOfSigma(self, bloblen):
        """
        Returns an np vertical stackobject showing how the "sigma" parameter
        varies over blob-sized regions along the sequence
        """

        self.__check_window_to_length(bloblen)

        nblobs = self.len - bloblen + 1

        # determine the flanking positions over which we don't calculate 
        # sigma values BUT still want to include to build a profile which is
        # the correct length. 
        flank = int(bloblen/2)
        
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank
        else:
            flank_start = flank - 1
            flank_end   = flank 

        blobsig = [0] * nblobs

        # for each overlapping blob in the sequence calculate the sigma
        # parameter
        for i in np.arange(0, nblobs):
            blob = self.chargePattern[i:(i + bloblen)]
            bpos = len(np.where(blob > 0)[0])
            bneg = len(np.where(blob < 0)[0])

            bncpr = (bpos - bneg) / (bloblen + 0.0)
            bfcr = (bpos + bneg) / (bloblen + 0.0)

            if(bfcr == 0):
                bsig = 0
            else:
                bsig = bncpr**2 / bfcr
            blobsig[i] = bsig

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blobsig + [0]*flank_end))
        
    #...................................................................................#
    def linearDistOfHydropathy(self, bloblen):
        """
        Returns an np vertical stackobject showing how the 0 to 1 normalized hydropathy
        varies over blob-sized regions along the sequence
        """
        
        self.__check_window_to_length(bloblen)

        nblobs = self.len - bloblen + 1

        # determine the flanking positions over which we don't calculate 
        # hydropathy values BUT still want to include to build a profile which is
        # the correct length. 
        flank = int(bloblen/2)
        
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank
        else:
            flank_start = flank - 1
            flank_end   = flank         

        blobhydro = [0] * nblobs

        # construct a vector of the hydropathy values of each residue in
        # in the sequence (using the Uversky-normalized KD scale which runs from 0 to 1
        # with 1 being the most hydrophobic.

        # get Uversky hydropathy lookup table
        KDU = aminoacids.get_KD_uversky()

        hydrochain = []
        for i in self.seq:
            hydrochain.append(KDU[aminoacids.ONE_TO_THREE[i]])

        # for each overlapping blob in the sequence calculate the hydropathy
        for i in np.arange(0, nblobs):
            blob = hydrochain[i:(i + bloblen)]
            blobhydro[i] = sum(blob) / float(bloblen)

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blobhydro + [0]*flank_end))

    #...................................................................................#
    def linearDistOfHydropathy_2(self, bloblen):
        """
        Returns an np vertical stackobject showing how the 0 to 9 normallized hydropathy
        varies over blob-sized regions along the sequence
        """

        self.__check_window_to_length(bloblen)
        
        nblobs = self.len - bloblen + 1

        flank = int(bloblen/2)
        
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank
        else:
            flank_start = flank - 1
            flank_end   = flank

        blobhydro = [0] * nblobs

        # construct a vector of the hydropathy values of each residue in
        # in the sequence (using the re-set KD scale which runs from 0 to 9
        # with 9 being the most hydrophobic. This scaling may be changed
        # in future versions...
        hydrochain = [lkupTab.lookUpHydropathy(res) for res in self.seq]

        # for each overlapping blob in the sequence calculate the hydropathy
        for i in np.arange(0, nblobs):
            blob = hydrochain[i:(i + bloblen)]
            blobhydro[i] = sum(blob)

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blobhydro + [0]*flank_end))

    #...................................................................................#
    def linearDenistyOfAAs(self,
                                  bloblen,
                                  targetAAs):
        """
        Function which returns a np.vertical stackobject showing the 0 to 1 normalized density
        of the AAs in the targetAAs list. 


        bloblen defines the window size used, while targetAAs should be a list of one or 
        more amino acids of interest (upper case).

        """

        self.__check_window_to_length(bloblen)
        
        nblobs = self.len - bloblen + 1

        flank = int(bloblen/2)
        
        if 2*flank+nblobs == self.len:
            flank_start = flank
            flank_end   = flank
        else:
            flank_start = flank - 1
            flank_end   = flank

        blob_density = [0] * nblobs

        # construct a binary list where residues in the targetAA list
        # are set to 1 and all others are set to 0 
        target_seq=[]
        for res in self.seq:
            if res in targetAAs:
                target_seq.append(1.0)
            else:
                target_seq.append(0.0)
            
        # for each overlapping blob in the sequence calculate the local density of amino acid
        for i in np.arange(0, nblobs):
            blob = target_seq[i:(i + bloblen)]
            blob_density[i] = sum(blob)/float(bloblen)

        return np.vstack((np.arange(1, self.len + 1), [0]*flank_start + blob_density + [0]*flank_end))
    
    #...................................................................................#
    def linearCompositions(self,
                           bloblen,
                           grps=[]):

        """
        Function which returns an n by X matrix showing the density of different groups of amino
        acids running along the sequence. The amino acid group can be defined by the user by providing
        a list of lists, where each sub-list contains amino acids for a distinct group. Alternativly, the
        default grouping provides 7 groups which are

        acidic    (E and D)
        basic     (K and R)
        charged   (E, D, K, and R)
        polar     (Q, N, S, T, G, C, and H)
        aliphatic (A, L, M, I, and V)
        aromatic  (F, Y, and W)
        proline   (P)       

        """        
        
        # if we've defined our own groups (NOTE no sanity checking...)
        if len(grps) > 0:

            # then skip ahead...

            sanitized_groups = []
            for group in grps:
                sanitized_groups.append(self.__parse_group(group))

            grps = sanitized_groups

        # else use the standard AA grouping
        else:                
            grps.append(['E','D'])                        # acidic
            grps.append(['R','K'])                        # basic
            grps.append(['R','K','E','D'])                # charged
            grps.append(['Q','N','S','T','G','H', 'C'])   # polar
            grps.append(['A','L','M','I','V'])            # hydrophobic
            grps.append(['F','Y','W'])                    # aromatic
            grps.append(['P'])                            # proline
            
        # first create a np vector object with the linear density of the
        # first group
        tmp = self.linearDenistyOfAAs(bloblen, grps[0])
        density = tmp[1]
        
        # for for each remaining group just keep adding to it
        for group in grps[1:]:
            tmp = self.linearDenistyOfAAs(bloblen, group)
            density = np.vstack((density, tmp[1]))

        # finally give the         
        return (tmp[0], density)
            
            
        
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                         SEQUENCE COMPLEXITY FUNCTIONS
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #...................................................................................#
    def get_linear_WF_complexity(
            self,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1):
        """
        Returns the Wooton-Federhen (WF) vectorial complexity using a sliding window approach

        """
        
        self.__check_window_to_length(windowSize)

        return self.ComplexityObject.get_WF_complexity(
            self.seq, alphabetSize, userAlphabet, windowSize, stepSize)

    #...................................................................................#
    def get_linear_LC_complexity(
            self,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1,
            wordSize=3):
        """
        Returns the Linguistic Complexity (LC) vectorial complexity using a sliding window approach

        """

        self.__check_window_to_length(windowSize)

        return self.ComplexityObject.get_LC_complexity(
            self.seq, alphabetSize, userAlphabet, windowSize, stepSize, wordSize)

    #...................................................................................#
    def get_linear_LZW_complexity(
            self,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1):
        """
        Returns the Lempel-Ziv-Welch (LZW) vectorial complexity using a sliding window approach.

        """

        self.__check_window_to_length(windowSize)

        return self.ComplexityObject.get_LZW_complexity(
            self.seq, alphabetSize, userAlphabet, windowSize, stepSize)


    #...................................................................................#
    def get_reducedAlphabetSequence(self, alphabetSize=20, userAlphabet={}):
        """
        Returns the amino acid sequence after being translated to a reduced alphabet

        """
        return self.ComplexityObject.reduce_alphabet(
            self.seq, alphabetSize, userAlphabet)

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                         KAPPA CALCULATING FUNCTIONS
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    # The functions in the following section focus on parameters involved in the kappa
    # calculation.
    #


    #...................................................................................#
    def sigma(self):
        """ Returns the sigma value for a sequence

           \sigma = \dfrac{NCPR^2}{FCR}

           When the sequence has one or more charged residues
           sigma = (NCPR^2)/FCR

           When the sequence has no charged residues
           sigma = 0
        """
        if(self.countNeut() == self.len):
            return 0
        else:
            return self.NCPR()**2 / self.FCR()

    #...................................................................................#
    def deltaForm(self, bloblen):
        """ Calculate the delta value as defined in REF 1
        """

        sigma = self.sigma()
        nblobs = self.len - bloblen + 1
        ans = 0

        for i in xrange(0, nblobs):

            # get the blob charge pattern list
            blob = self.chargePattern[i:(i + bloblen)]

            # calculate a bunch of parameters for the blob
            # with the blob sigma value being the ultimate
            # goal

            bpos = np.where(blob > 0)[0].size
            bneg = np.where(blob < 0)[0].size

            bncpr = (bpos - bneg) / (bloblen + 0.0)
            bfcr = (bpos + bneg) / (bloblen + 0.0)

            if(bfcr == 0):
                bsig = 0
            else:
                bsig = bncpr**2 / bfcr

            # calculate the square deviation of the
            # blob sigma from the sequence sigma and
            # weight by the number of blobs in the sequence
            ans += (sigma - bsig)**2 / nblobs

        return ans

    #...................................................................................#
    def delta(self):
        """
        Return the average delta value from blobs of size 5 and 6, as
        arrived upon as being optimal in REF 1
        """

        return (self.deltaForm(5) + self.deltaForm(6)) / 2

    #...................................................................................#
    def deltaMax(self):
        """
        Return the maximum possible delta value for the sequence
        """

        # If this has been computed already, then return it
        if(self.dmax != -1):
            return self.dmax
        elif(self.FCR() == 0):
            self.dmax = 0

        #################################################################
        # FIRST computational trick - if only positive or negative
        elif self.countPos() == 0 or self.countNeg() == 0:
            nneuts = self.countNeut()

            # construct a neutral  block
            neutralBlock = '0' * nneuts

            # construct a charged block
            if self.countPos() == 0:
                # make a negative charge block
                chargedBlock = '-' * self.countNeg()
                chargeV = "-"
                ncharge = self.countNeg()

            else:
                # make a positive charge block
                chargedBlock = '+' * self.countPos()
                chargeV = "+"
                ncharge = self.countPos()

            # if the charge block is shorter than the neutral block
            if len(neutralBlock) > len(chargedBlock):

                for position in xrange(0, (self.len - ncharge) + 1):

                    setupSequence = position * "0" + \
                        chargedBlock + "0" * (nneuts - position)

                    if not len(setupSequence) == self.len:
                        raise SequenceException(
                            "Error in DeltaMax calculation")

                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])
            # if the neutral block is shorter
            else:

                for position in xrange(0, (self.len - nneuts) + 1):

                    setupSequence = position * chargeV + \
                        neutralBlock + chargeV * (ncharge - position)

                    if not len(setupSequence) == self.len:
                        raise SequenceException(
                            "Error in DeltaMax calculation")

                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])

        #################################################################
        # Second computational trick (Maximum Charge Separation)
        # If we have no neutral residues this is easy
        elif(self.countNeut() == 0):

            nNeg = self.countNeg()
            nPos = self.countPos()

            posBlock = "+" * nPos
            negBlock = "-" * nNeg

            if len(posBlock) > len(negBlock):
                for position in xrange(0, (self.len - nNeg) + 1):
                    setupSequence = position * "+" + \
                        negBlock + "+" * (nPos - position)

                    if not len(setupSequence) == self.len:
                        raise SequenceException(
                            "Error in DeltaMax calculation")

                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])
            else:
                for position in xrange(0, (self.len - nPos) + 1):
                    setupSequence = position * "-" + \
                        posBlock + "-" * (nNeg - position)

                    if not len(setupSequence) == self.len:
                        raise SequenceException(
                            "Error in DeltaMax calculation")

                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])

        #################################################################
        # Third computational trick (Maximization of # of Charged Blobs)
        # relevant if we have 18 or more neutral residues
        elif(self.countNeut() >= 18):
            nneuts = self.countNeut()

            # construct positive and negative blocks

            posBlock = '+' * self.countPos()
            negBlock = '-' * self.countNeg()

            # There are three regions where you can have neutral residues
            # to maximize delta: start/middle/end
            # This set of loops iterates through the various relevant combinations
            # of neutral residues at those locations. Importantly, this depends
            # on a blob length of 5-6!

            for startNeuts in xrange(0, 7):
                for endNeuts in xrange(0, 7):
                    setupSequence = ''
                    midBlock = ''
                    endBlock = ''

                    setupSequence = '0' * startNeuts
                    endBlock = '0' * endNeuts
                    midBlock = '0' * (nneuts - startNeuts - endNeuts)

                    # now we construct *an* optimal sequence based on the start/mid/end
                    # permutation of neutral residues
                    setupSequence += posBlock
                    setupSequence += midBlock
                    setupSequence += negBlock
                    setupSequence += endBlock
                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])

        #################################################################
        # Fourth computational trick (Search through set of sequences that fit
        # DMAX pattern)
        else:
            posBlock = ''
            negBlock = ''
            nneuts = self.countNeut()

            # construct positive and negative blocks
            posBlock = '+' * self.countPos()
            negBlock = '-' * self.countNeg()

            # iterate through different permutations where
            #  the size of the middle number of neutral residues varies
            for midNeuts in xrange(0, nneuts + 1):
                midBlock = '0' * midNeuts

                for startNeuts in xrange(0, nneuts - midNeuts + 1):
                    setupSequence = ''

                    # construct some permutation of the sequence
                    setupSequence = '0' * startNeuts

                    setupSequence += posBlock
                    setupSequence += midBlock
                    setupSequence += negBlock

                    setupSequence = setupSequence + '0' * \
                        (nneuts - startNeuts - midNeuts)

                    nseq = Sequence(setupSequence)

                    # update self.dmax if relevant
                    self.dmax = max([nseq.delta(), self.dmax])

        return self.dmax

    #...................................................................................#
    def swapRes(self, index1, index2):
        """
        Causes the two residues indexed by index1 and index2 in the sequence
        to be swapped. This returns a new sequence object with pre-calculated dmax
        """

        if(index1 == index2):
            return Sequence(self.seq)

        # we want index 2 to be greater than index 1
        elif(index2 < index1):
            index1, index2 = index2, index1
        else:
            pass

        # create a new tempseq as a charlist and then
        # swap the characters and positions index1 and index2
        tempseq = list(self.seq)
        tempseq[index1], tempseq[index2] = tempseq[index2], tempseq[index1]

        charge1 = self.chargePattern[index1]
        charge2 = self.chargePattern[index2]

        # creates a totally new array with new values - do we need this?
        tempChargeSeq = cp.deepcopy(self.chargePattern)
        tempChargeSeq[index1] = charge2
        tempChargeSeq[index2] = charge1

        # returns a new sequence with dmax and the chargeSequence already
        # generated
        return Sequence(''.join(tempseq), self.dmax, tempChargeSeq)

    #...................................................................................#
    def swapRandChargeRes(self, frozen=set()):
        """ Function which randomly selects two residues and swaps them if that
            swap would change the kappa value
        """

        # get a random number
        rand = rng.Random()
        rand.seed(time.time())

        # determine the indices from which we can swap

        # (i.e. all positive indices which do not overlap with the set of frozen
        # residues)
        posInd = set(np.where(self.chargePattern > 0)[0]) - frozen
        negInd = set(np.where(self.chargePattern < 0)[0]) - frozen
        neutInd = set(np.where(self.chargePattern == 0)[0]) - frozen

        if(len(neutInd) == 0):
            if(len(posInd) == 0 or len(negInd) == 0):
                status_message(
                    'swap will not change kappa, only one charge type in sequence')
                return self
            else:
                chargeType = [1, 2]
        elif(len(negInd) == 0):
            if(len(posInd) == 0 or len(neutInd) == 0):
                status_message(
                    'swap will not change kappa, only one charge type in sequence')
                return self
            else:
                chargeType = [1, 3]
        elif(len(posInd) == 0):
            if(len(negInd) == 0 or len(neutInd) == 0):
                status_message(
                    'swap will not change kappa, only one charge type in sequence')
                return self
            else:
                chargeType = [2, 3]
        else:
            chargeType = rand.sample([1, 2, 3], 2)

        if(chargeType[0] == 1):
            swapPair1 = rand.sample(posInd, 1)
        elif(chargeType[0] == 2):
            swapPair1 = rand.sample(negInd, 1)
        elif(chargeType[0] == 3):
            swapPair1 = rand.sample(neutInd, 1)

        if(chargeType[1] == 1):
            swapPair2 = rand.sample(posInd, 1)
        elif(chargeType[1] == 2):
            swapPair2 = rand.sample(negInd, 1)
        elif(chargeType[1] == 3):
            swapPair2 = rand.sample(neutInd, 1)
        return self.swapRes(swapPair1[0], swapPair2[0])

    #...................................................................................#
    def full_shuffle(self, frozen=set()):
        """
           Function which totally shuffles a sequences, but keeps the positions
           in the frozen set in their position

        """

        moveable_indicies = set(np.arange(0, self.len)) - set(frozen)

        lookup = {}
        index = 0
        for i in self.seq:
            lookup[index] = i
            index = index + 1

        # lookup now corresponds to an index->AA lookupyable with lookup
        # indexed

        newseq = list(moveable_indicies)
        rand = rng.Random()
        rand.seed(time.time())

        rand.shuffle(newseq)
        sequencelist = list(self.seq)

        new_seq = []

        for i in xrange(0, self.len):
            if i in frozen:
                new_seq.append(lookup[i])
            else:
                new_seq.append(lookup[newseq.pop()])

        return Sequence("".join(new_seq), self.dmax)

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                         PHOSPHORYLATION RELATED FUNCTIONS
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    # The functions below focus on phosphorylation and the like
    #

    #...................................................................................#
    def setPhosPhoSites(self, listOfPsites):
        """
        Set one or more sites on your sequence which can be phosphorylated. Note that
        this indexes from 1 (like all of bioinformatics) and not from 0 (like all of
        computer science).

        i.e. "KKKYKKK" the Y here is at position 4

        Internally we do translate to indexing from 0, but this is not something you
        should have to worry about

        Note that all data validation for the phosphosite list is done in this function.

        """

        # if we passed a single value not in a list then convert to a list of length
        # 1
        if isinstance(listOfPsites, int):
            tmp = listOfPsites
            listOfPsites = []
            listOfPsites.appned(tmp)

        # evaluate proposed phosphosites
        for site in listOfPsites:

            # check we can convert to an integer!
            site = int(site)

            # python indexes from 0 but humans from 1
            idx = site - 1

            # if we're outside our sequence
            if idx >= len(self.seq) or idx < 0:
                warning_message("Proposed phosphosite (" + str(idx + 1) +
                                " is outside sequence range. Skipping...")
                pass

            # grab the residue letter from the sequence
            res = self.seq[idx]

            status_message("Setting " + res + str(idx + 1))

            if res not in ["S", "T", "Y"]:
                # we skip it if it seems like an unphosphorylatable residue
                warning_message(
                    'Position ' +
                    str(site) +
                    ' in sequence is a non phosphorylatable residue [' +
                    str(res) +
                    ']')
            else:
                if idx in self.phosphosites:
                    # don't add the same residue twice, but no need to warn
                    # about it
                    pass
                else:
                    # let's add that bad-boy!
                    self.phosphosites.append(idx)

    #...................................................................................#
    def calculateNumberDifferentPhosphoStates(self):
        """ Function which returns the number of different phosphorylation
            states assuming each phosphosite can be phosphorylated independently
            of each other.

            This can be useful if you want to estimate how long a complete
            phoshostate distribution calculation might take (each calc is about
            0.2-1 second, depending on your CPU)
        """
        return np.power(2, len(self.phosphosites))

    #...................................................................................#
    def clear_phosphosites(self):
        """ Function which zeros out all the stored phosphosite data """
        self.phosphosites = []

    #...................................................................................#
    def calculateKappaDistOfPhosphoStates(self):
        """ Function which, using the previously defined possible phosphosites, calculates the cloud
            of possible kappa values for all combination of phosphorylation. This is a relativly long
            operation due to the fact that the FULL kappa must be recalculated everytime (i.e. delta max
            must be recomputed each time, unlike with sequence permutants where DMAX remains constant).

            To aid with this, every 50 calculations the progress is printed

        """
        ##
        #----------------------------------------------------------------------
        # local function
        def print_progress(count, total):
            if (count % 50) == 0:
                status_message("Done " + str(count) + " of " + str(total))

        num_calcs = self.calculateNumberDifferentPhosphoStates()
        #----------------------------------------------------------------------

        # itertools.product("01"..) produces all variants of a list where each
        # element is either 0 or 1 using an iterator
        phosphokappa = []
        count = 0
        for phosphostatus in itertools.product(
                "01", repeat=len(self.phosphosites)):
            newseq = list(self.seq)

            indx = 0

            # look over each element in our phosphosite on/off list
            for i in phosphostatus:
                # if that element is ON
                if int(i) == 1:
                    # set the AA at that position to a negative residue (we use E but
                    # could be D)
                    newseq[self.phosphosites[indx]] = "E"
                indx = indx + 1
            # now we've replaced some number of T/Y/R with E representing a different
            # phosphostate
            newseq = "".join(newseq)
            newseqObj = Sequence(newseq)

            # now add that new kappa to the list, update, and repeat!
            phosphokappa.append(
                (newseqObj.kappa(),
                 newseqObj.Fplus(),
                 newseqObj.Fminus(),
                 newseqObj.FCR(),
                 newseqObj.NCPR(),
                 newseqObj.meanHydropathy(),
                 phosphostatus))

            print_progress(count, num_calcs)
            count = count + 1

        return phosphokappa

    #...................................................................................#
    def get_phosphosites(self):
        """ Returns the list of phosphosites (indexed from 1
            not zero)
        """
        newSites = []
        for i in self.phosphosites:
            newSites.append(i + 1)
        return newSites

    #...................................................................................#
    def get_phosphosequence(self):
        """
        Returns a sequence with all phosphorylated residues converted to E.

        Phosphorylated residues defined by the phosphosites object variable

        """

        if len(self.phosphosites) == 0:
            warning_message(
                "No phosphosites defined - phosphosequence will be equivalent to the unphosphorylated sequence")

        # first defined the empty sequence
        pseq = ""
        idx = 0

        # for each position if that residue is phosphorylatable set it to 'E' instead
        # of the actual Y/S/T
        for i in self.seq:
            if idx in self.phosphosites:

                # extra level of checking!
                if i not in ["S", "Y", "T"]:
                    raise SequenceException(
                        "In get_phosphosequence - trying to replace non-phsophrylatable residue with GLU")

                pseq = pseq + "E"
            else:
                pseq = pseq + self.seq[idx]

            idx = idx + 1
        return pseq

    #...................................................................................#
    def kappa_at_maxPhos(self):
        """
        Returns the kappa if all the phosphosites
        were phosphorylated
        """

        # no phosphosites defined
        if len(self.phosphosites) == 0:
            return self.kappa()
        else:
            newseq = list(self.seq)
            for pos in self.phosphosites:
                newseq[pos] = "E"
            newseq = "".join(newseq)
            newseqObj = Sequence(newseq)

            return newseqObj.kappa()

    #...................................................................................#
    def get_STY_residues(self):
        """
        Where are all the Ser/Thr/Tyr residues in your
        sequence.

        Note we return positions which index from 1 not from 0

        """

        sites = []
        idx = 1
        for i in self.seq:
            if i in ["Y", "S", "T"]:
                sites.append(idx)
            idx = idx + 1
        return sites

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #                               FORMATTING FUNCTION
    #
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    #...................................................................................#
    def toString(self):
        """


        """
        s = "%i\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f\t%3.5f" % (self.len, self.Fminus(), self.Fplus(
        ), self.FCR(), self.NCPR(), self.sigma(), self.delta(), self.deltaMax(), self.kappa(), self.meanHydropathy())
        return s

    #...................................................................................#
    def set_HTMLColorResiduePalette(self, colorDict):
        """
        Function which lets you userdefine the color pallete being used with a dictionary of [AA]="color" mappings where
        color must be one of the 17 key HTML colours;

        ['aqua', 'black', 'blue', 'fuchsia', 'gray', 'green', 'lime', 'maroon', 'navy', 'olive', 'orange', 'purple',
        'red', 'silver', 'teal', 'white', 'yellow']

        """

        valid = {}

        # for each one letter amino acid code
        for i in aminoacids.ONE_TO_THREE:

            # check all the amino acids were there
            if i not in colorDict:
                raise SequenceException(
                    "When trying to set amino acid color palette the amino acid " +
                    str(i) +
                    " was missing from the input dictionary")

            # check if the color is real
            if colorDict[i] not in [
                    'aqua',
                    'black',
                    'blue',
                    'fuchsia',
                    'gray',
                    'green',
                    'lime',
                    'maroon',
                    'navy',
                    'olive',
                    'orange',
                    'purple',
                    'red',
                    'silver',
                    'teal',
                    'white',
                    'yellow']:
                raise SequenceException(
                    "When trying to set amino acid color palette the color " + str(
                        colorDict[i]) + " was used, which is not one of the 17 standard HTML colours")

            # if we get here update the 'valid' dictionary
            valid[i] = colorDict[i].lower()

        # if we got here we got through the full 20 amino acids succesfully
        self.aminoAcidColorMap = {}
        for i in valid:
            self.aminoAcidColorMap[i] = valid[i]

    #...................................................................................#
    def get_HTMLColorString(self):
        """
        Function which creates a <p> </p> wrapped HTML string with the sequence colored according to the defined
        aminoAcidColorMap.

        """
        colorString = '<p style="font-family:Courier;">'
        count = -1
        for residue in self.seq:
            count = count + 1
            if(np.mod(count, 10) == 0):
                colorString = colorString + " "
            if(np.mod(count, 50) == 0):
                colorString = colorString + "<br>"

            color = self.aminoAcidColorMap[residue]

            colorString = '%s<span style="color:%s">%s</span>' % (
                colorString, color, residue)
        colorString = colorString + "</p>"
        return colorString

    #...................................................................................#
    def toFileString(self):
        s = "Sequence   :\t%s\n\n" % (self.seq)
        s += "N          :\t%i\n" % (self.len)
        s += "f-         :\t%3.5f\n" % (self.Fminus())
        s += "f+         :\t%3.5f\n" % (self.Fplus())
        s += "FCR        :\t%3.5f\n" % (self.FCR())
        s += "NCPR       :\t%3.5f\n" % (self.NCPR())
        s += "Kappa      :\t%3.5f\n" % (self.kappa())
        s += "Sigma      :\t%3.5f\n" % (self.sigma())
        s += "Delta      :\t%3.5f\n" % (self.delta())
        s += "Max Delta  :\t%3.5f\n" % (self.deltaMax())
        s += "Hydropathy :\t%3.5f\n\n" % (self.meanHydropathy())
        s += "Phase Plot Region: %i\n" % (self.phasePlotRegion())
        s += "Phase Plot Annotation: %s\n" % (self.phasePlotAnnotation())
        return s
    

    #...................................................................................#
    def __check_window_to_length(self, bloblen):
        """
        Check to ensure that a window size is not greater than the sequence length
        which will cause issues...


        """
        if len(self.seq) < bloblen:
            raise SequenceException('Trying to use a window/blob size of %i but the sequence is only %i residues in length! Window size must be < the sequence length.' % (bloblen, len(self.seq)))

                                
    #...................................................................................#            
    def __parse_group(self, localgrp):
        """ Helper function to take a list of AAs and make sure they're upper case and valid AAs"""

        try:
            localgrp = set([x.upper() for x in localgrp])
        except AttributeError:
            raise SequenceException("ERROR: There was a problem when parsing the amino acid grouping set [%s] - this should be a list of of native AA one letter codes" % (localgrp))
                            
        for res in localgrp:
            if res not in aminoacids.TWENTY_AAs:
                raise SequenceException("ERROR: Found non-natural amino acid [%s] in grouping " % (res))

        return localgrp

