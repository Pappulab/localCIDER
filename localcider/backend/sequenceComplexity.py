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
   ! MAIN AUTHOR:   Mary Richardson and Alex Holehouse                        !
   !                                                                          !
   !--------------------------------------------------------------------------!


   File Description:
   ================

"""
import zlib
import numpy as np
import math
from data.highComplexitySequences import maxComplexity

from localciderExceptions import SequenceException, SequenceComplexityException
from data.aminoacids import TWENTY_AAs


class SequenceComplexity:
    """
    Dedicated class to dealing with sequence complexity issues.

    This is a stateless class, but contains all the logic for carrying out
    sequence complexity calculations.

    """

    def __init__(self):
        pass



    #...................................................................................#
    def Zlib_compressed_complexity(self, sequence):
        """ Return the normalized sequence complexity, where normalization
            occurs relative to a random string of the same length.

            Note that we are explicitly defining sequence complexity based on the ability
            to minimally encode that specific sequence.
        """

        if len(sequence) > 1000:
            raise SequenceException(
                "Currently only sequences less than 1000 residues are subject to sequence complexity analysis - also, this feature is actually not officially out yet!")

        return float(len(zlib.compress(sequence))) / \
            maxComplexity[len(sequence)]



    #...................................................................................#
    def reduce_alphabet(self, sequence, alphabetSize=20, userAlphabet={}):
        """
        Converts your sequence to either one of the 11 predefined alphabets (defined based
        on size) or allows you to define your own alphabet using the userAlphabet dictionary.

        userAlphabet is checked for sanity (i.e. that it can deal with all 20 amino acids)

        Predefined alphabets are shown below

        two      - [(LVIMCAGSTPFYW), (EDNQKRH)]
        three    - [(LVIMCAGSTP), (FYW), (EDNQKRH)]
        four     - [(LVIMC), (AGSTP), (FYW), (EDNQKRH)]
        five     - [(LVIMC), (ASGTP), (FYW), (EDNQ), (KRH)]
        six      - [(LVIM), (ASGT), (PHC), (FYW), (EDNQ), (KR)]
        eight    - [(LVIMC), (AG), (ST), (P), (FYW), (EDNQ), (KR), (H)]
        ten      - [(LVIM), (C), (A), (G), (ST), (P), (FYW), (EDNQ), (KR), (H)]
        eleven   - [(LVIM), (C), (A), (G), (ST), (P), (FYW), (ED), (NQ), (KR), (H)]
        twelve   - [(LVIM), (C), (A), (G), (ST), (P), (FY), (W), (EQ), (DN), (KR), (H)]
        fifteen  - [(LVIM), (C), (A), (G), (S), (T), (P), (FY), (W), (E), (Q), (D), (N), (KR), (H)]
        eighteen - [(LM), (VI), (C), (A), (G), (S), (T), (P), (F), (Y), (W), (E), (D), (N), (Q), (K), (R), (H)]
        twenty   - all twenty!


        """

        two    = ['L', 'E']
        three  = ['L', 'F', 'E']
        four   = ['L', 'A', 'F', 'E']
        five   = ['L', 'A', 'F', 'E', 'K']
        six    = ['L', 'A', 'P', 'F', 'E', 'K']
        eight  = ['L', 'A', 'S', 'P', 'F', 'E', 'K', 'H']
        ten    = ['L', 'C', 'A', 'G', 'S', 'P', 'F', 'E', 'K', 'H']
        eleven = ['L', 'C', 'A', 'G', 'S', 'P', 'F', 'E', 'K', 'H', 'Q']
        twelve = ['L', 'C', 'A', 'G', 'S', 'P', 'F', 'W', 'E', 'D', 'K', 'H']
        fifteen = [
            'L',
            'C',
            'A',
            'G',
            'S',
            'T',
            'P',
            'F',
            'W',
            'E',
            'Q',
            'D',
            'N',
            'K',
            'H']
        eighteen = ['L', 'V', 'C', 'A', 'G', 'S', 'T', 'P',
                    'F', 'Y', 'W', 'E', 'D', 'N', 'Q', 'K', 'R', 'H']
        twenty = [
            'R',
            'H',
            'K',
            'D',
            'E',
            'S',
            'T',
            'N',
            'Q',
            'C',
            'G',
            'P',
            'A',
            'I',
            'L',
            'M',
            'F',
            'W',
            'Y',
            'V']

        aa = []
        alphabet = []

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # USER DEFINED ALPHABET
        #
        # Note that if we're using a user defined alphabet then we ignore the alphabetSize
        # parameter entirly

        if len(userAlphabet) > 0:
            if not isinstance(userAlphabet, dict):
                raise SequenceComplexityException(
                    'Invalid user alphabet supplied, must be a dictionary or dictionary derived object')

            # this is all sanity checking for our userdefined alphabet, but is important
            # because most of the other functions assume the sequence being
            # passed is valid
            for x in TWENTY_AAs:
                try:
                    converted = userAlphabet[x]
                except KeyError:
                    raise SequenceComplexityException(
                        'Invalid user alphabet supplied - does not allow mapping of amino acid %s' % x)

                if converted not in TWENTY_AAs:
                    raise SequenceComplexityException(
                        'Invalid user alphabet supplied - amino acid %s maps to %s, which is not a valid amino acid (must be upper case)' %
                        (x, converted))

            # build the reduced sequence
            for x in sequence:
                aa.append(userAlphabet[x])

            # finally build the alphabet
            alphabet = []
            for x in TWENTY_AAs:
                val = userAlphabet[x]
                if val not in alphabet:
                    alphabet.append(val)

            return ("".join(aa), alphabet)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # PREDEFINED DEFINED ALPHABET
        #
        # NOTE - we don't get here if we've used a userdefined alphabet
        # check the type of alphabetSize is OK
        try:
            alphabetSize = int(alphabetSize)
        except ValueError:
            raise SequenceComplexityException(
                "alphabetSize should be an number, or a string which can be converted into an integer")

        # check that one of the known alphabets is being used
        if alphabetSize not in [2, 3, 4, 5, 6, 8, 10, 11, 12, 15, 18, 20]:
            raise SequenceComplexityException(
                'Predefined alphabet sizes must be one of (2,3,4,5,6,8,10,11 12,15,18,20)\nGot [%i]'%alphabetSize)

        # Now we figure out which alphabet we're using and then generate a reduced resolution
        # sequence
        # 2: [(LVIMCAGSTPFYW), (EDNQKRH)]
        if (alphabetSize == 2):
            alphabet = two
            for x in sequence:
                if x in (
                        'L',
                        'V',
                        'I',
                        'M',
                        'C',
                        'A',
                        'G',
                        'S',
                        'T',
                        'P',
                        'F',
                        'Y',
                        'W'):
                    aa.append('L')
                else:
                    aa.append('E')

        # 3: [(LVIMCAGSTP), (FYW), (EDNQKRH)]
        elif (alphabetSize == 3):
            alphabet = three
            for x in sequence:
                if x in ('L', 'V', 'I', 'M', 'C', 'A', 'G', 'S', 'T', 'P'):
                    aa.append('L')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                else:
                    aa.append('E')

        # 4: [(LVIMC), (AGSTP), (FYW), (EDNQKRH)]
        elif (alphabetSize == 4):
            alphabet = four
            for x in sequence:
                if x in ('L', 'V', 'I', 'M', 'C'):
                    aa.append('L')
                elif x in ('A', 'G', 'S', 'T', 'P'):
                    aa.append('A')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                else:
                    aa.append('E')

        # 5: [(LVIMC), (ASGTP), (FYW), (EDNQ), (KRH)]
        elif (alphabetSize == 5):
            alphabet = five
            for x in sequence:
                if x in ('L', 'V', 'I', 'M', 'C'):
                    aa.append('L')
                elif x in ('A', 'S', 'G', 'T', 'P'):
                    aa.append('A')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                elif x in ('E', 'D', 'N', 'Q'):
                    aa.append('E')
                else:
                    aa.append('K')

        # 6: [(LVIM), (ASGT), (PHC), (FYW), (EDNQ), (KR)]
        elif (alphabetSize == 6):
            alphabet = six
            for x in sequence:
                if x in ('L', 'V', 'I', 'M'):
                    aa.append('L')
                elif x in ('A', 'S', 'G', 'T'):
                    aa.append('A')
                elif x in ('P', 'H', 'C'):
                    aa.append('P')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                elif x in ('E', 'D', 'N', 'Q'):
                    aa.append('E')
                else:
                    aa.append('K')

        # 8: [(LVIMC), (AG), (ST), (P), (FYW), (EDNQ), (KR), (H)]
        elif (alphabetSize == 8):
            alphabet = eight
            for x in sequence:
                if x in ('L', 'V', 'I', 'M', 'C'):
                    aa.append('L')
                elif x in ('A', 'G'):
                    aa.append('A')
                elif x in ('S', 'T'):
                    aa.append('S')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                elif x in ('E', 'D', 'N', 'Q'):
                    aa.append('E')
                elif x in ('H'):
                    aa.append('H')
                elif x in ('P'):
                    aa.append('P')
                else:
                    aa.append('K')

        # 10: [(LVIM), (C), (A), (G), (ST), (P), (FYW), (EDNQ), (KR), (H)]
        elif (alphabetSize == 10):
            alphabet = ten
            for x in sequence:
                if x in ('L', 'V', 'I', 'M'):
                    aa.append('L')
                elif x in ('S', 'T'):
                    aa.append('S')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                elif x in ('E', 'D', 'N', 'Q'):
                    aa.append('E')
                elif x in ('K', 'R'):
                    aa.append('K')
                else:
                    aa.append(x)

        # 11: [(LVIM), (C), (A), (G), (ST), (P), (FYW), (ED), (NQ), (KR), (H)]
        elif (alphabetSize == 11):
            alphabet = eleven
            for x in sequence:
                if x in ('L', 'V', 'I', 'M'):
                    aa.append('L')
                elif x in ('S', 'T'):
                    aa.append('S')
                elif x in ('F', 'Y', 'W'):
                    aa.append('F')
                elif x in ('E', 'D'):
                    aa.append('E')
                elif x in ('N', 'Q'):
                    aa.append('Q')
                elif x in ('K', 'R'):
                    aa.append('K')
                else:
                    aa.append(x)

        # 12: [(LVIM), (C), (A), (G), (ST), (P), (FY), (W), (EQ), (DN), (KR),
        # (H)]
        elif (alphabetSize == 12):
            alphabet = twelve
            for x in sequence:
                if x in ('L', 'V', 'I', 'M'):
                    aa.append('L')
                elif x in ('S', 'T'):
                    aa.append('S')
                elif x in ('F', 'Y'):
                    aa.append('F')
                elif x in ('E', 'Q'):
                    aa.append('E')
                elif x in ('D', 'N'):
                    aa.append('D')
                elif x in ('K', 'R'):
                    aa.append('K')
                else:
                    aa.append(x)

        # 15: [(LVIM), (C), (A), (G), (S), (T), (P), (FY), (W), (E), (Q), (D),
        # (N), (KR), (H)]
        elif (alphabetSize == 15):
            alphabet = fifteen
            for x in sequence:
                if x in ('L', 'V', 'I', 'M'):
                    aa.append('L')
                elif x in ('F', 'Y'):
                    aa.append('F')
                elif x in ('K', 'R'):
                    aa.append('K')
                else:
                    aa.append(x)

        # 18: [(LM), (VI), (C), (A), (G), (S), (T), (P), (F), (Y), (W), (E),
        # (D), (N), (Q), (K), (R), (H)]
        elif (alphabetSize == 18):
            alphabet = eighteen
            for x in sequence:
                if x in ('L', 'M'):
                    aa.append('L')
                elif x in ('V', 'I'):
                    aa.append('V')
                else:
                    aa.append(x)
        elif (alphabetSize == 20):
            alphabet = twenty
            aa = sequence
        else:
            print (
                "ERROR: invalid aa alphabet selected, using 20-letter alphabet instead")
            alphabet = twenty
            aa = sequence

        return ("".join(aa), alphabet)
    
        

    #...................................................................................#
    ###########################################################
    # calculate Wootton-Federhen complexity
    ###########################################################
    def CWF(self, sequence, alphabet, windowSize, stepSize):
        """
        Function to calculate the Wootton-Federhen complexity

        Requires four parameters
        1) Amino acid sequence
        2) Alphabet being used
        3) windowSize for sliding window
        4) Stepsize for moving along the sliding window

        """

        # Initialize the current step
        step = 0

        CWF_array = []

        # for each position
        while (step <= len(sequence) - windowSize):

            # restart complexity calculation for this window
            CWF = 0

            # get the current window
            window = sequence[step:step + windowSize]

            # for every residue in the alphabet
            for x in alphabet:

                # p is between 0 and 1
                p = float(window.count(x)) / windowSize

                if p > 0:
                    CWF = p * (math.log(p, len(alphabet))) + CWF

            # store the complexity score for this window
            CWF_array.append(-CWF)
            step = step + stepSize  # increment the step

        return CWF_array



    #...................................................................................#
    ###########################################################
    # calculate linguistic complexity
    ###########################################################
    def LC(self, sequence, alphabet, windowSize, stepSize, wordSize):
        # Simple explanation of linguistic complexity
        #
        # Linguistic complexity is basically asking, given a window of some size,
        # lets figure out the number of unique 'words' of some (specific) size in
        # that window, and divide that by the maximum possible number of unique words
        # in the window. We do this for each window, to create a vectorial
        # representation of the complexity

        # the current step
        step = 0
        LC_array = []

        while (step <= len(sequence) - windowSize):

            # restart complexity calculation for this window
            LC = 0

            # reset the position for this window
            i = 0

            ngrams = set()

            ngram = ''

            # for each position in the window we step through with a word of
            # size wordSize
            for i in range(0, windowSize - wordSize):

                # extract the current ngram
                position = step + i

                ngram = ''.join(sequence[position:position + wordSize])

                # if this ngram is not already present in the set
                if ngram not in ngrams:
                    # add it to the ngrams set
                    ngrams.add(ngram)

            # size of ngrams set
            v = len(ngrams)

            # max possible number of words
            vmax = min(len(alphabet)**wordSize, windowSize - 1 + wordSize)

            # linguistic complexity associated with this window
            LC = float(v) / vmax

            # add this to an array of the complexity profile scores and
            # increment the step
            LC_array.append(LC)
            step = step + stepSize

        return LC_array



    #...................................................................................#
    ###########################################################
    # calculate Lempel-Ziv-Welch complexity
    ###########################################################
    def LZW(self, sequence, alphabet, windowSize, stepSize):
        # the current step
        step = 0

        LZW_array = []

        while (step <= len(sequence) - windowSize):
            LZW = 0  # restart complexity calculation for this window
            i = 0   # reset the position for this window
            w = ''  # reset the word

            ngrams = set()

            # for each position in the window
            for i in range(0, windowSize):
                position = step + i

                # if we find a word which is already in
                # the ngrams set, update w to that
                if (w + sequence[position] in ngrams):

                    w = sequence[position] + w
                    #w = w+sequence[position]

                # else we found a new word
                else:
                    ngrams.add(w + sequence[position])
                    w = sequence[position]

            # size of ngrams set
            n = len(ngrams)

            if (windowSize > 0):
                LZW = float(n) / windowSize
                # add this to an array of the complexity profile scores
                LZW_array.append(LZW)
            step += stepSize         # increment the step

        return LZW_array



    #...................................................................................#
    def get_WF_complexity(
            self,
            sequence,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1):

        # reduce alphabet complexity
        (reduced_sequence, alphabet) = self.reduce_alphabet(
            sequence, alphabetSize, userAlphabet)

        complexity_vector = self.CWF(reduced_sequence, alphabet, windowSize, stepSize)

        return self.get_indexed_complexity_vector(complexity_vector, len(sequence))

        

    def get_LC_complexity(
            self,
            sequence,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1,
            wordSize=3):

        # reduce alphabet complexity
        (reduced_sequence, alphabet) = self.reduce_alphabet(
            sequence, alphabetSize, userAlphabet)

        complexity_vector =  self.LC(
            reduced_sequence,
            alphabet,
            windowSize,
            stepSize,
            wordSize)

        return self.get_indexed_complexity_vector(complexity_vector, len(sequence))

        
    def get_LZW_complexity(
            self,
            sequence,
            alphabetSize=20,
            userAlphabet={},
            windowSize=10,
            stepSize=1):

        # reduce alphabet complexity
        (reduced_sequence, alphabet) = self.reduce_alphabet(
            sequence, alphabetSize, userAlphabet)


        complexity_vector = self.LZW(reduced_sequence, alphabet, windowSize, stepSize)

        return self.get_indexed_complexity_vector(complexity_vector, len(sequence))


    def get_indexed_complexity_vector(self, complexity_vector, seq_len):
        """
        Takes a complexity vector of length N and returns a 2xN matrix containing
        the complexity vector and the corresponding residue positions distributed
        equally along the sequence of length seq_len

        """

        # having retrieved the vector we want to also return indicies which represent
        # an even distribution of points across the sequence (this is always going
        # to be the correct representation of the vector)
        complexity_vector_len = len(complexity_vector)
        spacing = seq_len/complexity_vector_len

        # based on the spacing calculate the necessary remainder (i.e. what's left
        # over..)
        remainder = (seq_len - (spacing*complexity_vector_len))

        if remainder % 2 == 0:
            flank_start = remainder/2
            flank_end = remainder/2
        else:
            flank_start = (remainder-1)/2
            flank_end = (remainder+1)/2

        # finally set the start and end positions as half-way through the spacing
        # unit + the appropriate flanking sequence offset from either the first
        # or last residue

        index_start = (flank_start+1) + spacing/2
        index_end   = ((seq_len + 1)-flank_end) + spacing/2
        indices = np.arange(index_start, index_end, spacing, dtype=int)

        return np.vstack((indices, complexity_vector))

