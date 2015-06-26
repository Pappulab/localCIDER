""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.7                                                         !
   !                                                                          !
   !    Copyright (C) 2014, The localCIDER development team (current and      !
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
   ! MAIN AUTHOR:   Mary Richardson and Alex Holehouse                        !
   !                                                                          !
   !--------------------------------------------------------------------------!

   
   File Description:
   ================

"""
import zlib
import math
from data.highComplexitySequences import maxComplexity

from localciderExceptions import SequenceException, SequenceComplexityException

class SequenceComplexity:
    """
    Dedicated class to dealing with sequence complexity issues.

    This is a stateless class, but contains all the logic for carrying out
    sequence complexity calculations.

    """

    def __init__(self):
        pass
        

        

    #...................................................................................#
    def raw_compressed_complexity(self, sequence):
        """ Return the normalized sequence complexity, where normalization
            occurs relative to a random string of the same length.

            Note that we are explicitly defining sequence complexity based on the ability
            to minimally encode that specific sequence.
        """
        
        if len(sequence) > 1000:
            raise SequenceException("Currently only sequences less than 1000 residues are subject to sequence complexity analysis - also, this feature is actually not officially out yet!")
            
        return float(len(zlib.compress(sequence)))/maxComplexity[len(sequence)]

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
        two      = ['L','E']
        three    = ['L','F','E']        
	four     = ['L','A','F','E']
	five     = ['L','A','F','E','K']
	six      = ['L','A','P','F','E','K']
	eight    = ['L','A','S','P','F','E','K','H']
	ten      = ['L','C','A','G','S','P','F','E','K','H']
        eleven   = ['L','C','A','G','S','P','F','E','K','H','Q']
	twelve   = ['L','C','A','G','S','P','F','W','E','D','K','H']
	fifteen  = ['L','C','A','G','S','T','P','F','W','E','Q','D','N','K','H']
	eighteen = ['L','V','C','A','G','S','T','P','F','Y','W','E','D','N','Q','K','R','H']
	twenty   = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']

    
	aa = []
	alphabet = []


        ##
        ## USER DEFINED ALPHABET
        ##
        
        if len(userAlphabet) > 0:            
            # user defined alphabet being used
            for x in ['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']:
                if not x in userAlphabet:
                    raise SequenceComplexityException('Invalid user alphabet supplied - does not allow mapping of amino acid %s'%x)

            # build the reduced sequence
            for x in sequence:
                aa.append(userAlphabet[x])
            
            # finally build the alphabet
            alphabet = []
            for x in ['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']:
                val = userAlphabet[x]
                if val not in alphabet:
                    alphabet.append(val)

            return (aa, alphabet)
        
        ##
        ## PREDEFINED DEFINED ALPHABET
        ##

        if alphabetSize not in [2,3,4,5,6,8,10,11,12,15,18,20]:
            raise SequenceComplexityException('Predefined alphabet sizes must be one of (2,3,4,5,6,8,10,11 12,15,18,20)')

        
        alphabetSize = int(alphabetSize)
            
	#2: [(LVIMCAGSTPFYW), (EDNQKRH)]
        if (alphabetSize == 2):
            alphabet = two
            for x in sequence:
                if x in ('L','V','I','M','C','A','G','S','T','P','F','Y','W'):
                    aa.append('L')
                else:
                    aa.append('E')

        #3: [(LVIMCAGSTP), (FYW), (EDNQKRH)]
        elif (alphabetSize == 3):
            alphabet = three
            for x in sequence:
                if x in ('L','V','I','M','C','A','G','S','T','P'):
                    aa.append('L')
                elif x in ('F','Y','W'):
                    aa.append('F')
                else:
                    aa.append('E')

        #4: [(LVIMC), (AGSTP), (FYW), (EDNQKRH)]
	elif (alphabetSize == 4):
            alphabet = four
            for x in sequence:
                if x in ('L','V','I','M','C'):
                    aa.append('L')
                elif x in ('A','G','S','T','P'):
                    aa.append('A')
                elif x in ('F','Y','W'):
                    aa.append('F')
                else:
                    aa.append('E')

	#5: [(LVIMC), (ASGTP), (FYW), (EDNQ), (KRH)]
        elif (alphabetSize == 5):
            alphabet = five
            for x in sequence:
                if x in ('L','V','I','M','C'):
                    aa.append('L')
                elif x in ('A','S','G','T','P'):
                    aa.append('A')
                elif x in ('F','Y','W'):
                    aa.append('F')
                elif x in ('E','D','N','Q'):
                    aa.append('E')
                else:
                    aa.append('K')
        
	#6: [(LVIM), (ASGT), (PHC), (FYW), (EDNQ), (KR)]
	elif (alphabetSize == 6):
            alphabet = six
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('A','S','G','T'):
                    aa.append('A')
                elif x in ('P','H','C'):
                    aa.append('P')
                elif x in ('F','Y','W'):
                    aa.append('F')
                elif x in ('E','D','N','Q'):
                    aa.append('E')
                else:
                    aa.append('K')
        
	#8: [(LVIMC), (AG), (ST), (P), (FYW), (EDNQ), (KR), (H)]
	elif (alphabetSize == 8):
            alphabet = eight
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('A','G'):
                    aa.append('A')
                elif x in ('S','T'):
                    aa.append('S')
                elif x in ('F','Y','W'):
                    aa.append('F')
                elif x in ('E','D','N','Q'):
                    aa.append('E')
                else:
                    aa.append('K')
        
        #10: [(LVIM), (C), (A), (G), (ST), (P), (FYW), (EDNQ), (KR), (H)]
	elif (alphabetSize == 10):
            alphabet = ten
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('S','T'):
                    aa.append('S')
                elif x in ('F','Y','W'):
                    aa.append('F')
                elif x in ('E','D','N','Q'):
                    aa.append('E')
                elif x in ('K','R'):
                    aa.append('K')
                else:
                    aa.append(x)
                    
        #11: [(LVIM), (C), (A), (G), (ST), (P), (FYW), (ED), (NQ), (KR), (H)]
	elif (alphabetSize == 11):
            alphabet = ten
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('S','T'):
                    aa.append('S')
                elif x in ('F','Y','W'):
                    aa.append('F')
                elif x in ('E','D'):                
                    aa.append('E')
                elif x in ('N','Q'):
                    aa.append('Q')
                elif x in ('K','R'):
                    aa.append('K')
                else:
                    aa.append(x)
                    
	#12: [(LVIM), (C), (A), (G), (ST), (P), (FY), (W), (EQ), (DN), (KR), (H)]
        elif (alphabetSize == 12):
            alphabet = twelve
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('S','T'):
                    aa.append('S')
                elif x in ('F','Y'):
                    aa.append('F')
                elif x in ('E','Q'):
                    aa.append('E')
                elif x in ('D','N'):
                    aa.append('D')
                elif x in ('K','R'):
                    aa.append('K')
                else:
                    aa.append(x)
                    
        #15: [(LVIM), (C), (A), (G), (S), (T), (P), (FY), (W), (E), (Q), (D), (N), (KR), (H)]
	elif (alphabetSize == 15):
            alphabet = fifteen
            for x in sequence:
                if x in ('L','V','I','M'):
                    aa.append('L')
                elif x in ('F','Y'):
                    aa.append('F')
                elif x in ('K','R'):
                    aa.append('K')
                else:
                    aa.append(x)
                    
	#18: [(LM), (VI), (C), (A), (G), (S), (T), (P), (F), (Y), (W), (E), (D), (N), (Q), (K), (R), (H)]
        elif (alphabetSize == 18):
            alphabet = eighteen
            for x in sequence:
                if x in ('L','M'):
                    aa.append('L')
                elif x in ('V','I'):
                    aa.append('V')
                else:
                    aa.append(x)
	elif (alphabetSize == 20):
            alphabet = twenty
            aa = sequence
	else:
            print ("ERROR: invalid aa alphabet selected, using 20-letter alphabet instead")
            alphabet = twenty
            aa = sequence


	return aa, alphabet

    ###########################################################
    # calculate Wootton-Federhen complexity
    ###########################################################
    def CWF(self, sequence, alphabet, window_size, step_size):
        """
        Function to calculate the Wootton-Federhen complexity

        Requires four parameters
        1) Amino acid sequence
        2) Alphabet being used
        3) Window_size for sliding window
        4) Stepsize for moving along the sliding window

        
        """

        # Initialize the current step
	step = 0 

	CWF_array = []
        
        # for each position
	while (step <= len(sequence)-window_size):

            # restart complexity calculation for this window
            CWF = 0 #
        
            # get the current window
            window = sequence[step:step+window_size]

            # for every residue in the alphabet
            for x in alphabet: 

                # p is between 0 and 1
                p = float(window.count(x))/window_size
                
                if p > 0:
                    CWF = p*(math.log(p,len(alphabet))) + CWF
                    
            CWF_array.append(-CWF) #store the complexity score for this window
            step = step + step_size #increment the step

	return CWF_array


    def get_WF_complexity(self, sequence, alphabetSize=20, userAlphabet={}, window_size=10,step_size=1):

        # reduce alphabet complexity
        (reduced_sequence, alphabet) = self.reduce_alphabet(sequence, alphabetSize, userAlphabet)

        return self.CWF(reduced_sequence, alphabet, window_size, step_size)



    
