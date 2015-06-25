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
from data.highComplexitySequences import maxComplexity

class SequenceComplexity:
    """
    Dedicated class to dealing with sequence complexity issues.

    This is a stateless class, but contains all the logic for carrying out
    sequence complexity calculations.

    """

    def __init__(self, alphabet=20):
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

    
