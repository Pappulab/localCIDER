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
   ! MAIN AUTHOR:   James Ahad and Alex Holehouse                             !
   !                                                                          !
   !--------------------------------------------------------------------------!


   File Description:
   ================

   Class constructs a querayable ResTable object which provides rapid
   for each amino acid's properties


"""
import csv
import numpy as np

from residue import Residue
import data
from data.aminoacids import THREE_TO_ONE, ONE_TO_THREE

from localciderExceptions import ResTableException


class ResTable:
    """ Class which holds data on amino acids as read from the data.aminoacids file

        """

    #...................................................................................#
    def __init__(self):
        """ Takes a filename and initializes the self.residue_table list
            with the contents of a residue table file
        """

        self.residue_table = {}

        residue_list = data.aminoacids.buildTable()

        for r in residue_list:
            res = Residue(r[0], r[1], r[2], r[3], r[4], r[5])
            self.residue_table[res.letterCode3] = res

    #...................................................................................#
    def lookForRes(self, resCode):
        """ Returns a Residue objected defined by the resCode (can be
            three letter or one letter - we don't care!
        """

        # if we're using a three letter code
        if len(resCode) == 1:
            if resCode in ONE_TO_THREE.keys():
                return self.residue_table[ONE_TO_THREE[resCode]]

        # if we're usin a one letter code
        elif len(resCode) == 3:
            if resCode in THREE_TO_ONE.keys():
                return self.residue_table[resCode]

        # if we got here we had an invalid AA code
        raise ResTableException(
            "Invalid amino acid code provided [" + str(resCode) + "]")

    #...................................................................................#
    def lookUpHydropathy(self, resCode):
        """
        Get an amino acid's hydropathy
        """

        res = self.lookForRes(resCode)
        return res.hydropathy

    #...................................................................................#
    def lookUpCharge(self, resCode):
        """
        Get the charge on a residue (1/-1/0)
        """

        # first check for reduced residue types (aka charge +/-/0)
        if(resCode == '+'):
            return 1
        elif(resCode == '-'):
            return -1
        elif(resCode == '0'):
            return 0
        else:
            res = self.lookForRes(resCode)
            return res.charge

    #...................................................................................#
    def lookUpPPII(self, resCode):
        """
        Get an amino acid's PPII propensity as defined by Hilser.
        """

        res = self.lookForRes(resCode)
        return res.PPII

        
