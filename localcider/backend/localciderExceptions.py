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

   Contains exceptions associated with localcider

"""


######################
# backend.plotting
#...................................................................................#
class PlottingException(Exception):
    """
    Exception for the plotting functions
    """
    pass

######################
# backend.keyfile
#...................................................................................#


class KeyFileException(Exception):
    """
    Exception for the KeyFile class
    """
    pass


######################
# backend.sequence
#...................................................................................#
class SequenceException(Exception):
    """
    Exception for the sequence class

    """
    pass

######################
# backend.seqfileparser
#...................................................................................#


class SequenceFileParserException(Exception):
    pass


######################
# backend.restables
#...................................................................................#
class ResTableException(Exception):
    pass


######################
# backend.wang_landau
#...................................................................................#
class WLException(Exception):
    pass


######################
# backend.wang_landau
#...................................................................................#
class SequenceComplexityException(Exception):
    pass
