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



   File Description:
   ================

   backendtools contains various utilities and tools used by localCIDER. Especially
   useful is the fact all STDOUT is done through the various *_message methods
   defined here. This allows easy control over the degree of verbosity localCIDER
   displays. Verbocity is hard-coded into the config file



"""
import sys
import os
from config import HUSH_WARNINGS, HUSH_STATUS, HUSH_ALL

#...................................................................................#


def warning_message(message):
    if not HUSH_WARNINGS and not HUSH_ALL:
        print "WARNING: " + message


#...................................................................................#
def status_message(message):
    """
    Unless the HUSH_STATUS variable is set to True
    this function prints $message to STDOUT
    """
    if not HUSH_STATUS and not HUSH_ALL:
        print message


#...................................................................................#
def running_dotdotdot():
    """
    Print a single dot (period) to STDOUT without a newline.
    Useful for indicating progress.
    """
    if not HUSH_STATUS and not HUSH_ALL:
        sys.stdout.write('.')
        sys.stdout.flush()


#...................................................................................#
def warn_thisWillBeRemoved():
    """
    Prints a standard depreciation warning
    """
    warning_message("<> TO BE DEPRECATED IN THE NEAR FUTURE <>")


#...................................................................................#
def warn_notReadyYet():
    """
    Prints a standard 'this feature isn't ready yet' warning
    """
    warning_message("<> THIS FEATURE IS COMING IN THE NEXT VERSION <>")


#...................................................................................#
def return_absolute_datafile_path(filename):
    """
    Function which returns the absolute path
    of a file in the package's data directory.

    """
    warn_thisWillBeRemoved()

    # get the absolute path of where we are now
    absolute_path = os.path.realpath('__file__')
    splitted = os.path.split(absolute_path)

    print absolute_path

    # cut off the first two path
    for i in range(0, 2):
        splitted = os.path.split(splitted[0])

    return(os.path.join(splitted[0], "data", filename))


#...................................................................................#
def verifyType(obj, typeHere):
    """
    Function which takes an object and some type and ensures the
    object

    """

    try:
        objclass = obj.__class__
        if objclass == typeHere:
            return True
        else:
            return False
    except SyntaxError as e:
        # a primitive upon which .__class__ is invalid
        if isinstance(obj, typeHere):
            return True
        else:
            return False
