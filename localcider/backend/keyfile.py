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

   keyfile parses a localCIDER keyfile for the more complex calculations. Right
   now this is just Wang-Landau based density of kappa state calculations, but
   there's no reason why more advanced analytical techniques can't be added in
   the future.

   With version 0.1.0 the keyfile functionality is redundant, as we don't provide
   stand-alone binaries. However, that is a main priority of version 0.3.0

"""


import os
import numpy as np
from backendtools import warning_message, status_message
from seqfileparser import SequenceFileParser

from localciderExceptions import KeyFileException

#=======================================================
# Default values for non-crucial keyfile parameters
DEFAULT_VALS = {
    'BIN_MIN': 0.05,
    'BIN_MAX': 0.95,
    'NUMBER_OF_BINS': 20,
    'FLATCHECK_FREQ': 10000,
    'CONVERGENCE': np.exp(0.0001),
    'FLATNESS_CRITERION': 0.7,
    'WL_TYPE': 'ZOOM'}

# list of expected keywords
KEYWORD_LIST = [
    'SEQFILE',
    'OUTDIR',
    'FREEZE_FILE',
    'BIN_MIN',
    'BIN_MAX',
    'NUMBER_OF_BINS',
    'FLATCHECK_FREQ',
    'CONVERGENCE',
    'FLATNESS_CRITERION',
    'WL_TYPE']

# types of WL
WL_TYPES = ['NORMAL', 'ZOOM']

#======================================================


class KeyFile():
    """
        KeyFile objects take a filename and do all the reuqired keyfile parsing during initializtaion
        creating a single object from which you can get any WL keyfile parameters of interest
    """

    # SEQUENCE_FILE  - file containing sequence
    # FREEZE_FILE    - file defining regions to hold constant
    # BIN_START      - starting bin value (optional)
    # BIN_END        - ending bin value (optional)
    # NUMBER_OF BINS - number of bins
    # OUTDIR         - directory for output to be written to
    # WL             -

    #...................................................................................#
    def __init__(self, filename):
        """
           Initializer defines the keywords we expect to see, and takes a filename and parses it

        """

        # initialize and zero out the keywords in the object dictionary
        self.KEYWORDS = {}
        for i in KEYWORD_LIST:
            self.KEYWORDS[i] = ''

        # parse the keyfile
        self.parse_keyfile(filename)

    #...................................................................................#
    def parse_keyfile(self, filename):
        """
        Function which takes a filename and parses it into the keyfile object for easy
        interaction with the file's content
        """

        status_message("Parsing keyfile...")
        status_message("---------------------------------------")

        SeqFileParser = SequenceFileParser()  # create a sequence file parsing object

        # read file to end
        with open(filename) as filehandle:
            content = filehandle.readlines()

        # [PHASE 1 START]
        # PARSE THE KEYFILE
        for line in content:
            line = line.strip()

            # if empty line
            if len(line) == 0:
                continue

            # comments in the keyfile
            if line[0] == "#":
                continue

            # if inline comment kill everything after the comment
            # character
            if len(line.split("#")) > 1:
                line = line.split("#")[0]

            # finally remove any other trailing whitespace
            line = line.strip()

            # split the remaining by whitespace
            line_list = line.split(" ")

            # now cycle over the first value in the whitelist splitted lits
            # and check if it matches one of the predefined KEYWORDS
            if line_list[0].strip() in self.KEYWORDS:

                # if we find a keyword and there's a single string after the keyword load it into the
                # KEYWORDS dictionary (i.e. this is what we expect!)
                if len(line_list) == 2:
                    self.KEYWORDS[line_list[0].strip()] = line_list[1].strip()
                # there was more than one whitespace seperated string after the
                # keyword - we basically fail at this
                else:
                    raise KeyFileException("Error: Found keyword " +
                                           str(line_list[0].strip()) +
                                           " but unable to parse associated value")
            else:
                warning_message(
                    "Found unexpected keyword [" + str(line_list[0].strip()) + "] - ignorning...")

        # Now add default for the sequene, which will hopefully be set in the
        # next section by parsing the sequencefile
        self.KEYWORDS["SEQUENCE"] = ""
        # [PHASE 1 END]

        status_message("---------------------------------------")
        status_message("Keyfile parsed!\n")
        status_message("Validating keyfile contents")
        status_message("---------------------------------------")

        # Having parsed the keyfile we now validate the keyfile so
        # we don't have to worry about validation later on

        # VALIDATE the parsed values
        # [PHASE 2 START]
        #

        for keyword in KEYWORD_LIST:
            # extract the value associated with each keyword in turn
            value = self.KEYWORDS[keyword]

            ##
            # SEQUENCE FILE VALIDATION AND PARSING
            ##
            if keyword == "SEQFILE":
                if value == "":
                    raise KeyFileException(
                        "ERROR: No sequence file provided in keyfile (expecting keyword [SEQFILE])")
                else:
                    if not os.path.isfile(value):
                        raise KeyFileException(
                            "Expected " + str(value) + " to be file")

                    # if its a file lets try and extract a sequence from it!
                    self.KEYWORDS[
                        "SEQUENCE"] = SeqFileParser.parseSeqFile(value)

                    # if we get here we *should* now have a sequence...
                    if self.KEYWORDS["SEQUENCE"] == "" or self.KEYWORDS[
                            "SEQUENCE"] is None:
                        raise KeyFileException(
                            "ERROR: No sequence was parsed from the sequence file...")

            ##
            # OUTPUT DIRECTORY VALIDATION
            ##
            elif keyword == "OUTDIR":
                if value == "":
                    raise KeyFileException(
                        "ERROR: No output directory provided in keyfile (expecting keyword [OUTDIR]")
                else:
                    # creates the output directory if it doesn't already exist
                    if not os.path.exists(value):
                        status_message(
                            "Creating output directory " + str(value))
                        try:
                            os.makedirs(value)
                        except OSError as e:
                            print "----------------------------"
                            print ""
                            print "ERROR Creating output directory - do you have permission to create the directory [" + str(value) + "]"
                            print ""
                            print "----------------------------"
                            raise e
                    # if it does exist raise a quick warning
                    else:
                        # check if its empty
                        if len(os.listdir(value)) > 0:
                            warning_message(
                                "Output directory exists already and is not empty [RISK OF OVERWRITING!]")
                        else:
                            pass  # empty directory already exists - brilliant!

            ##
            # FREEZE FILE VALIDATION
            ##
            elif keyword == "FREEZE_FILE":
                if value == "":
                    pass  # no freeze file, no problem
                else:
                    if not os.path.isfile(value):
                        raise KeyFileException(
                            "Expected " + str(value) + " to be file")
                    status_message("Using freeze file")
                    self.KEYWORDS[keyword] = value

            ##
            # WL TYPE
            ##
            elif keyword == "WL_TYPE":
                if value == "":
                    self.KEYWORDS[keyword] = DEFAULT_VALS[keyword]
                    status_message(
                        "Setting WL type to default [" + str(DEFAULT_VALS[keyword]) + "]")
                else:
                    if value in WL_TYPES:
                        self.KEYWORDS[keyword] = value
                        status_message(
                            "Setting WL type to keyfile defined [" + str(value) + "]")
                    else:
                        raise KeyFileException(
                            "Unexpected WL algorithm type selected " + str(value) + " ")

                self.KEYWORDS[keyword] = value

            ##
            # SET NUMERIC VALUES
            ##

            elif keyword == "BIN_MIN":
                self.__set_numeric(keyword, value)

            elif keyword == "BIN_MAX":
                self.__set_numeric(keyword, value)

            elif keyword == "NUMBER_OF_BINS":
                self.__set_numeric(keyword, value)

            elif keyword == "FLATCHECK_FREQ":
                self.__set_numeric(keyword, value)

            elif keyword == "CONVERGENCE":
                self.__set_numeric(keyword, value)

            elif keyword == "FLATNESS_CRITERION":
                self.__set_numeric(keyword, value)

            else:
                raise KeyFileException("SHOULD NOT BE GETTING HERE...")

    #...................................................................................#
    def __set_numeric(self, keyword, value):
        """
           Function which sets the KEYWORDS dictionary $keyword value to $value if $value can be treated
           as a numerical value or uses the default if it wasn't set (BUT DOES NOT use the default
           if an in-parsable value was set - we want to know when things are going wrong, silent errors
           cost lives. Maybe.).
        """
        if value == "":
            status_message("Setting " + keyword +
                           " to default [" + str(DEFAULT_VALS[keyword]) + "]")
            self.KEYWORDS[keyword] = DEFAULT_VALS[keyword]
        else:
            try:
                float(value)
                status_message("Setting " + keyword +
                               " to keyfile defined [" + str(value) + "]")
                self.KEYWORDS[keyword] = value
            except ValueError:
                raise KeyFileException(
                    "\n\nERROR: Invalid value for " +
                    keyword +
                    " - unable to convert [" +
                    value +
                    "] into a number\n")

    #...................................................................................#
    def get_bin_min(self):
        return self.KEYWORDS["BIN_MIN"]

    #...................................................................................#
    def get_bin_max(self):
        return self.KEYWORDS["BIN_MAX"]

    #...................................................................................#
    def get_num_bins(self):
        return self.KEYWORDS["NUMBER_OF_BINS"]

    #...................................................................................#
    def get_flatcheck_freq(self):
        return self.KEYWORDS["FLATCHECK_FREQ"]

    #...................................................................................#
    def get_convergence(self):
        return self.KEYWORDS["CONVERGENCE"]

    #...................................................................................#
    def get_sequence(self):
        return self.KEYWORDS["SEQUENCE"]

    #...................................................................................#
    def get_freezefile(self):
        return self.KEYWORDS["FREEZE_FILE"]

    #...................................................................................#
    def get_outdir(self):
        return self.KEYWORDS["OUTDIR"]

    #...................................................................................#
    def get_flatness_criterion(self):
        return self.KEYWORDS["FLATNESS_CRITERION"]

    #...................................................................................#
    def get_WL_type(self):
        return self.KEYWORDS["WL_TYPE"]
