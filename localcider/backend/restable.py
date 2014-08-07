import csv
import numpy as np

from residue import Residue
import data
from data.aminoacids import THREE_TO_ONE, ONE_TO_THREE

class ResTableException(Exception):
    pass

class ResTable:
    """ Class which holds data on amino acids as read from an amino acid
        data table. Each line in such a file should be of the following format

        <name>, <3 letter code>, <1 letter code>, <hydropathy>, <charge>

        e.g.

        Isoleucine,ILE,I,9,0

        """
    def __init__(self):
        """ Takes a filename and initializes the self.residue_table list
            with the contents of a residue table file        
        """
        
        self.residue_table = {}

        residue_list = data.aminoacids.buildTable()

        for r in residue_list:
            res = Residue(r[0],r[1],r[2],r[3],r[4])
            self.residue_table[res.letterCode3] = res

        """
        # set the filename
        self.filename = filename

        # create the empty residue_table dictionary
        
        
        # read in your residue file
        with open(self.filename, 'rb') as f:
            reader = csv.reader(f)
            
            for row in reader:

                # NOTE there's no validation here so make sure this file is correctly formatte
                res = Residue(row[0],row[1],row[2],float(row[3]),int(row[4]))
                
                # residue_table dictionary is keyed by the 3 letter code
                self.residue_table[res.letterCode3] = res
        """
    def lookForRes(self,resCode):
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
        raise ResTableException("Invalid amino acid code provided [" + str(resCode) +"]")
            

    def lookUpHydropathy(self,resCode):
        res = self.lookForRes(resCode)
        return res.hydropathy
        
    def lookUpCharge(self,resCode):
        
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


