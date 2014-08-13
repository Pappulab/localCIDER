""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.1                                                         !
   !                                                                          !
   !    Copyright (C) 2014, The localCIDER development team (current and      !
   !                        former contributors): Alex Holehouse, James       !
   !                        Ahad, Rahul K. Das.                               !
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
   ! MAIN AUTHOR:    Alex Holehouse                                           !
   !                                                                          !
   !--------------------------------------------------------------------------!

   Amino acids data

"""


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <>          
# <>          THREE-TO-ONE LETTER CODE TRANSLATION
# <> 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


THREE_TO_ONE = {'ALA':'A', 
                'CYS':'C',
                'ASP':'D',
                'GLU':'E',
                'PHE':'F',
                'GLY':'G',
                'HIS':'H', 
                'ILE':'I',
                'LYS':'K',
                'LEU':'L',
                'MET':'M',
                'ASN':'N',
                'PRO':'P',
                'GLN':'Q',
                'ARG':'R',
                'SER':'S',
                'THR':'T',
                'VAL':'V',
                'TRP':'W',
                'TYR':'Y'}

ONE_TO_THREE = {'A':'ALA', 
                'C':'CYS',
                'D':'ASP',
                'E':'GLU',
                'F':'PHE',
                'G':'GLY',
                'H':'HIS', 
                'I':'ILE',
                'K':'LYS',
                'L':'LEU',
                'M':'MET',
                'N':'ASN',
                'P':'PRO',
                'Q':'GLN',
                'R':'ARG',
                'S':'SER',
                'T':'THR',
                'V':'VAL',
                'W':'TRP',
                'Y':'TYR'}





# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# <>          
# <>          HYDROPHOBICITY DEFINING FUNCTIONS
# <> 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 
# KYTE-DOOLITTLE SCALES
#
# References
#
# Main hydrophobicity scale
# ...............................
#
# A simple method for displaying the hydropathic character of a protein. 
# Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
#
#
#
# 0 to 1 normalized hydrophobicity scale
# ...............................................
#
# Why are "natively unfolded" proteins unstructured under physiological conditions?
# Valdimir N. Uversky, Joel R. Gillespie, and Anthony L. Frink
# Protines: Structure, function, and genetics 41:415-427 (2000)
#

def get_KD_original():
    """ Function which returns the original KD hydropathy lookup table
    """    
    
    return  {'ILE': 4.5,                
             'VAL': 4.2,
             'LEU': 3.8,
             'PHE': 2.8,
             'CYS': 2.5,
             'MET': 1.9,
             'ALA': 1.8,
             'GLY': -0.4,
             'THR': -0.7,
             'SER': -0.8,
             'TRP': -0.9,
             'TYR': -1.3,
             'PRO': -1.6,
             'HIS': -3.2,
             'GLU': -3.5,
             'GLN': -3.5,
             'ASP': -3.5,
             'ASN': -3.5,
             'LYS': -3.9,
             'ARG': -4.5}

def get_KD_shifted():
    """ 
    Function which returns the shifted KD hydropathy lookup table (such that
    it runs from 0 to 9 instead of -4.5 to 4.5)
    """
    
    original  = get_KD_original()
    
    shifted = {}
    for i in original:
        shifted[i] = original[i]+4.5

    return shifted
        

    """ Should look like this!
    {'ALA': 6.3,
    'ARG': 0.0,
    'ASN': 1.0,
    'ASP': 1.0,
    'CYS': 7.0,
    'GLN': 1.0,
    'GLU': 1.0,
    'GLY': 4.1,
    'HIS': 1.3,
    'ILE': 9.0,
    'LEU': 8.3,
    'LYS': 0.6,
    'MET': 6.4,
    'PHE': 7.3,
    'PRO': 2.9,
    'SER': 3.7,
    'THR': 3.8,
    'TRP': 3.6,
    'TYR': 3.2,
    'VAL': 8.7
    }
    """

def get_KD_uversky():
    """
    Returns a 0-to-1 normalized KD scale

    """
    shifted = get_KD_shifted()
    
    uversky = {}
    for i in shifted:
        uversky[i] = shifted[i]/9.0

    return uversky 

           


""""
EMPTY       = {'ILE': 
               'VAL': 
               'LEU': 
               'PHE': 
               'CYS': 
               'MET': 
               'ALA': 
               'GLY': 
               'THR': 
               'SER': 
               'TRP': 
               'TYR': 
               'PRO': 
               'HIS': 
               'GLU': 
               'GLN': 
               'ASP': 
               'ASN': 
               'LYS': 
               'ARG': }
"""


def buildTable():
    return build_amino_acids_skeleton()


def build_amino_acids_skeleton():
    """ 
    The default hydrophobicity is an augmented Kyte-Doolitle where
    0 = most hydrophilic and 9 is most hydrophobic

    """
    
    return [["Alanine",        "ALA", "A", 0,  0],
            ["Cysteine",       "CYS", "C", 7.0,  0],
            ["Aspartic_Acid",  "ASP", "D", 1.0, -1],
            ["Glutamic_Acid",  "GLU", "E", 1.0, -1],
            ["Phenylalanine",  "PHE", "F", 7.3,  0],
            ["Glycine",        "GLY", "G", 4.1,  0],
            ["Histidine",      "HIS", "H", 1.3,  0],
            ["Isoleucine",     "ILE", "I", 9.0,  0],
            ["Lysine",         "LYS", "K", 0.6,  1],
            ["Leucine",        "LEU", "L", 8.3,  0],
            ["Methionine",     "MET", "M", 6.4,  0],
            ["Asparagine",     "ASN", "N", 1.0,  0],
            ["Proline",        "PRO", "P", 2.9,  0],
            ["Glutamine",      "GLN", "Q", 1.0,  0],
            ["Arginine",       "ARG", "R", 0.0,  1],
            ["Serine",         "SER", "S", 3.7,  0],
            ["Threonine",      "THR", "T", 3.8,  0],
            ["Valine",         "VAL", "V", 8.7,  0],
            ["Tryptophan",     "TRP", "W", 3.6,  0],
            ["Tyrosine",       "TYR", "Y", 3.2,  0]]
    
    

def update_hydrophobicity(aalist, scale):
    
    
    index=0
    for i in aalist:
        i[3] = scale[i[1]]

    return aalist
        
        
