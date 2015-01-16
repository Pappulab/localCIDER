""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.4                                                         !
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

DEFAULT_COLOR_PALETTE = {'A':'black', 
                         'C':'black',
                         'D':'red',
                         'E':'red',
                         'F':'black',
                         'G':'black',
                         'H':'green', 
                         'I':'black',
                         'K':'blue',
                         'L':'black',
                         'M':'black',
                         'N':'green',
                         'P':'fuchsia',
                         'Q':'green',
                         'R':'blue',
                         'S':'green',
                         'T':'green',
                         'V':'black',
                         'W':'black',
                         'Y':'black'}






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

def get_residue_charge():
    """ Function which returns the original KD hydropathy lookup table
    """    
    
    return  {'ILE': 0,                
             'VAL': 0,
             'LEU': 0,
             'PHE': 0,
             'CYS': 0,
             'MET': 0,
             'ALA': 0,
             'GLY': 0,
             'THR': 0,
             'SER': 0,
             'TRP': 0,
             'TYR': 0,
             'PRO': 0,
             'HIS': 0,
             'GLU': -1,
             'GLN': 0,
             'ASP': -1,
             'ASN': 0,
             'LYS': 1,
             'ARG': 1}



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

    This is a weird way of doing this but essentially it means we define
    the charge and hydrophobicity in one place (in separate functions)
    and then construct the skeleton before appending those values to the
    amino acid skeleton.

    The advantage of this is we can add further information to induvidual 
    amino acids. 

    """
    
    # get a dictionary of 3 letter to charge
    charge_Dict = get_residue_charge()

    # get a dictionary of 3 letter to KD hydrophobicity
    KD_Dict     = get_KD_shifted()
    
    # build the initial skeleton of amnino acids 
    skeleton=[["Alanine",        "ALA", "A"],
            ["Cysteine",       "CYS", "C"],
            ["Aspartic_Acid",  "ASP", "D"],
            ["Glutamic_Acid",  "GLU", "E"],
            ["Phenylalanine",  "PHE", "F"],
            ["Glycine",        "GLY", "G"],
            ["Histidine",      "HIS", "H"],
            ["Isoleucine",     "ILE", "I"],
            ["Lysine",         "LYS", "K"],
            ["Leucine",        "LEU", "L"],
            ["Methionine",     "MET", "M"],
            ["Asparagine",     "ASN", "N"],
            ["Proline",        "PRO", "P"],
            ["Glutamine",      "GLN", "Q"],
            ["Arginine",       "ARG", "R"],
            ["Serine",         "SER", "S"],
            ["Threonine",      "THR", "T"],
            ["Valine",         "VAL", "V"],
            ["Tryptophan",     "TRP", "W"],
            ["Tyrosine",       "TYR", "Y"]]


    # for each residue update with further information
    for res in skeleton:

        # update this residue with the charge value
        res.append(KD_Dict[res[1]])

        # update this residue with the KD hydrophobicity
        res.append(charge_Dict[res[1]])
    
    return skeleton
    

def update_hydrophobicity(aalist, scale):
    """
    Function which lets you update the hydrophobicity of an amino
    acid to your own scale.

    """
    
    index=0
    for i in aalist:
        i[3] = scale[i[1]]

    return aalist
