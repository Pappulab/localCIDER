""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Copyright (C) 2014 - 2021                                             !
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

TWENTY_AAs = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','I','L','M','F','W','Y','V']


# Update on 0.1.11 gly went from black to green
# Aromatics became orange
#
#
DEFAULT_COLOR_PALETTE = {'A':'black', 
                         'C':'black',
                         'D':'red',
                         'E':'red',
                         'F':'orange',
                         'G':'green',
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
                         'W':'orange',
                         'Y':'orange'}






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

def get_WW_original():
    """ Function which returns the original Wimley and White hydropathy lookup table
    """
    return {'ILE':  0.31,
            'VAL': -0.07,
            'LEU':  0.56,
            'PHE':  1.13,
            'CYS':  0.24,
            'MET':  0.23,
            'ALA': -0.17,
            'GLY': -0.01,
            'THR': -0.14,
            'SER': -0.13,
            'TRP':  1.85,
            'TYR':  0.94,
            'PRO': -0.45,
            'HIS': -0.96,
            'GLU': -2.02,
            'GLN': -0.58,
            'ASP': -1.23,
            'ASN': -0.42,
            'LYS': -0.99,
            'ARG': -0.81}
           
def get_PPII_Hilser():
    """
    Returns an amino acid dictionary with the PPII propensity of
    each residue as calculated by Elam et al [1] (Taken from [2]). 

    [1] - Elam WA, Schrank TP, Campagnolo AJ, Hilser VJ. Evolutionary 
    conservation of the polyproline II conformation surrounding intrinsically 
    disordered phosphorylation sites. 
    Protein Sci. 2013; 22: 405 - 417. doi: 10.1002/pro.2217 PMID: 23341186

    [2] - Tomasso, M. E., Tarver, M. J., Devarajan, D. & Whitten, S. T. 
    Hydrodynamic Radii of Intrinsically Disordered Proteins Determined 
    from Experimental Polyproline II Propensities. 
    PLoS Comput. Biol. 12, e1004686 (2016).

    """
    return  {'ILE': 0.39,                
             'VAL': 0.39,
             'LEU': 0.24,
             'PHE': 0.17,
             'CYS': 0.25,
             'MET': 0.36,
             'ALA': 0.37,
             'GLY': 0.13,
             'THR': 0.32,
             'SER': 0.24,
             'TRP': 0.25,
             'TYR': 0.25,
             'PRO': 1.00,
             'HIS': 0.20,
             'GLU': 0.42,
             'GLN': 0.53,
             'ASP': 0.30,
             'ASN': 0.27,
             'LYS': 0.56,
             'ARG': 0.38}

def get_PPII_Creamer():
    """
    Returns an amino acid dictionary with the PPII propensity of
    each residue as calculated by Rucker et al [1] (Taken from [2]). 

    Note that Trp and Try do not have values reported by Rucker 
    et al., so we followed the convention uesed by Tomasso et al. 
    and set both to the mean value (0.58)
    

    [1] - Rucker, A.L., Pager, C.T., Campbell, M.N., Qualls, J.E., 
    and Creamer, T.P. (2003). Host-guest scale of left-handed 
    polyproline II helix formation. Proteins 53, 68-75. 

    [2] - Tomasso, M. E., Tarver, M. J., Devarajan, D. & Whitten, S. T. 
    Hydrodynamic Radii of Intrinsically Disordered Proteins Determined 
    from Experimental Polyproline II Propensities. 
    PLoS Comput. Biol. 12, e1004686 (2016).

    """
    return  {'ILE': 0.50,
             'VAL': 0.49,
             'LEU': 0.58,
             'PHE': 0.58,
             'CYS': 0.55,
             'MET': 0.55,
             'ALA': 0.61, 
             'GLY': 0.58,
             'THR': 0.53,
             'SER': 0.58,
             'TRP': 0.58,
             'TYR': 0.58,
             'PRO': 0.67,
             'HIS': 0.55,
             'GLU': 0.61,
             'GLN': 0.66,
             'ASP': 0.63,
             'ASN': 0.55,
             'LYS': 0.59,
             'ARG': 0.61}


def get_PPII_Kallenbach():
    """
    Returns an amino acid dictionary with the PPII propensity of
    each residue as calculated by Shi et al [1] (Taken from [2]).

    Note that Gly and Pro do not have values reported by Shi et al.,
    so we followed the convention uesed by Tomasso et al. and 
    set Gly = 0.5 and Pro = 1.0.

    [1] - Shi, Z., Chen, K., Liu, Z., Ng, A., Bracken, W.C., and 
    Kallenbach, N.R. (2005). Polyproline II propensities from 
    GGXGG peptides reveal an anticorrelation with beta-sheet scales. 
    Proc. Natl. Acad. Sci. U. S. A. 102, 17964-17968.

    [2] - Tomasso, M. E., Tarver, M. J., Devarajan, D. & Whitten, S. T. 
    Hydrodynamic Radii of Intrinsically Disordered Proteins Determined 
    from Experimental Polyproline II Propensities. 
    PLoS Comput. Biol. 12, e1004686 (2016).

    """
    return  {'ILE': 0.519,
             'VAL': 0.743,
             'LEU': 0.574,
             'PHE': 0.639,
             'CYS': 0.557,
             'MET': 0.498, 
             'ALA': 0.818,
             'GLY': 0.500,
             'THR': 0.553,
             'SER': 0.774,
             'TRP': 0.764,
             'TYR': 0.630,
             'PRO': 1.000,
             'HIS': 0.428,
             'GLU': 0.684,
             'GLN': 0.654,
             'ASP': 0.552,
             'ASN': 0.667,
             'LYS': 0.581,
             'ARG': 0.638}


def get_pKa():
    """
    Return pKa values for titratable residues. Values taken here correspond to values from
    EMBOSS, but in reality pKa values are pretty context-dependent so these are 'model'
    values, but the true value in the context of a protein is likely to be pretty seriously
    different.

    """
    return {'C': 8.5,
            'Y': 10.1,
            'H': 6.5,
            'E': 4.1,
            'D': 3.9,
            'K': 10.0,
            'R': 12.5}
    

def get_molecular_weight_Da():
    """
    Returns a lookup table with the molecular weights of each
    amino acid in Da

    """
    return  {'I': 131.2,
             'V': 117.1,
             'L': 131.2,
             'F': 165.2,
             'C': 121.2,
             'M': 149.2,
             'A': 89.1,
             'G': 75.1,
             'T': 119.1,
             'S': 105.1,
             'W': 204.2,
             'Y': 181.2,
             'P': 115.1,
             'H': 155.2,
             'E': 147.1,
             'Q': 146.2,
             'D': 133.1,
             'N': 132.1,
             'K': 146.2,
             'R': 174.2,
             }



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

    PPII_Dict_Hilser       = get_PPII_Hilser()
    PPII_Dict_Creamer      = get_PPII_Creamer()
    PPII_Dict_Kallenbach   = get_PPII_Kallenbach()
    
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
        
        # update the residue with the PPII content
        res.append(PPII_Dict_Hilser[res[1]])
        res.append(PPII_Dict_Creamer[res[1]])
        res.append(PPII_Dict_Kallenbach[res[1]])

    # each residue is now defined by a list of length 8 which are
    # [full name, 3 letter, one letter, Hydrophobicity, charge, PPII_Hilser, PPII_Creamer, PPII_Kallenbach]
    
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
