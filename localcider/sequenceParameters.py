""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.0                                                         !
   !--------------------------------------------------------------------------!
   
   File Description:
   ================
   
   This is one of the key user interface/API files from which users (AKA you!)
   should use to access localCIDER.

   For a full description please see the documentation!

"""

from backend.sequence import Sequence
from backend.seqfileparser import SequenceFileParser 
from backend.backendtools import status_message
from backend import plotting

class SequenceParameters:

    def __init__(self, sequence="", sequenceFile=""):

        # provide the flexibility to submit either a sequence
        # file or an actual sequence as a string
        if sequence=="" and sequenceFile=="":
            return None

        # if the sequence isn't empty constuct a local 
        # sequence object using the sequence
        if not sequence == "":
            self.SeqObj = Sequence(sequence)
        else:
            parserMachine = SequenceFileParser()
            self.SeqObj = Sequence(parserMachine.parseSeqFile(sequenceFile))

    # ============================================ #
    # ============= SETTER FUNCTIONS ============= #  
    #...................................................................................#
    def set_phosphosites(self, phosphosites):
        """ 
        Set putativly phosphorylated sites on your sequence.
        You can call this multiple times to update, though it
        will only ever add sites to the list. To clear sites
        use the clear_phosphosites() method.
            
        INPUT: list of site positions, indexed from 1

        OUTPUT: None, but the underlying Sequence object
        has its properties updated
        """
        
        # note that we do all the relevant data validation
        # in the calling function
        self.SeqObj.setPhosPhoSites(phosphosites)


    #...................................................................................#
    def clear_phosphosites(self):
        """
        Clears the putative phosphosites on the underlying'
        sequence object. Useful if you want to reset the phosphostatus.
        
        OUTPUT: None, but the underlying Sequence object
        has its properties updated
        """

        self.SeqObj.clear_phosphosites()


    # ============================================ #
    # ============= GETTER FUNCTIONS ============= #  
    #...................................................................................#
    def get_sequence(self):
        """
        Get the protein's primary amino acid sequence

        OUTPUT: Single string with the protein sequence
        """
        
        return self.SeqObj.seq


    #...................................................................................#
    def get_mean_hydropathy(self):
        """
        Get a proteins mean hydropphobiicity

        OUTPUT: Float with the sequence's mean hydrophobicity
        """

        return self.SeqObj.meanHydropathy()


    #...................................................................................#
    def get_uversky_hydrophobicity(self):
        """
        Get a proteins mean hydrophobicity as defined by Uversky,
        using a normalized Kyte-Doolittle hydrophobicity table
        """
        
        return self.SeqObj.uverskyHydropathy()
        
        
    #...................................................................................#
    def get_fraction_disorder_promoting(self):
        """
        Get a proteins D to O ratio (ratio of disorder promiting residues to order promoting residues)
        """
        
        return self.SeqObj.fraction_disorder_promoting()
        

    #...................................................................................#
    def get_amino_acid_fractions(self):
        """
        Returns a dictionary with the fractions of each amino acid in your sequence
        """

        return self.SeqObj.amino_acid_fraction()


    #...................................................................................#
    def get_kappa(self):
        """ 
        Get the kappa value for a sequence

        OUTPUT: Float with the sequence's kappa
        """

        return self.SeqObj.kappa()


    #...................................................................................#
    def get_countPos(self):
        """ 
        Get the number of positive residues in the sequence 
        
        OUTPUT: Integer with number of positive residues in your sequence
        """

        return self.SeqObj.countPos()


    #...................................................................................#
    def get_countNeg(self):
        """ 
        Get the number of negative residues in the sequence

        OUTPUT: Integer with number of negative residues in your sequence
        """

        return self.SeqObj.countNeg() 


    #...................................................................................#
    def get_countNeut(self):
        """ 
        Get the number of neutral residues in the sequence

        OUTPUT: Integer with number of neutral residues in your sequence
        """
        return self.SeqObj.countNeut() 


    #...................................................................................#
    def get_fraction_positive(self):
        """ 
        Get the fraction of positive residues in the sequence 

        OUTPUT: Float with the sequence's fraction of positive residues (F+)
        """

        return self.SeqObj.Fplus() 


    #...................................................................................#
    def get_fraction_negative(self):
        """ 
        Get the fraction of negative residues in the sequence 

        OUTPUT: Float with the sequence's fraction of positive residues (F+)
        """

        return self.SeqObj.Fminus() 


    #...................................................................................#
    def get_FCR(self):      
        """ 
        Get the fraction of charged residues in the sequence 

        OUTPUT: Float with the sequence's fraction of charged residues
        """

        return self.SeqObj.FCR() 


    #...................................................................................#
    def get_NCPR(self):
        """ 
        Get the net charge per residue of the sequence 

        OUTPUT: Float with the sequence's net charge per residue
        """

        return self.SeqObj.NCPR() 

    #...................................................................................#
    def get_mean_net_charge(self):
        """ 
        Get the absolute magnitude of the mean net charge 


        OUTOUT: Float equal to the [absolute magnitude] of the mean net charge of the sequence
        """
        return self.SeqObj.mean_net_charge()


    #...................................................................................#
    def get_phasePlotRegion(self):
        """ 

        Returns the IDP diagram of states [REF 1] region based on the FCR/NCPR 
        rules. Possible return values are;

        

            1 +-----------------------+---------------+
              |                     +                 |
              |                   +                   |
              |      4          +                     |
        F     |               +                       |
        R     |             +                         | 
        A     |           +                           +
        C     |         +                           + |
        T     |       +                           +   |
        .     |     +            3              +     |
          0.5 |   +                           +       |
        N     | +                           +         |
        E     +                           +           |
        G     | +                       +             |
              + 2 +                   +               |
        R     | +   +               +       5         |
        E     |   + 2 +           +                   |
        S     |     +   +       +                     |
              | 1     + 2 +   +                       |
           0  +---------+---+-------------------------+

              0                  0.5                  1
                    FRACTION OF POSITIVE RESIDUES
                    
                    
        1) Weak polyampholytes and polyelectrolytes
        2) Janus sequences
        3) Strong polyampholytes
        4 and 5) Strong polyelectrolytes

        OUTPUT: Returns an integer describing the region on the density 
        of states diagmra (above)

        """

        return self.SeqObj.phasePlotRegion()


    #...................................................................................#
    def get_phosphosites(self):
        """
        Function which returns a list of currently assigned
        phosphosites on your sequence.

        OUTPUT: Returns a list of integers where each
        integer represents a site currently defined as 
        phosphorylatable. Such sites *must*, by definition
        by T/Y/S.

        OUTPUT: returns a list of integers corresponding to the sites
        which are currently defined as being phosphorylatable based on
        user input
        """
        
        return self.SeqObj.get_phosphosites()
        
        
    #...................................................................................#
    def get_kappa_after_phosphorylation(self):
        """
        Function which recomputes kappa after complete
        phosphorylation based on the currently defined
        phosphosites.

        OUTPUT: returns a float corresponding to the sequence's kappa
        value if all the currently defined phosphosites were phosphorylated
        """
        
        if len(self.get_phosphosites()) == 0:
            status_message("Be aware that there are no phosphosites currently set - getting 'naked' kappa")
        return self.SeqObj.kappa_at_maxPhos()


    #...................................................................................#
    def get_all_phosphorylatable_sites(self):
        """
        Function which returns a list of all the positions which *could* be
        phosphorylated (i.e. are T/S/Y). NOTE this does not use any kind of 
        smart lookup, metadata, or analysis. It's literally, where are the Y/T/S
        residues.

        Note positions are returned as indexed from 1 (so you can feed these positions
        directly into the set_phosphosites function.
        
        OUTPUT: Returns a list of integers corresponding to S/T/Y positions
        in your sequence
        """

        return self.SeqObj.get_STY_residues()

    
    #...................................................................................#
    def get_full_phosphostatus_kappa_distribution(self):
        """
        This function calculates the kappa value of all possible phosphorylation
        statuses, given the defined phosphosites. 

        """
        
        ncalcs = self.SeqObj.calculateNumberDifferentPhosphoStates()
        
        status_message("This function will now make " + str(ncalcs) + " independent kappa calculations\nIf this is a big number you may want to investigate a subset of possible phosphosites or\nuse a Monte Carlo approach to subsample")
        
        return self.SeqObj.calculateKappaDistOfPhosphoStates()


    #...................................................................................#
    def get_phosphosequence(self):
        """
        Returns the sequence with phosphorylated sites set to E instead of S/Y/T
        """

        return self.SeqObj.get_phosphosequence()

        
    
    # ============================================ #
    # ======= PLOTTING DIAGRAM FUNCTIONS ========= #  
    #...................................................................................#
    def show_phaseDiagramPlot(self,label=""):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and creates it on the screen

        <<REQUIRES MATPLOTLIB>>>
       
        """

        plotting.show_single_phasePlot(self.get_fraction_positive(), self.get_fraction_negative(),label)


    #...................................................................................#
    def save_phaseDiagramPlot(self, filename,label=""):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and saves it at the <filename> location

        <<REQUIRES MATPLOTLIB>>>

        INPUT: Writeable filename
        OUTPUT: Nothing, but creates a .png file at the filename location

        """
        
        plotting.save_single_phasePlot(self.get_fraction_positive(), self.get_fraction_negative(), filename, label)


    #...................................................................................#
    def show_uverskyPlot(self,label=""):
        """ 
        Generates the Uversky phase diagram (hydropathy vs NCPR), places
        this sequence on that plot, and creates it on the screen

        <<REQUIRES MATPLOTLIB>>>
       
        """

        plotting.show_single_uverskyPlot(self.get_uversky_hydrophobicity(), self.get_mean_net_charge(), label)


    #...................................................................................#
    def save_uverskyPlot(self, filename,label=""):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and saves it at the <filename> location

        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but creates a .png file at the filename location

        """
        
        plotting.save_single_uverskyPlot(self.get_uversky_hydrophobicity(), self.get_mean_net_charge(), filename, label)


    

    #...................................................................................#
    def save_linearNCPR(self, blobLen, filename):
        """ 
        Generates a plot of how the NCPR (net charge per residue) changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.
        
        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but creates a .png file at the filename location

        """

        plotting.save_linearplot(plotting.build_NCPR_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def save_linearSigma(self, blobLen, filename):
        """ 
        Generates a plot of how the Sigma parameter changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Recall that sigma is defined as

        NCPR^2/FCR (net charge per residue squared divided by the fraction of charged residues)

        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but creates a .png file at the filename location

        """

        plotting.save_linearplot(plotting.build_sigma_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def save_linearHydropathy(self, blobLen, filename):
        """ 
        Generates a plot of how the mean hydropathy changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Hydropathy here is calculated using a NORMALIZED Kyte-Doolittle scale, where 1 is
        the most hydrophobic and 0 the least.

        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but creates a .png file at the filename location

        """

        plotting.save_linearplots(plotting.build_hydropathy_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def show_linearNCPR(self, blobLen):
        """ 
        Generates a plot of how the NCPR (net charge per residue) changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.
        
        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but the plot is displayed on screen

        """

        plotting.show_linearplot(plotting.build_NCPR_plot, self.SeqObj, blobLen)


    #...................................................................................#
    def show_linearSigma(self, blobLen):
        """ 
        Generates a plot of how the sigma parameter changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Recall that sigma is defined as

        NCPR^2/FCR (net charge per residue squared divided by the fraction of charged residues)

        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but the plot is displayed on screen

        """

        plotting.show_linearplot(plotting.build_sigma_plot, self.SeqObj, blobLen)


    #...................................................................................#
    def show_linearHydropathy(self, blobLen):
        """ 
        Generates a plot of how the mean hydropathy changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Hydropathy here is calculated using a NORMALIZED Kyte-Doolittle scale, where 1 is
        the most hydrophobic and 0 the least.

        <<REQUIRES MATPLOTLIB>>>

        INPUT:  Writeable filename
        OUTPUT: Nothing, but the plot is displayed on screen
        """

        plotting.show_linearplot(plotting.build_hydropathy_plot, self.SeqObj, blobLen)




    # ============================================ #
    # ============ IMPLICIT FUNCTIONS ============ #  
        
    #...................................................................................#
    def __unicode__(self):
        """ Returns the sequences """
        return "SequenceParameter [len="+str(len(self.SeqObj.seq)) + "]" 

        
    #...................................................................................#
    def __str__(self):
        """ Returns the sequences """
        return self.__unicode__(self)


    
        
        
            



