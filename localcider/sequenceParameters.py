""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.3                                                         !
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
from backend.localciderExceptions import SequenceException

class SequenceParameters:

    def __init__(self, sequence="", sequenceFile=""):

        # provide the flexibility to submit either a sequence
        # file or an actual sequence as a string        
        if sequence=="" and sequenceFile=="":
            raise SequenceException("Empty sequence/sequence file")

        # if the sequence isn't empty constuct a local 
        # sequence object using the sequence
        if not sequence == "":
            
            # note we set the validate flag to check the sequence is valid and deal
            # with whitespace
            self.SeqObj = Sequence(sequence,validateSeq=True)
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
            
        INPUT: 
        --------------------------------------------------------------------------------
        phosphosites    |  list of site positions, indexed from 1


        OUTPUT: 
        --------------------------------------------------------------------------------
        None, but the underlying Sequence object has its properties
        updated
        """
        
        # note that we do all the relevant data validation
        # in the calling function
        self.SeqObj.setPhosPhoSites(phosphosites)


    #...................................................................................#
    def clear_phosphosites(self):
        """
        Clears the putative phosphosites on the underlying
        sequence object. Useful if you want to reset the phosphostatus.
        
        OUTPUT: 
        --------------------------------------------------------------------------------
        None, but the underlying Sequence object has its properties updated

        """

        self.SeqObj.clear_phosphosites()


    # ============================================ #
    # ============= GETTER FUNCTIONS ============= #  
    #...................................................................................#
    def get_sequence(self):
        """
        Get the protein's primary amino acid sequence

        OUTPUT: 
        --------------------------------------------------------------------------------
        Single string with the protein sequence
        """
        
        return self.SeqObj.seq


    #...................................................................................#
    def get_length(self):
        """
        Get the protein's amino acid sequence length

        OUTPUT: 
        --------------------------------------------------------------------------------
        Integer equal to the sequence length
        """
        
        return len(self.SeqObj.seq)


    #...................................................................................#
    def get_mean_hydropathy(self):
        """
        Get a protein's mean hydropathy (a value between 0 and 9, where 0 is the least 
        hydrophobic and 9 is the most). This is simply a re-based Kyte-Doolittle scale,
        which runs from 0 to 9 instead of from -4.5 to 4.5, as in the original paper.

        OUTPUT: 
        --------------------------------------------------------------------------------
        Float with the sequence's mean hydropathy
        """

        return self.SeqObj.meanHydropathy()


    #...................................................................................#
    def get_uversky_hydropathy(self):
        """
        Get a protein's mean hydropathy as defined by Uversky,
        using a normalized Kyte-Doolittle hydropathy table. Here, values range between
        0 and 1, with 0 being the least hydrophobic and 1 the most hydrophobic.

        OUTPUT:
        --------------------------------------------------------------------------------
        Float with the normalized Kyte-Doolittle hydropathy (between 0 and 1)
        """
        
        return self.SeqObj.uverskyHydropathy()
        
        
    #...................................................................................#
    def get_fraction_disorder_promoting(self):
        """
        Get a protein's fraction of residues which are considered '[D]isorder promoting'.
        
        For more details see the reference below;
        
        ********************************************************************************
        Reference:
        TOP-IDP-scale: a new amino acid scale measuring propensity for intrinsic disorder.
        Protein Pept Lett. 2008;15(9):956-63.
        Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK.
        ********************************************************************************
        
        
        OUTPUT:
        --------------------------------------------------------------------------------
        Float with the fraction of disorder promoting residues
        
        """
        
        return self.SeqObj.fraction_disorder_promoting()
        

    #...................................................................................#
    def get_amino_acid_fractions(self):
        """
        Returns a dictionary with the fractions of each amino acid in your sequence
        
        OUTPUT:
        --------------------------------------------------------------------------------
        Dictionary of amino acids, where the keys are each of the 20 amino acids and the
        values represents the fraction of that amino acid
        
        """

        return self.SeqObj.amino_acid_fraction()


    #...................................................................................#
    def get_kappa(self):
        """ 
        Get the kappa value for a sequence.  

        OUTPUT: 
        --------------------------------------------------------------------------------
        Float with the sequence's kappa value
        """

        return self.SeqObj.kappa()

    
    #...................................................................................#
    def get_deltaMax(self):
        """
        Get the maximum delta value for a sequence of this composition. Note kappa is 
        delta/deltaMax.

        OUTPUT: 
        --------------------------------------------------------------------------------
        Float with the sequence's delta max (identical for all permutants)
        """
        
        return self.SeqObj.deltaMax()


    #...................................................................................#
    def get_delta(self):
        """
        Get the delta value for this specific sequence. Note kappa is delta/deltaMax.

        OUTPUT:
        --------------------------------------------------------------------------------
        Float with the sequence's delta max (will vary with permutants)
        """
        
        return self.SeqObj.delta()
        

    #...................................................................................#
    def get_countPos(self):
        """ 
        Get the number of positive residues in the sequence 
        
        OUTPUT:
        --------------------------------------------------------------------------------
        Integer with number of positive residues in your sequence
        """

        return self.SeqObj.countPos()


    #...................................................................................#
    def get_countNeg(self):
        """ 
        Get the number of negative residues in the sequence

        OUTPUT:
        --------------------------------------------------------------------------------
        Integer with number of negative residues in your sequence
        """

        return self.SeqObj.countNeg() 


    #...................................................................................#
    def get_countNeut(self):
        """ 
        Get the number of neutral residues in the sequence

        OUTPUT: 
        --------------------------------------------------------------------------------
        Integer with number of neutral residues in your sequence
        """
        return self.SeqObj.countNeut() 


    #...................................................................................#
    def get_fraction_positive(self):
        """ 
        Get the fraction of positive residues in the sequence 

        OUTPUT:
        --------------------------------------------------------------------------------
        Float with the sequence's fraction of positive residues (F+)
        """

        return self.SeqObj.Fplus() 


    #...................................................................................#
    def get_fraction_negative(self):
        """ 
        Get the fraction of negative residues in the sequence 

        OUTPUT
        --------------------------------------------------------------------------------: 
        Float with the sequence's fraction of positive residues (F+)
        """

        return self.SeqObj.Fminus() 


    #...................................................................................#
    def get_FCR(self):      
        """ 
        Get the fraction of charged residues in the sequence 

        OUTPUT: 
        --------------------------------------------------------------------------------
        Float with the sequence's fraction of charged residues
        """

        return self.SeqObj.FCR() 


    #...................................................................................#
    def get_NCPR(self):
        """ 
        Get the net charge per residue of the sequence 

        OUTPUT: 
        --------------------------------------------------------------------------------
        Float with the sequence's net charge per residue
        """

        return self.SeqObj.NCPR() 

    #...................................................................................#
    def get_mean_net_charge(self):
        """ 
        Get the absolute magnitude of the mean net charge 

        OUTOUT:
        --------------------------------------------------------------------------------
        Float equal to the [absolute magnitude] of the mean net charge of the sequence
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
        R     | +   +               +       4         |
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

        OUTPUT: 
        --------------------------------------------------------------------------------
        Returns an integer describing the region on the density of states 
        diagram (above)

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

        OUTPUT: 
        --------------------------------------------------------------------------------
        Returns a list of integers corresponding to the sites
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

        OUTPUT: 
        --------------------------------------------------------------------------------
        returns a float corresponding to the sequence's kappa value if 
        all the currently defined phosphosites were phosphorylated
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
        
        OUTPUT: 
        --------------------------------------------------------------------------------
        Returns a list of integers corresponding to S/T/Y positions in your sequence
        """

        return self.SeqObj.get_STY_residues()

    
    #...................................................................................#
    def get_full_phosphostatus_kappa_distribution(self):
        """
        This function calculates the kappa value of all possible phosphorylation
        statuses, given the defined phosphosites. 

        OUTPUT:
        --------------------------------------------------------------------------------
        
        

        """
        
        ncalcs = self.SeqObj.calculateNumberDifferentPhosphoStates()
        status_message("Running exaustive kappa distribution analysis based on phosphorylation states")
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
    def show_phaseDiagramPlot(self,label="", title="Diagram of states",legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and creates it on the screen

        INPUT: 
        --------------------------------------------------------------------------------
        label     | A label for the point on the phase diagram
        title     | Plot title (DEFAULT = 'Diagram of states')
        legendOn  | Boolean for if the figure legend should be displayed or not
        xLim      | Max value for the x axis (fract. positive charge) (DEFAULT = 1)
        yLim      | Max value for the y axis (fract. negative charge) (DEFAULT = 1)
        fontSize  | Size of font for label (DEFAULT = 10)
        getFig    | Returns a matplotlib figure object instead of simply displaying the 
                  | plot on the screen (DEFAULT = False)


        OUTPUT: 
        --------------------------------------------------------------------------------
        If the argument getFig is False (which it is by default) then the Uversky plot appears on
        the screen. If getFig is set to True then the function returns a matplotlib plt object, which
        can be further manipulated.
        
        """

        if getFig:
            return plotting.show_single_phasePlot(self.get_fraction_positive(), self.get_fraction_negative(),label, legendOn, xLim, yLim, fontSize, getFig)
        else:
            plotting.show_single_phasePlot(self.get_fraction_positive(), self.get_fraction_negative(),label, title, legendOn, xLim, yLim, fontSize, getFig)


    #...................................................................................#
    def save_phaseDiagramPlot(self, filename,label="", title="Diagram of states", legendOn=True, xLim=1, yLim=1, fontSize=10):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and saves it at the <filename> location

        INPUT: 
        --------------------------------------------------------------------------------
        filename  | Writeable filename  
        
        label     | A label for the point on the phase diagram
        title     | Plot title (DEFAULT = 'Diagram of states')
        legendOn  | Boolean for if the figure legend should be displayed or not
        xLim      | Max value for the x axis (fract. positive charge) (DEFAULT = 1)
        yLim      | Max value for the y axis (fract. negative charge) (DEFAULT = 1)
        fontSize  | Size of font for label (DEFAULT = 10)


        OUTPUT:
        --------------------------------------------------------------------------------
        Nothing, but creates a .png file at the filename location

        """
        
        plotting.save_single_phasePlot(self.get_fraction_positive(), self.get_fraction_negative(), filename, label, title, legendOn, xLim, yLim, fontSize)


    #...................................................................................#
    def show_uverskyPlot(self, label="", title="Uversky plot", legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False):
        """ 
        Generates the Uversky phase diagram (hydropathy vs NCPR), places
        this sequence on that plot, and creates it on the screen

        INPUT:
        --------------------------------------------------------------------------------
        label     | A label for the point on the phase diagram
        title     | Plot title (DEFAULT = 'Uversky plot')
        legendOn  | Boolean for if the figure legend should be displayed or not
        xLim      | Max value for the x axis (mean net charge) (DEFAULT = 1)
        yLim      | Max value for the y axis (hydropathy) (DEFAULT = 1)
        fontSize  | Size of font for label (DEFAULT = 10)
        getFig    | Returns a matplotlib figure object instead of simply displaying the 
                    plot on the screen (DEFAULT = False)
        
        
        OUTPUT: 
        --------------------------------------------------------------------------------
        If the argument getFig is False (which it is by default) then the Uversky plot appears on
        the screen. If getFig is set to True then the function returns a matplotlib plt object, which
        can be further manipulated.
        
       
        """
        if getFig:
            return plotting.show_single_uverskyPlot(self.get_uversky_hydropathy(), self.get_mean_net_charge(), label, title, legendOn, xLim, yLim, fontSize, getFig)
        else:
            plotting.show_single_uverskyPlot(self.get_uversky_hydropathy(), self.get_mean_net_charge(), label, title, legendOn, xLim, yLim, fontSize, getFig)


    #...................................................................................#
    def save_uverskyPlot(self, filename, label="", title="Uversky plot", legendOn=True, xLim=1, yLim=1, fontSize=10):
        """ 
        Generates the Pappu-Das phase diagram (diagram of states), places
        this sequence on that plot, and saves it at the <filename> location

        INPUT:   
        --------------------------------------------------------------------------------
        filename  | A writeable filename

        label     | A label for the point on the phase diagram
        title     | Plot title (DEFAULT = 'Uversky plot')
        legendOn  | Boolean for if the figure legend should be displayed or not
        xLim      | Max value for the x axis (mean net charge) (DEFAULT = 1)
        yLim      | Max value for the y axis (hydropathy) (DEFAULT = 1)
        fontSize  | Size of font for label (DEFAULT = 10)
        
        
        OUTPUT: 
        --------------------------------------------------------------------------------
        If the argument getFig is False (which it is by default) then the Uversky plot appears on
        the screen. If getFig is set to True then the function returns a matplotlib plt object, which
        can be further manipulated.
        
        """
        
        plotting.save_single_uverskyPlot(self.get_uversky_hydropathy(), self.get_mean_net_charge(), filename, label, title, legendOn, xLim, yLim, fontSize)


    

    #...................................................................................#
    def save_linearNCPR(self, filename, blobLen=5):
        """ 
        Generates a plot of how the NCPR (net charge per residue) changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.
        
        INPUT:  
        --------------------------------------------------------------------------------
        filename | Name of the file to write
        bloblen  | Set the windowsize (DEFAULT = 5)
        

        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but creates a .png file at the filename location


        """

        plotting.save_linearplot(plotting.build_NCPR_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def save_linearSigma(self, filename, blobLen=5):
        """ 
        Generates a plot of how the Sigma parameter changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Recall that sigma is defined as

        NCPR^2/FCR (net charge per residue squared divided by the fraction of charged residues)

        INPUT:  
        --------------------------------------------------------------------------------
        filename | Name of the file to write
        bloblen  | Set the windowsize (DEFAULT = 5)


        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but creates a .png file at the filename location

        """

        plotting.save_linearplot(plotting.build_sigma_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def save_linearHydropathy(self, filename, blobLen=5):
        """ 
        Generates a plot of how the mean hydropathy changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Hydropathy here is calculated using a NORMALIZED Kyte-Doolittle scale, where 1 is
        the most hydrophobic and 0 the least.

        INPUT:  
        --------------------------------------------------------------------------------
        filename | Name of the file to write
        bloblen  | Set the windowsize (DEFAULT = 5)


        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but creates a .png file at the filename location

        """

        plotting.save_linearplot(plotting.build_hydropathy_plot, self.SeqObj, blobLen, filename)


    #...................................................................................#
    def show_linearNCPR(self, blobLen=5, getFig=False):
        """ 
        Generates a plot of how the NCPR (net charge per residue) changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        INPUT:  
        --------------------------------------------------------------------------------
        bloblen  | Set the windowsize (DEFAULT = 5)
        getFig   | Do you want to get the matplotlib figure object (DEFAULY = FALSE)


        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but the plot is displayed on screen

        """

        return plotting.show_linearplot(plotting.build_NCPR_plot, self.SeqObj, blobLen, getFig)


    #...................................................................................#
    def show_linearSigma(self, blobLen=5, getFig=False):
        """ 
        Generates a plot of how the sigma parameter changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Recall that sigma is defined as

        NCPR^2/FCR (net charge per residue squared divided by the fraction of charged residues)

        INPUT:  
        --------------------------------------------------------------------------------
        bloblen  | Set the windowsize (DEFAULT = 5)
        getFig   | Do you want to get the matplotlib figure object (DEFAULY = FALSE)


        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but the plot is displayed on screen

        """

        return plotting.show_linearplot(plotting.build_sigma_plot, self.SeqObj, blobLen, getFig)


    #...................................................................................#
    def show_linearHydropathy(self, blobLen=5, getFig=False):
        """ 
        Generates a plot of how the mean hydropathy changes as we move
        along the linear amino acid sequence in blobLen size steps. This uses a sliding window 
        approach and calculates the average within that window.

        Hydropathy here is calculated using a NORMALIZED Kyte-Doolittle scale, where 1 is
        the most hydrophobic and 0 the least.

        INPUT:  
        --------------------------------------------------------------------------------
        bloblen  | Set the windowsize (DEFAULT = 5)
        getFig   | Do you want to get the matplotlib figure object (DEFAULY = FALSE)


        OUTPUT: 
        --------------------------------------------------------------------------------
        Nothing, but the plot is displayed on screen

        """

        return plotting.show_linearplot(plotting.build_hydropathy_plot, self.SeqObj, blobLen, getFig)
        

    # ============================================== #
    # ============ FORMATTING FUNCTIONS ============ #  


    #...................................................................................#
    def set_HTMLColorResiduePalette(self, colorDict):
        """
        Allows the user to define the colormapping used for the get_HTMLColorString.

        The input parameter is a dictionary which contains amino acid single letter codes 
        as keys and colours as values. The colours must be one of the 17 standard HTML
        color names. These are;

        aqua, black, blue, fuchsia, gray, green, lime, maroon, navy, olive, orange, 
        purple, red, silver, teal, white, and yellow.

        By default, polar residues are green, negative residues red, positive blue,
        proline fuschia and all others are black.

        The function carrys out checking to ensure that
        1) Colors are valid (case insensitive)
        2) All 20 amino acids are accounted for
        
        No return value is provided, but an exception is raised if the operation
        cannot be completed.
        
        """
        self.SeqObj.set_HTMLColorResiduePalette(colorDict)



    def get_HTMLColorString(self):
        """
        Returns a string which contains a correctly formatted HTML string of the sequence
        defined with specific coloring.

        The sequence is divided into 10 residue blocks, with a newline (<br>) every
        50 residues.

        """
        return self.SeqObj.get_HTMLColorString()

        


    # ============================================ #
    # ============ IMPLICIT FUNCTIONS ============ #  
        
    #...................................................................................#
    def __unicode__(self):
        """ Returns the sequences """
        return "SequenceParameter [len="+str(len(self.SeqObj.seq)) + "], [seq='" + self.SeqObj.seq +"']" 

        
    #...................................................................................#
    def __str__(self):
        """ Returns the sequences """
        return self.__unicode__()

    def __len__(self):
        """ Returns the sequence length """
        return len(self.SeqObj.seq)


    
        
        
            



