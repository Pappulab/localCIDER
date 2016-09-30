"""
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.9                                                         !
   !--------------------------------------------------------------------------!

   File Description:
   ================

   This is a file containing various plotting functions. Note that localCIDER
   doesn't require these functions to work - the ability to plot things is
   nice but not essential.

   There are eight plotting functions which can be called, and are
   as follows

   >>> Show or save a single sequence on the diagam of states plot
   - show_single_phasePlot
   - save_single_phasePlot

   >>> Show or save multiple sequences on the diagram of states plot
   - show_multiple_phasePlot
   - save_multiple_phasePlot
   - show_multiple_phasePlot2
   - save_multiple_phasePlot2

   >>> Show or save a single sequence on the Uversky plot
   - show_single_uverskyPlot
   - save_single_uverskyPlot

   >>> Show or save multiple sequences on the Uversky plot
   - show_multiple_uverskyPlot
   - save_multiple_uverskyPlot
   - show_multiple_uverskyPlot2
   - save_multiple_uverskyPlot2


   These functions can be called independently of other localCIDER
   functionality, or you could calculate the fraction of positive and
   negative residues for a bunch of sequences and use these functions
   independent of the localCIDER SequenceParameter objects.

   Note that for plotting linear sequence data (e.g. linear hydropathy,
   linear NCPR etc) you must use a SequenceParameter object.

"""

from backend import plotting


#...................................................................................#
def show_single_phasePlot(
        fp,
        fn,
        label="",
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plot a single sequence on the Pappu-Das phase plot (diagram of states).

    INPUT:
    --------------------------------------------------------------------------------
    fp         | Fraction of positive residues
    fn         | Fraction of negative residues

    label      | On-plot label of sequence (DEFAULT = no label)
    title      | Plot title (DEFAULT = "Diagram of states")
    legendOn   | Include the phase diagram region legend (DEFAULT = True)
    xLim       | Set upper limit for the x axis (DEFAULT = 1)
    yLim       | Set upper limit for the y axis (DEFAULT = 1)
    fontSize   | Set font size for the point labels (DEFAULT = 10)
    getFig     | Returns a matplotlib figure object instead of simply displaying the
               | plot on the screen (DEFAULT = False)

    OUTPUT:
    --------------------------------------------------------------------------------

    Nothing but a plot should be generated on screen
    """

    return plotting.show_single_phasePlot(
        fp, fn, label, title, legendOn, xLim, yLim, fontSize, getFig)


#...................................................................................#
def save_single_phasePlot(
        fp,
        fn,
        filename,
        label="",
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plot a single sequence on the Pappu-Das phase plot (diagram of states).

    INPUT:
    --------------------------------------------------------------------------------
    fp         | Fraction of positive residues
    fn         | Fraction of negative residues
    filename   | Path/name of file to save plot (.png is appended)

    label      | On-plot label of sequence (DEFAULT = no label)
    title      | Plot title (DEFAULT = "Diagram of states")
    legendOn   | Include the phase diagram region legend (DEFAULT = True)
    xLim       | Set upper limit for the x axis (DEFAULT = 1)
    yLim       | Set upper limit for the x axis (DEFAULT = 1)
    fontSize   | Set font size for the point labels (DEFAULT = 10)
    saveFormat | Defines the file formal to save plots as. This parameter
                 is passed to matplotlibs savefig command which supports 
                 the following filetypes: emf, eps, pdf, png, ps, raw, 
                 rgba, svg, svgz. (DEFAULT = png)

    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing but a plot should be generated on screen

    """

    plotting.save_single_phasePlot(
        fp,
        fn,
        filename,
        label,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        saveFormat)


#...................................................................................#
def show_multiple_phasePlot(
        fp_list,
        fn_list,
        label=[""],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plot multiple sequences on the same Pappu-Das phase plot (diagram of states).

    INPUT:
    --------------------------------------------------------------------------------
    fp_list    | Fraction of positive residues
    fn_list    | Fraction of negative residues

    label_list | On-plot label of sequence (DEFAULT = no label)
    title      | Plot title (DEFAULT = "Diagram of states")
    legendOn   | Include the phase diagram region legend (DEFAULT = True)
    xLim       | Set upper limit for the x axis (DEFAULT = 1)
    yLim       | Set upper limit for the x axis (DEFAULT = 1)
    fontSize   | Set font size for the point labels (DEFAULT = 10)
    getFig     | Returns a matplotlib figure object instead of simply displaying the
               | plot on the screen (DEFAULT = False)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a plot with multiple points should appear on the screen

    """

    return plotting.show_multiple_phasePlot(
        fp_list,
        fn_list,
        label,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        getFig)


#...................................................................................#
def show_multiple_phasePlot2(
        SeqParam_list,
        label_list=[],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plot a single sequence on the Pappu-Das phase plot (diagram of states). This function takes
    SequenceParameter objects rather than raw values

    INPUT:
    --------------------------------------------------------------------------------
    SeqParam_list | list of sequence parameter objects

    label         | On-plot label of sequence (DEFAULT = no label)
    title         | Plot title (DEFAULT = "Diagram of states")
    legendOn      | Include the phase diagram region legend (DEFAULT = True)
    xLim          | Set upper limit for the x axis (DEFAULT = 1)
    yLim          | Set upper limit for the y axis (DEFAULT = 1)
    fontSize      | Set font size for the point labels (DEFAULT = 10)
    getFig        | Returns a matplotlib figure object instead of simply displaying the
                  | plot on the screen (DEFAULT = False)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a plot with multiple points should appear on the screen

    """

    # construct the fraction positive and negative vectors
    fp_list = []
    fn_list = []
    for seq in SeqParam_list:
        fp_list.append(seq.get_fraction_positive())
        fn_list.append(seq.get_fraction_negative())

    return plotting.show_multiple_phasePlot(
        fp_list,
        fn_list,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        getFig)


#...................................................................................#
def save_multiple_phasePlot(
        fp_list,
        fn_list,
        filename,
        label_list=[],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plot multiple sequences on the same Pappu-Das phase plot (diagram of states) and save that file

    INPUT:
    --------------------------------------------------------------------------------
    fp_list    | Fraction of positive residues
    fn_list    | Fraction of negative residues
    filename   | name of file to save (.png is appended)

    label_list | On-plot label of sequence (DEFAULT = no label)
    title      | Plot title (DEFAULT = "Diagram of states")
    legendOn   | Include the phase diagram region legend (DEFAULT = True)
    xLim       | Set upper limit for the x axis (DEFAULT = 1)
    yLim       | Set upper limit for the y axis (DEFAULT = 1)
    fontSize   | Set font size for the point labels (DEFAULT = 10)
    saveFormat | Defines the file formal to save plots as. This parameter
                 is passed to matplotlibs savefig command which supports 
                 the following filetypes: emf, eps, pdf, png, ps, raw, 
                 rgba, svg, svgz. (DEFAULT = png)



    OUTPUT:
    --------------------------------------------------------------------------------
    No output, but if succesful a file with filename is generated with the associated plot

    """

    plotting.save_multiple_phasePlot(
        fp_list,
        fn_list,
        filename,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        saveFormat)


#...................................................................................#
def save_multiple_phasePlot2(
        SeqParam_list,
        filename,
        label_list=[],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plot multiple sequences on the same Pappu-Das phase plot (diagram of states) and save that file.
    This function takes SequenceParameter objects rather than raw values.

    INPUT:
    --------------------------------------------------------------------------------
    SeqParam_list | list of sequence parameter objects
    filename      | name of file to save (.png is appended)

    label_list    | On-plot label of sequence (DEFAULT = no label)
    title         | Plot title (DEFAULT = "Diagram of states")
    legendOn      | Include the phase diagram region legend (DEFAULT = True)
    xLim          | Set upper limit for the x axis (DEFAULT = 1)
    yLim          | Set upper limit for the y axis (DEFAULT = 1)
    fontSize      | Set font size for the point labels (DEFAULT = 10)
    saveFormat    | Defines the file formal to save plots as. This parameter
                    is passed to matplotlibs savefig command which supports 
                    the following filetypes: emf, eps, pdf, png, ps, raw, 
                    rgba, svg, svgz. (DEFAULT = png)

    OUTPUT:
    --------------------------------------------------------------------------------
    No output, but if succesful a file with filename is generated with the associated plot

    """

    # construct the fraction positive and negative vectors
    fp_list = []
    fn_list = []
    for seq in SeqParam_list:
        fp_list.append(seq.get_fraction_positive())
        fn_list.append(seq.get_fraction_negative())

    plotting.save_multiple_phasePlot(
        fp_list,
        fn_list,
        filename,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        saveFormat)


####################################################
##
# UVERSKY PLOTS BELOW
##
####################################################

#...................................................................................#
def show_single_uverskyPlot(
        hydropathy,
        mean_net_charge,
        label="",
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plots a single sequence on the Uversky plot (hydropathy vs. mean net charge)

    INPUT:
    --------------------------------------------------------------------------------
    hydropathy      | Mean hydropathy for sequence
    mean_net_charge | Absolute magnitude of the protein's net charge divided by sequence length

    label           | On-plot label of sequence (DEFAULT = no label)
    title           | Plot title (DEFAULT = "Uversky plot")
    legendOn        | Include the phase diagram region legend (DEFAULT = True)
    xLim            | Set upper limit for the x axis (DEFAULT = 1)
    yLim            | Set upper limit for the y axis (DEFAULT = 1)
    fontSize        | Set font size for the point labels (DEFAULT = 10)
    getFig          | Returns a matplotlib figure object instead of simply displaying the
                    | plot on the screen (DEFAULT = False)



    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but an annotated Uversky plot should be generated on screen

    """

    return plotting.show_single_uverskyPlot(
        hydropathy,
        mean_net_charge,
        label,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        getFig)


#...................................................................................#
def save_single_uverskyPlot(
        hydropathy,
        mean_net_charge,
        filename,
        label="",
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plots a single sequence on the Uversky plot (hydropathy vs. mean net charge) and save it to 'filename' (.png is
    appended).

    INPUT:
    --------------------------------------------------------------------------------
    hydropathy      | Mean hydropathy for sequence
    mean_net_charge | Absolute magnitude of the protein's net charge divided by sequence length
    filename        | Path/name of file to save plot (.png is appended)

    label           | On-plot label of sequence (DEFAULT = no label)
    title           | Plot title (DEFAULT = "Uversky plot")
    legendOn        | Include the phase diagram region legend (DEFAULT = True)
    xLim            | Set upper limit for the x axis (DEFAULT = 1)
    yLim            | Set upper limit for the y axis (DEFAULT = 1)
    fontSize        | Set font size for the point labels (DEFAULT = 10)
    saveFormat      | Defines the file formal to save plots as. This parameter
                      is passed to matplotlibs savefig command which supports 
                      the following filetypes: emf, eps, pdf, png, ps, raw, 
                      rgba, svg, svgz. (DEFAULT = png)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing but single uversky plot should be saved to disk

    """

    plotting.save_single_uverskyPlot(
        hydropathy,
        mean_net_charge,
        filename,
        label,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        saveFormat)


#...................................................................................#
def show_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        label_list=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plots multiple sequences on the Uversky plot (hydropathy vs. mean net charge) and shows the plot on the screen.

    INPUT:
    --------------------------------------------------------------------------------
    hydropathy_list      | List of proteins' mean hydropathy
    mean_net_charge_list | List of the absolute magnitude of the protein's net charge divided by sequence length

    label_list           | List of labels for each sequence (empty = no list)
    title                | Plot title (DEFAULT = "Uversky plot")
    legendOn             | Include the phase diagram region legend (DEFAULT = True)
    xLim                 | Set upper limit for the x axis (DEFAULT = 1)
    yLim                 | Set upper limit for the y axis (DEFAULT = 1)
    fontSize             | Set font size for the point labels (DEFAULT = 10)
    getFig               | Returns a matplotlib figure object instead of simply displaying the
                         | plot on the screen (DEFAULT = False)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a Uversky plot with multiple points should appear on the screen

    """
    return plotting.show_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        getFig)


#...................................................................................#
def show_multiple_uverskyPlot2(
        SeqParam_list,
        label_list=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Plots multiple sequences on the Uversky plot (hydropathy vs. mean net charge) and shows the plot on the screen.
    This function takes a list of SequenceParameter objects instead of the actual values

    INPUT:
    --------------------------------------------------------------------------------
    SeqParam_list        | List of sequence parameter objects

    label_list           | List of labels for each sequence (empty = no list)
    title                | Plot title (DEFAULT = "Uversky plot")
    legendOn             | Include the phase diagram region legend (DEFAULT = True)
    xLim                 | Set upper limit for the x axis (DEFAULT = 1)
    yLim                 | Set upper limit for the y axis (DEFAULT = 1)
    fontSize             | Set font size for the point labels (DEFAULT = 10)
    getFig               | Returns a matplotlib figure object instead of simply displaying the
                         | plot on the screen (DEFAULT = False)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a Uversky plot with multiple points should appear on the screen

    """

    # construct the fraction positive and negative vectors
    hydropathy_list = []
    mean_net_charge_list = []
    for seq in SeqParam_list:

        hydropathy_list.append(seq.get_uversky_hydropathy())
        mean_net_charge_list.append(seq.get_mean_net_charge())

    return plotting.show_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        getFig)


#...................................................................................#
def save_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        filename,
        label_list=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plots multiple sequences on the Uversky plot (hydropathy vs. mean net charge) and saves that plot to 'filename' (.png
    is appended).

    INPUT:
    --------------------------------------------------------------------------------
    hydropathy           | List of mean hydropathies for sequences
    mean_net_charge_list | List of the absolute magnitude of the protein's net charge divided by sequence length
    filename             | Path/name of file to save plot (.png is appended)

    label_list           | List of labels for each sequence (empty = no list)
    title                | Plot title (DEFAULT = "Uversky plot")
    legendOn             | Include the phase diagram region legend (DEFAULT = True)
    xLim                 | Set upper limit for the x axis (DEFAULT = 1)
    yLim                 | Set upper limit for the y axis (DEFAULT = 1)
    fontSize             | Set font size for the point labels (DEFAULT = 10)
    saveFormat           | Defines the file formal to save plots as. This parameter
                           is passed to matplotlibs savefig command which supports 
                           the following filetypes: emf, eps, pdf, png, ps, raw, 
                           rgba, svg, svgz. (DEFAULT = png)



    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a Uversky plot with multiple points should be saved to the 'filename', where filename
    defines the name of a file/path for saving

    """

    plotting.save_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        filename,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize,
        saveFormat)


#...................................................................................#
def save_multiple_uverskyPlot2(
        SeqParam_list,
        filename,
        label_list=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Plots multiple sequences on the Uversky plot (hydropathy vs. mean net charge) and saves that plot to 'filename' (.png
    is appended). This function takes SequenceParameter objects instead of list of hydropathy and mean_net_charge
    values

    INPUT:
    --------------------------------------------------------------------------------
    SeqParam_list   | List of sequence parameter objects
    filename        | Path/name of file to save plot (.png is appended)

    label           | On-plot label of sequence (DEFAULT = no label)
    title           | Plot title (DEFAULT = "Uversky plot")
    legendOn        | Include the phase diagram region legend (DEFAULT = True)
    xLim            | Set upper limit for the x axis (DEFAULT = 1)
    yLim            | Set upper limit for the y axis (DEFAULT = 1)
    fontSize        | Set font size for the point labels (DEFAULT = 10)
    saveFormat      | Defines the file formal to save plots as. This parameter
                      is passed to matplotlibs savefig command which supports 
                      the following filetypes: emf, eps, pdf, png, ps, raw, 
                      rgba, svg, svgz. (DEFAULT = png)


    OUTPUT:
    --------------------------------------------------------------------------------
    Nothing, but a Uversky plot with multiple points should be saved to the 'filename', where filename
    defines the name of a file/path for saving

    """

    # construct the fraction positive and negative vectors
    hydropathy_list = []
    mean_net_charge_list = []
    for seq in SeqParam_list:
        hydropathy_list.append(seq.get_uversky_hydropathy())
        mean_net_charge_list.append(seq.get_mean_net_charge())

    plotting.save_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        filename,
        label_list,
        title,
        legendOn,
        xLim,
        yLim,
        fontSize)
