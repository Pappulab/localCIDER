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

   Plotting is the backened for all localCIDER's plotting functionality.

   As with all functions in backend, nothing should be called directly, but in this case
   through the localcider.plots API module.

   Given plots are their own stateless function the plotting module is purely
   functional, no classes or state is maintained.

   For adding additional plotting features please see the paradigm defined with
   the linear plots - specifically using a hierarchical approach which passes
   functions. This makes for much cleaner code and reduces maintainance.


"""
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline


from sequence import Sequence
from backendtools import verifyType
from localciderExceptions import PlottingException

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
    Display a single-sequence Das-Pappu phase diagram plot on the screen
    """

    phaseplot_validate(fp, fn)

    # Create a plotting object
    initial_plottingObject = single_plot(fp, fn, label, fontSize)
    finalized_plottingObject = finalize_DasPappu(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # show the plot, or return the matplotlib fig object
    if getFig:
        return finalized_plottingObject
    else:
        finalized_plottingObject.show()


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
        saveFormat='pdf'):
    """
    Save a single-sequence Das-Pappu phase diagram to file
    """

    phaseplot_validate(fp, fn)

    # Delete any file in the filename location
    if(os.path.exists(filename)):
        os.remove(filename)

     # Create a plotting object
    initial_plottingObject = single_plot(fp, fn, label, fontSize)
    finalized_plottingObject = finalize_DasPappu(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # save the plot and then close the matloblib plotting object
    # note we default to 200 dpi for PNG files - this value
    # is hardcoded
    if saveFormat == 'png':
        finalized_plottingObject.savefig(filename, format=saveFormat, dpi=200)
    else:
        finalized_plottingObject.savefig(filename, format=saveFormat)

    finalized_plottingObject.close()


#...................................................................................#
def show_multiple_phasePlot(
        fp_list,
        fn_list,
        label=[],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Display multiple-sequences on a Das-Pappu phase diagram plot on the screen
    """

    # validate the various points
    for fp, fn in zip(fp_list, fn_list):
        phaseplot_validate(fp, fn)

    # Create a plotting object
    initial_plottingObject = multiple_plot(fp_list, fn_list, label, fontSize)
    finalized_plottingObject = finalize_DasPappu(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # show the plot, or return the matplotlib fig object
    if getFig:
        return finalized_plottingObject
    else:
        finalized_plottingObject.show()


#...................................................................................#
def save_multiple_phasePlot(
        fp_list,
        fn_list,
        filename,
        label=[],
        title="Diagram of states",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):
    """
    Save multiple-sequences on a Das-Pappu phase diagram to file
    """

    # validate the various points
    for fp, fn in zip(fp_list, fn_list):
        phaseplot_validate(fp, fn)

    # Create a plotting object
    initial_plottingObject = multiple_plot(fp_list, fn_list, label, fontSize)
    finalized_plottingObject = finalize_DasPappu(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # save the plot and then close the matloblib plotting object
    # save the plot and then close the matloblib plotting object
    # note we default to 200 dpi for PNG files - this value
    # is hardcoded
    if saveFormat == 'png':
        finalized_plottingObject.savefig(filename, format=saveFormat, dpi=200)
    else:
        finalized_plottingObject.savefig(filename, format=saveFormat)

    
    finalized_plottingObject.close()


# =========================
###
# UVERSKY PLOT FUNCTIONS
###
# =========================
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

    # Create a plotting object
    initial_plottingObject = single_plot(
        mean_net_charge, hydropathy, label, fontSize)
    finalized_plottingObject = finalize_uversky(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # show the plot, or return the matplotlib fig object
    if getFig:
        return finalized_plottingObject
    else:
        finalized_plottingObject.show()


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
    Function to save a single point on a Uversky plots

    """

    # Delete any file in the filename location
    if(os.path.exists(filename)):
        os.remove(filename)

     # Create a plotting object
    initial_plottingObject = single_plot(
        mean_net_charge, hydropathy, label, fontSize)
    finalized_plottingObject = finalize_uversky(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # save the plot and then close the matloblib plotting object
    if saveFormat == 'png':
        finalized_plottingObject.savefig(filename, format=saveFormat, dpi=200)
    else:
        finalized_plottingObject.savefig(filename, format=saveFormat)

    finalized_plottingObject.close()


#...................................................................................#
def show_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        label=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        getFig=False):
    """
    Function to show multiple points on a Uversky plots

    """

    # Create a plotting object
    initial_plottingObject = multiple_plot(
        mean_net_charge_list, hydropathy_list, label, fontSize)
    finalized_plottingObject = finalize_uversky(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # show the plot, or return the matplotlib fig object
    if getFig:
        return finalized_plottingObject
    else:
        finalized_plottingObject.show()


#...................................................................................#
def save_multiple_uverskyPlot(
        hydropathy_list,
        mean_net_charge_list,
        filename,
        label=[],
        title="Uversky plot",
        legendOn=True,
        xLim=1,
        yLim=1,
        fontSize=10,
        saveFormat='png'):

    # Create a plotting object
    initial_plottingObject = multiple_plot(
        mean_net_charge_list, hydropathy_list, label, fontSize)
    finalized_plottingObject = finalize_uversky(
        initial_plottingObject, legendOn, title, xLim, yLim)

    # save the plot and then close the matloblib plotting object
    if saveFormat == 'png':
        finalized_plottingObject.savefig(filename, format=saveFormat, dpi=200)
    else:
        finalized_plottingObject.savefig(filename, format=saveFormat)

    plt.close()


## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ##
##                                                              ##
##   INTERNAL FUNCTIONS FOR Uversky/Phase diagram plots         ##
##                                                              ##
## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ##

#...................................................................................#
def single_plot(x, y, label="", fontSize=10):
    """
    Internal function for creating a single sequence MATPLOTLIB object which can
    either then be saved or be displayed.

    INPUT:
    x          | x-coordinate to be plotted
    y          | y-coordinate to be plotted
    label      | label associated with that point (DEFAULT = "")
    fontSize  | size of label font (DEFAULT = 10)

    OUTPLOT:
    matplotlib object which can be plotted or saved

    """

    # set the x and y position of the mark
    x = float(x)
    y = float(y)

    # draw the actual plot
    plt.scatter(x, y, s=50, marker='o', color='Black', zorder=5)

    # if no label return plt object
    if label == "":
        return plt

    # we scale the location of the label in the event of extreme values
    if x > 0.8:
        x_lab = x - (0.01 * len(label) + 0.03)
    else:
        x_lab = x + 0.01

    if y > 0.9:
        y_lab = y - 0.025
    else:
        y_lab = y + 0.01

    # add the plot label
    plt.annotate(label, xy=(x_lab, y_lab), fontsize=fontSize)

    return plt


#...................................................................................#
def multiple_plot(x_list, y_list, label_list, fontSize):
    """
    Internal function for creating multiple induvidual points on a single set of axes.
    This plot can then be saved or displayed

    NOTE: This function ASSUMES x_list and y_list are lists of floats OF THE SAME
    LENGTH - this should be checked before!!

    INPUT:
    x_list     | list of x values to plot
    y_list     | list of y values to plot
    label_list | list of labels to associate with the points
    fontSize  | size of label font (DEFAULT = 10)

    OUTPLOT:
    matplotlib object which can be plotted or saved

    """

    # if we have no labels then construct the empty label list
    if len(label_list) == 0:
        # construct an empty label list
        label_list = []
        for i in xrange(0, len(x_list)):
            label_list.append("")

    # check that the three lists are the same length
    if not (len(x_list) == len(y_list) == len(label_list)):
        raise PlottingException(
            "Unequal length of positive fraction list, negative fraction list, and label list")

    # plot all the points
    for x, y, label in zip(x_list, y_list, label_list):
        plt.scatter(x, y, s=10, marker='o', color='Black', zorder=2)
        plt.annotate(label, xy=(x, y + 0.01), fontsize=fontSize)

    # annotate, set it all up, and return the plot object
    return plt


#...................................................................................#
def finalize_DasPappu(plt, legendOn, title, xLim, yLim):
    """
    Common function which finalizes up a plot by drawing on the regions 1-5, adding
    the legend and title. Used by both single and multiple phasePlot function

    """

    # define the five regions by filling the plot
    alphaval = 1
    reg1, = plt.fill([0, 0, 0.25], [0, 0.25, 0],
                     color='Chartreuse', alpha=alphaval, zorder=1)
    reg2, = plt.fill([0, 0, 0.35, 0.25], [0.25, 0.35, 0, 0],
                     color='MediumSeaGreen', alpha=alphaval, zorder=1)
    reg3, = plt.fill([0, 0.325, 0.675, 0.35], [0.35, 0.675,
                                               0.325, 0], color='DarkGreen', alpha=alphaval, zorder=1)

    reg4, = plt.fill([0, 0, 0.325], [0.35, 1, 0.675],
                     color='Red', alpha=alphaval, zorder=1)

    reg5, = plt.fill([0.35, 0.675, 1], [0, 0.325, 0],
                     color='Blue', alpha=alphaval, zorder=1)

    # set the plot limits
    plt.xlim([0, xLim])
    plt.ylim([0, yLim])

    # label the axes and set the title
    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')
    plt.xlabel(
        'Fraction of positively charged residues',
        fontproperties=axes_pro)
    plt.ylabel(
        'Fraction of negatively charged residues',
        fontproperties=axes_pro)

    # update the font property for the title
    axes_pro.set_size('x-large')
    plt.title(title, fontproperties=axes_pro)

    # if we the legend is on add the annotation
    if legendOn:

        # create and set set legend font options
        fontP = FontProperties()
        fontP.set_size('small')
        plt.legend([reg1, reg2, reg3, reg4, reg5],
                   ['Weak polyampholytes & polyelectrolytes:\nGlobules & tadpoles',
                    'Janus sequences:\nCollapsed or expanded - context dependent',
                    'Strong polyampholytes:\nCoils, hairpins, & chimeras',
                    'Negatively charged strong polyelectrolytes:\nSwollen coils',
                    'Positively charged strong polyelectrolytes:\nSwollen coils'],
                   prop=fontP)
    return plt


#...................................................................................#
def phaseplot_validate(fp, fn):
    """
    Validate if the fp and fn are resonable for making the diagram of states
    phase plots

    """

    # first check we can convert these guys into floating point numbers
    try:
        fp = float(fp)
    except ValueError as e:
        raise PlottingException(
            "Unable to convert " +
            str(fp) +
            " into a float")

    try:
        fn = float(fn)
    except ValueError as e:
        raise PlottingException(
            "Unable to convert " +
            str(fn) +
            " into a float")

    # next check they're both between 0 and 1
    if (fp < 0) or (fp > 1):
        raise PlottingException(
            "Fraction of positive residues outside of appropriate range [" + str(fp) + "]")
    if (fn < 0) or (fn > 1):
        raise PlottingException(
            "Fraction of positive residues outside of appropriate range [" + str(fn) + "]")

    # if we get here everything looks OK!


#...................................................................................#
def finalize_uversky(plt, legendOn, title, xLim, yLim):
    """
    Common function which finalizes up a plot by drawing on the regions 1-5, adding
    the legend and title. Used by both single and multiple phasePlot function

    """

    # define the five regions by filling the plot
    alphaval = 0.15

    # folded region
    reg1, = plt.fill([0.0, 0, 0.772], [1, 0.413, 1],
                     color='Chartreuse', alpha=0.25, zorder=1)

    # unfolded region
    reg2, = plt.fill([0, 0, 0.772, 1, 1],
                     [0, 0.413, 1, 1, 0],
                     color='Red',
                     alpha=0.15,
                     zorder=1)

    # set the plot limits
    plt.xlim([0, xLim])
    plt.ylim([0, yLim])

    # label the axes and set the title
    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')
    plt.xlabel('Mean net charge', fontproperties=axes_pro)
    plt.ylabel('Mean hydropathy <H>', fontproperties=axes_pro)

    # update the font property for the title
    axes_pro.set_size('x-large')
    plt.title(title, fontproperties=axes_pro)

    # if we the legend is on add the annotation
    if legendOn:

        # create and set set legend font options
        fontP = FontProperties()
        fontP.set_size('small')
        plt.legend([reg1, reg2], ['Folded proteins',
                                  'Natively unfolded'], prop=fontP)

    return plt


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#                           Linear sequence plots
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#...................................................................................#
def save_linearNCPR(SeqObj, blobLen, filename, saveFormat='png'):
    save_linearplot(build_NCPR_plot, SeqObj, blobLen, filename, saveFormat)


#...................................................................................#
def save_linearFCR(SeqObj, blobLen, filename, saveFormat='png'):
    save_linearplot(build_FCR_plot, SeqObj, blobLen, filename, saveFormat)


#...................................................................................#
def save_linearSigma(SeqObj, blobLen, filename, saveFormat='png'):
    save_linearplot(build_sigma_plot, SeqObj, blobLen, filename, saveFormat)


#...................................................................................#
def save_linearHydropathy(SeqObj, blobLen, filename, saveFormat='png'):
    save_linearplots(build_hydropathy_plot, SeqObj, blobLen, filename, saveFormat)

#...................................................................................#
def save_linearComplexity(complexityVector, complexityType, seqlen, filename, saveFormat='png'): 
    # Different structure to the others because here we have to pass a built complexity
    # vector into this function from the SeqObj
    

    plt = show_linearComplexity(complexityVector, complexityType, seqlen, getFig=True)
    
    if saveFormat == 'png':
        plt.savefig(filename, format=saveFormat, dpi=200)
    else:
        plt.savefig(filename, format=saveFormat)

    plt.close()

#...................................................................................#
def show_linearNCPR(SeqObj, blobLen, getFig=False):
    if getFig:
        return show_linearplot(build_NCPR_plot, SeqObj, blobLen, getFig)
    else:
        show_linearplot(build_NCPR_plot, SeqObj, blobLen, getFig)

#...................................................................................#
def show_linearFCR(SeqObj, blobLen, getFig=False):
    if getFig:
        return show_linearplot(build_FCR_plot, SeqObj, blobLen, getFig)
    else:
        show_linearplot(build_FCR_plot, SeqObj, blobLen, getFig)

#...................................................................................#
def show_linearSigma(SeqObj, blobLen, getFig=False):

    if getFig:
        return show_linearplot(build_sigma_plot, SeqObj, blobLen, getFig)
    else:
        show_linearplot(build_sigma_plot, SeqObj, blobLen, getFig)

#...................................................................................#
def show_linearHydropathy(SeqObj, blobLen, getFig=False):

    if getFig:
        return show_linearplot(build_hydropathy_plot, SeqObj, blobLen)
    else:
        show_linearplot(build_hydropathy_plot, SeqObj, blobLen)

#...................................................................................#
def show_linearComplexity(complexityVector, complexityType, seqlen, getFig=False):
    """
    The complexity plotting functions opperate outside of the general linear sequence
    framework as there are types of options/behaviours specific enough to the
    complexity plots that trying to shoe-horn them into the existing code would not
    be a good design decision.


    """
    
    # first generate the bar-plot and save the list of bars
    barlist = plt.bar(complexityVector[0,:], 
                      complexityVector[1,:],
                      width=1,
                      linewidth=1.0,
                      edgecolor='k',
                      color='#A8A8A8')

    # set the limits
    plt.ylim([0,1])
    plt.xlim([1, seqlen])

    # set the font properties
    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')

    # set the axis labels
    plt.xlabel('Residue', fontproperties=axes_pro)
    plt.ylabel('Complexity', fontproperties=axes_pro)
    
    # set the title (i.e. what type of complexity was calculated)
    axes_pro.set_size('x-large')
    if complexityType == 'WF':
        title='Wooton-Federhen complexity'
    elif complexityType == 'LC':
        title='Linguistic complexity'
    elif complexityType == 'LZW':
        title='Lempel-Ziv-Welch complexity'
    else:
        raise PlottingException('Unexpected complexity type passed - should never happen')

    plt.title(title, fontproperties=axes_pro)
    
    # finally either show the plot or return the plt object
    if getFig:
        return plt
    else:
        plt.show()


##
# Functions below allow construction of the various Linear Sequence Plots - i.e. should 
# be considered internal to this file
##

#...................................................................................#
def __build_linear_plot(
        data,
        title="",
        xlabel="Blob index",
        ylabel="",
        ylimits=[0,1],
        hline=None,
        setPositiveNegativeBars=False):
    """
    Internal function which expects data to be a Nx2 matrix (np.vstack) where column 1
    is the x values and column 2 is the y values. It also assumes the Y values
    are scaled to between 0 and 1 (so Y-axis limits are 0 and 1)

    hline defines a list of horizontal lines = so [0.2,-0.2] would draw horizontal lines


    """

    # plot the data
    barlist = plt.bar( data[0,:], 
                       data[1,:],
                       width=1,
                       linewidth=1.1,
                       edgecolor='k',
                       color='#A8A8A8')

    # this is really inefficient but means we have a consistent
    #
    if setPositiveNegativeBars:
        for bar in xrange(0, len(barlist)):
            if data[1, bar] < 0:
                barlist[bar].set_color('r')
                barlist[bar].set_edgecolor('k')
            else:
                barlist[bar].set_color('b')
                barlist[bar].set_edgecolor('k')

        # draw mao lines
        plt.plot([0, len(barlist)], [0.26, 0.26],
                 color='k', linewidth=1.5, linestyle="--")
        plt.plot([0, len(barlist)], [-0.26, -0.26],
                 color='k', linewidth=1.5, linestyle="--")

    # set Y lims
    plt.ylim(ylimits)
    plt.xlim([1, len(data[0, :])])

    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')

    # label
    plt.xlabel(xlabel, fontproperties=axes_pro)
    plt.ylabel(ylabel, fontproperties=axes_pro)

    axes_pro.set_size('x-large')
    plt.title(title, fontproperties=axes_pro)

    # return plot object
    return plt


#...................................................................................#
def build_NCPR_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the NCPR along your sequence divided into blobs of some size
    """

    try:
        plt = __build_linear_plot(SeqObj.linearDistOfNCPR(blobLen),
                                  title='NCPR distribution (blob ' + str(int(blobLen)) + ')',
                                  ylabel='NCPR', ylimits=[-1, 1], setPositiveNegativeBars=True)
    except PlottingException:
        raise PlottingException(
            "NCPR plot construction requires Sequence object")

    return plt

#...................................................................................#


def build_FCR_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the FCR along your sequence divided into blobs of some size
    """

    try:
        plt = __build_linear_plot(SeqObj.linearDistOfFCR(blobLen),
                                  title='FCR distribution (blob ' + str(int(blobLen)) + ')',
                                  ylabel='FCR')
    except PlottingException:
        raise PlottingException(
            "FCR plot construction requires Sequence object")

    return plt


#...................................................................................#
def build_sigma_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the sigma along your sequence divided into blobs of some size
    """

    try:
        plt = __build_linear_plot(SeqObj.linearDistOfSigma(blobLen),
                                  title='Sigma distribution (blob ' + str(int(blobLen)) + ')',
                                  ylabel='Sigma')

    except AttributeError as e:
        raise PlottingException(
            "Sigma plot construction requires Sequence object")

    return plt

#...................................................................................#


def build_hydropathy_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the Uversky-hydropathy along your sequence divided into blobs of some size
    """

    try:
        plt = __build_linear_plot(SeqObj.linearDistOfHydropathy(blobLen),
                                  title='Hydropathy distribution (blob ' + str(int(blobLen)) + ')',
                                  ylabel='Hydropathy')

    except AttributeError as e:
        raise PlottingException(
            "Hydropathy plot construction requires Sequence object")

    return plt


#...................................................................................#
def save_linearplot(build_fun, SeqObj, blobLen, filename, saveFormat='png'):
    """
    Internal function which builds and saves a linear sequence plot

    build_fun is the linear sequence building function we're using
    """

    plt = build_fun(SeqObj, blobLen)
    
    if saveFormat == 'png':
        plt.savefig(filename, format=saveFormat, dpi=200)
    else:
        plt.savefig(filename, format=saveFormat)

    plt.close()


#...................................................................................#
def show_linearplot(build_fun, SeqObj, blobLen, getFig=False):
    """
    Internal function which builds and shows a linear sequence plot

    build_fun is the linear sequence building function we're using
    """

    plt = build_fun(SeqObj, blobLen)
    if getFig:
        return plt
    else:
        plt.show()



#...................................................................................#
def save_local_composition_plot(residue_number, density_vectors, legend_color, legend_names, filename, saveFormat='png', max_val=-1, line_thickness=[], title='', plot_data=False):
    """
    Unlike the other plotting functions, the save_local_composition_plot is the only function for plotting local linear
    composition. It doesn't take advantage of the other builder functions, but instead exists as its own stand alone function.

    Arguments are as follows:

    residue_number - vector with 1-n numbering the residues

    density_vector - x by n np matrix where there are n values for each of the
                     n residues for x different residue groups

    legend_color - list of matplotlib compatible colors (must be the same length
                   as the number of groups)

    legend_names = list of names to assign for each group type (must be the same length
                   as the number of groups)

    filename - name of file to be saved

    save_format - format to be saved (default = 'png', 'pdf' also accepted)

    max_val - default maximum density value for residue group density (default =-1,
              which means the optimum value for dynamic range is calculated automatically
        

    """


    # determine the number of different groups
    n_groups = density_vectors.shape[0]

    # correct for bad starting values
    if max_val > 1 or max_val < 0:
        max_val = -1
        
    # determine the maximum density of any of the groups (note if a user-defined max_val is provided this
    # is skipped
    if max_val == -1:
        # this means max is 0.1 above max value OR 1.0, whicever is smaller
        max_val = min(np.max(density_vectors)+0.1, 1.0)

    # do some sanity checking to make sure all our ducks are in a row before we go all in
    if n_groups != len(legend_names):
        raise PlottingException('Mismatch in the number of of groups provided and the number of groups named for figure legend')

    # sanity check groups and legend color vector
    if  n_groups != len(legend_color):
        raise PlottingException('Mismatch in the number of of groups provided and the number of groups defined for line colors')
            

    ##
    ## The following section sets the figure size and font size, and is kind of an empyrical hack to make
    ## the generated figure look good and legible regardless of the number of residues. Works well!

    # if we're working with a sequence < 250 points
    if len(residue_number) < 250:
        width  = 20.0
        height = 12.5
        fs1    = 22
        fs2    = 28
        if len(line_thickness) == 0:
            for i in xrange(0, n_groups):
                line_thickness.append(1.5)

    # if we're working with a sequence > 250 points
    else:
        width  = 0.075*len(residue_number)
        height = 0.04*len(residue_number)
        fs1    = 0.08*len(residue_number)
        fs2    = 0.10*len(residue_number)
        if len(line_thickness) == 0:
            for i in xrange(0, n_groups):
                line_thickness.append(3.0)
        
    num_res = len(residue_number)
    
    handles_vector = []
    #fig = plt.figure(
    plt.figure(figsize=(width, height))

    font = {'family' : 'Bitstream Vera Sans',
            'weight' : 'normal',
            'size'   : fs1}

    matplotlib.rc('font', **font)

    
    for i in xrange(0, n_groups):
        
        # fit the density associated with group
        fx = UnivariateSpline(residue_number, density_vectors[i], s=1)

        if plot_data:
            h, = plt.plot(residue_number, density_vectors[i], color=legend_color[i], label=legend_names[i], linewidth=(line_thickness[i]*0.5))

        h, = plt.plot(residue_number, fx(residue_number), color=legend_color[i], label=legend_names[i], linewidth=line_thickness[i])
        handles_vector.append(h)


    # set the legend and bounding box positions so everything looks nice
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])                     
    plt.legend(handles=handles_vector, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=4)
               
    # set some asthetics stuff to make it look nice
    axes = plt.gca()
    axes.set_ylim([0,max_val])
    axes.set_xlim([1,num_res])
    plt.ylabel('Local Amino Acid Density', fontsize=fs2)

    # if we provided a title set that badboy
    if len(title) > 0:        
        plt.title(title, fontsize=fs2, loc='left')
    
    plt.savefig(filename, format=saveFormat)
    plt.close()




        

    
