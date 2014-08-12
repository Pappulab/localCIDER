""" 
   !--------------------------------------------------------------------------!
   ! LICENSE INFO:                                                            !
   !--------------------------------------------------------------------------!
   !    This file is part of localCIDER.                                      !
   !                                                                          !
   !    Version 0.1.0                                                         !
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
   ! MAIN AUTHOR:   James Ahad and Alex Holehouse                             !
   !                                                                          !
   !--------------------------------------------------------------------------!

   
   File Description:
   ================
   
   Plotting is the backened for all localCIDER's plotting functionality. 

   As with all functions in backend, nothing should be called directly, but
   instead through the localcider.plots API module.

   Given plots are their own stateless function the plotting module is purely
   functional, no classes or state is maintained.

   For adding additional plotting features please see the paradigm defined with
   the linear plots - specifically using a hierarchical approach which passes
   functions. This makes for much cleaner code and reduces maintainance.
  

"""
import os
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from sequence import Sequence
from backendtools import verifyType
from localciderExceptions import PlottingException

#...................................................................................#
def show_single_phasePlot(fp, fn, label="",title="Diagram of states",legendOn=True):
    """
    Display a single-sequence Das-Pappu phase diagram plot on the screen
    """
    
    phaseplot_validate(fp,fn)
    
    # Create a plotting object
    initial_plottingObject = single_plot(fp, fn, label)
    finalized_plottingObject = finalize_DasPappu(initial_plottingObject, legendOn, title)

    # show the plot
    plt.show()


#...................................................................................#
def save_single_phasePlot(fp, fn, filename, label="", title="Diagram of states",legendOn=True):
    """
    Save a single-sequence Das-Pappu phase diagram to file
    """
    
    phaseplot_validate(fp,fn)
     
    # Delete any file in the filename location
    if(os.path.exists(filename)):
        os.remove(filename)
        
     # Create a plotting object
    initial_plottingObject = single_plot(fp, fn, label)
    finalized_plottingObject = finalize_DasPappu(initial_plottingObject, legendOn, title)


    # save the plot and then close the matloblib plotting object
    plt.savefig(filename,dpi=200)
    plt.close()


#...................................................................................#
def show_multiple_phasePlot(fp_list, fn_list, label=[], title="Diagram of states",legendOn=True):
    """
    Display multiple-sequences on a Das-Pappu phase diagram plot on the screen
    """

    # validate the various points
    for fp,fn in zip(fp_list,fn_list):
        phaseplot_validate(fp,fn)

    # Create a plotting object
    initial_plottingObject = multiple_plot(fp_list, fn_list, label)
    finalized_plottingObject = finalize_DasPappu(initial_plottingObject, legendOn, title)

    # show the plot
    finalized_plottingObject.show()


#...................................................................................#
def save_multiple_phasePlot(fp_list, fn_list, filename, label=[], title="Diagram of states",legendOn=True):
    """
    Save multiple-sequences on a Das-Pappu phase diagram to file
    """
    
    # validate the various points
    for fp,fn in zip(fp_list,fn_list):
        phaseplot_validate(fp,fn)
    
    # Create a plotting object
    initial_plottingObject = multiple_plot(fp_list, fn_list, label)
    finalized_plottingObject = finalize_DasPappu(initial_plottingObject, legendOn, title)

    # save the plot and then close the matloblib plotting object
    plt.savefig(filename,dpi=200)
    plt.close()



### =========================
###
###   UVERSKY PLOT FUNCTIONS
###
### =========================
#...................................................................................#
def show_single_uverskyPlot(hydropathy, mean_net_charge, label="",title="Uversky plot",legendOn=True):
             
    #uverskyPlot_validate(hyd,fn)

    # Create a plotting object
    initial_plottingObject = single_plot(mean_net_charge, hydropathy, label)
    finalized_plottingObject = finalize_uversky(initial_plottingObject, legendOn, title)

    # show the plot
    plt.show()


#...................................................................................#
def save_single_uverskyPlot(hydropathy, mean_net_charge, filename, label="", title="Uversky plot",legendOn=True):

    #uverskyPlot_validate(fp,fn)
     
    # Delete any file in the filename location
    if(os.path.exists(filename)):
        os.remove(filename)
        
     # Create a plotting object
    initial_plottingObject = single_plot(mean_net_charge, hydropathy, label)
    finalized_plottingObject = finalize_uversky(initial_plottingObject, legendOn, title)


    # save the plot and then close the matloblib plotting object
    plt.savefig(filename,dpi=200)
    plt.close()


#...................................................................................#
def show_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, label=[], title="Uversky plot",legendOn=True):
    
    # validate the various points
    #for fp,fn in zip(fp_list,fn_list):
    #    phaseplot_validate(fp,fn)

    # Create a plotting object
    initial_plottingObject = multiple_plot(mean_net_charge_list, hydropathy_list, label)
    finalized_plottingObject = finalize_uversky(initial_plottingObject, legendOn, title)

    # show the plot
    plt.show()


#...................................................................................#
def save_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, filename, label=[], title="Uversky plot",legendOn=True):
    
    # validate the various points
    #for fp,fn in zip(fp_list,fn_list):
    #    phaseplot_validate(fp,fn)
    
    # Create a plotting object
    initial_plottingObject = multiple_plot(mean_net_charge_list, hydropathy_list, label)
    finalized_plottingObject = finalize_uversky(initial_plottingObject, legendOn, title)

    # save the plot and then close the matloblib plotting object
    plt.savefig(filename,dpi=200)
    plt.close()


        

## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ##
##                                                              ##
##   INTERNAL FUNCTIONS FOR Uversky/Phase diagram plots         ##
##                                                              ##
## <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> ##

#...................................................................................#
def single_plot(x,y,label=""):
    """
    Internal function for creating a single sequence MATPLOTLIB object which can
    either then be saved or be displayed.

    INPUT:
    fp         | Fraction of positive residues
    fn         | Fraction of negative residues
    label      | Lable to associated with sequence [DEFAULT = ""]
    legendOn   | Write legend in plot [DEFAULT = True]
    title      | Title for plot [DEFAULT = "Diagrams of states"]

    OUTPLOT:
    matplotlib object which can be plotted or saved
        
    """
    
    # see if the input are SequenceParameter objects

    
    
    # set the x and y position of the mark
    x=float(x)
    y=float(y)
    
    # draw the actual plot, with the sequence annotated
    plt.scatter(x,y,s=50,marker='o',color='Black',zorder=5)
    
    # we scale the location of the label in the event of extreme values
    if x > 0.8:
        x_lab=x-(0.01*len(label)+0.03)
    else:
        x_lab=x+0.01

    if y > 0.9:
        y_lab=y-0.025
    else:
        y_lab=y+0.01

    plt.annotate(label,xy=(x_lab,y_lab))

    return plt


#...................................................................................#
def multiple_plot(x_list,y_list,label_list=[],legendOn=True,title="Diagram of states"):
    """
    Internal function for creating a single sequence MATPLOTLIB object which can
    either then be saved or be displayed.

    NOTE: This function ASSUMES fp_list and fn_list are lists of floats OF THE SAME
    LENGTH - this should be checked before!!

    INPUT:
    fp         | List of fraction of positive residues
    fn         | List of fraction of negative residues

    
    label      | List of lables to associated with sequence [DEFAULT = ""]
    legendOn   | Write legend in plot [DEFAULT = True]
    title      | Title for plot [DEFAULT = "Diagrams of states"]

    OUTPLOT:
    matplotlib object which can be plotted or saved
        
    """

    # if we have no labels then construct the empty label list
    if len(label_list) == 0:
        # construct an empty label list
        label_list = []
        for i in xrange(0,len(x_list)):
            label_list.append("")

    # check that the three lists are the same length
    if not (len(x_list) == len(y_list) == len(label_list)):
        raise PlottingException("Unequal length of positive fraction list, negative fraction list, and table list")

                    
    # plot all the points
    for x,y,label in zip(x_list,y_list,label_list):
        plt.scatter(x,y,s=20,marker='o',color='Black',zorder=2)
        plt.annotate(label,xy=(x+.01,y+.01))

    # annotate, set it all up, and return the plot object
    return plt
    

#...................................................................................#
def finalize_DasPappu(plt, legendOn, title):
    """
    Common function which finalizes up a plot by drawing on the regions 1-5, adding
    the legend and title. Used by both single and multiple phasePlot function

    """
    
    # define the five regions by filling the plot
    alphaval=1
    reg1, = plt.fill([0,0,.25],[0,.25,0],color = 'Chartreuse',alpha=alphaval, zorder=1)
    reg2, = plt.fill([0,0,.35,.25],[.25,.35,0,0],color = 'MediumSeaGreen',alpha=alphaval,zorder=1)
    reg3, = plt.fill([0,.35,.65,.35],[.35,.65,.35,0],color = 'DarkGreen',alpha=alphaval,zorder=1)
    reg4, = plt.fill([0,0,.35],[.35,1,.65],color = 'Red',alpha=alphaval,zorder=1)
    reg5, = plt.fill([.35,.65,1],[0,.35,0],color = 'Blue',alpha=alphaval,zorder=1)
        
    # set the plot limits
    plt.ylim([0,1])
    plt.xlim([0,1])
    
    # label the axes and set the title
    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')
    plt.xlabel('Fraction of positively charged residues', fontproperties = axes_pro)
    plt.ylabel('Fraction of negatively charged residues', fontproperties = axes_pro)

    # update the font property for the title
    axes_pro.set_size('x-large')
    plt.title(title,fontproperties = axes_pro)
    
    # if we the legend is on add the annotation
    if legendOn:
        
        # create and set set legend font options
        fontP = FontProperties()
        fontP.set_size('small')
        plt.legend([reg1,reg2,reg3,reg4,reg5],
                   ['Weak polyampholytes & polyelectrolytes:\nGlobules & tadpoles',
                    'Janus sequences:\nCollapsed or expanded - context dependent',
                    'Strong polyampholytes:\nCoils, hairpins, & chimeras',
                    'Negatively charged strong polyelectrolytes:\nSwollen coils',
                    'Positively charged strong polyelectrolytes:\nSwollen coils'],
                   prop = fontP)
    return plt
    

#...................................................................................#
def phaseplot_validate(fp,fn):
    """
    Validate if the fp and fn are resonable for making the diagram of states
    phase plots

    """

    # first check we can convert these guys into floating point numbers
    try:
        fp=float(fp)
    except ValueError, e:
        raise PlottingException("Unable to convert " + str(fp) + " into a float")

    try:
        fn=float(fn)
    except ValueError, e:
        raise PlottingException("Unable to convert " + str(fn) + " into a float")
                            
    # next check they're both between 0 and 0.5
    if (fp < 0) or (fp > 1):
        raise PlottingException("Fraction of positive residues outside of appropriate range [" + str(fp) + "]")
    if (fn < 0) or (fn > 1):
        raise PlottingException("Fraction of positive residues outside of appropriate range [" + str(fn) + "]")

    # if we get here everything looks OK!


#...................................................................................#
def finalize_uversky(plt, legendOn, title):
    """
    Common function which finalizes up a plot by drawing on the regions 1-5, adding
    the legend and title. Used by both single and multiple phasePlot function

    """
    
    # define the five regions by filling the plot
    alphaval=0.15
    
    # folded region
    reg1, = plt.fill([0.0,0,0.772],[1,0.413,1],color = 'Chartreuse',alpha=0.25, zorder=1)

    # unfolded region
    reg2, = plt.fill([0,  0,  0.772,  1,  1], 
                     [0,  0.413,  1,  1,  0],
                     color = 'Red',
                     alpha=0.15,
                     zorder=1)
        
    # set the plot limits
    plt.ylim([0,1])
    plt.xlim([0,1])
    
    # label the axes and set the title
    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')
    plt.xlabel('Mean net charge', fontproperties = axes_pro)
    plt.ylabel('Mean hydropathy <H>', fontproperties = axes_pro)

    # update the font property for the title
    axes_pro.set_size('x-large')
    plt.title(title,fontproperties = axes_pro)
    
    # if we the legend is on add the annotation
    if legendOn:
        
        # create and set set legend font options
        fontP = FontProperties()
        fontP.set_size('small')
        plt.legend([reg1,reg2],['Folded proteins','Natively unfolded'], prop=fontP)


    return plt


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#                           Linear sequence plots 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#...................................................................................#
def save_linearNCPR(SeqObj, blobLen, filename):
    save_linearplot(build_NCPR_plot, SeqObj, blobLen, filename)


#...................................................................................#
def save_linearSigma(SeqObj, blobLen, filename):
    save_linearplot(build_sigma_plot, SeqObj, blobLen, filename)


#...................................................................................#
def save_linearHydropathy(SeqObj, blobLen, filename):
    save_linearplots(build_hydropathy_plot, SeqObj, blobLen, filename)


#...................................................................................#
def show_linearNCPR(SeqObj, blobLen):
    show_linearplot(build_NCPR_plot, SeqObj, blobLen)


#...................................................................................#
def show_linearSigma(SeqObj, blobLen):
    show_linearplot(build_sigma_plot, SeqObj, blobLen)


#...................................................................................#
def show_linearHydropathy(SeqObj, blobLen):
    show_linearplot(build_hydropathy_plot, SeqObj, blobLen)

##
## Functions below allow construction of the various Linear Sequence Plots #DRY
##


def __build_linear_plot(data, title="", xlabel="Blob index", ylabel="", ylimits=[-0.1, 1.1]):
    """
    Internal plot which expects data to be a Nx2 matrix (np.vstack) where column 1
    is the x values and column 2 is the y values. It also assumes the Y values
    are scaled to between 0 and 1 (so Y-axis limits are -0.1 and +1.1)
    """
    
    # plot the data
    plt.plot(data[0,:], data[1,:])

    # set Y lims
    plt.ylim( ylimits)

    axes_pro = FontProperties()
    axes_pro.set_size('large')
    axes_pro.set_weight('bold')
  

    # label
    plt.xlabel(xlabel,fontproperties = axes_pro)
    plt.ylabel(ylabel,fontproperties = axes_pro)

    axes_pro.set_size('x-large')
    plt.title(title, fontproperties = axes_pro)

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
                                  title='NCPR distribution (blob ' + str(int(blobLen))+')',
                                  ylabel='NCPR',ylimits=[-1.1,1.1])                              
    except PlottingException:
        raise PlottingException("NCPR plot construction requires Sequence object")
    
    return plt


#...................................................................................#
def build_sigma_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the sigma along your sequence divided into blobs of some size
    """

    try:
        plt = __build_linear_plot(SeqObj.linearDistOfSigma(blobLen), 
                                  title='Sigma distribution (blob ' + str(int(blobLen))+')',
                                  ylabel='Sigma' )                              

    except AttributeError as e:
        raise PlottingException("Sigma plot construction requires Sequence object")

    return plt

#...................................................................................#
def build_hydropathy_plot(SeqObj, blobLen):
    """
    Function which returns a matplotlib.plt object ready for saving/plotting
    of the Uversky-hydropathy along your sequence divided into blobs of some size
    """


    try:
        plt = __build_linear_plot(SeqObj.linearDistOfHydropathy(blobLen),
                                  title='Hydropathy distribution (blob ' + str(int(blobLen))+')',
                                  ylabel='Hydropathy')

    except AttributeError as e:
        raise PlottingException("Hydropathy plot construction requires Sequence object")

    return plt


#...................................................................................#
def save_linearplot(build_fun, SeqObj, blobLen, filename):
    """
    Internal function which builds and saves a linear sequence plot 

    build_fun is the linear sequence building function we're using
    """

    plt = build_fun(SeqObj, blobLen)
    plt.savefig(filename,dpi=200)
    plt.close()


#...................................................................................#
def show_linearplot(build_fun, SeqObj, blobLen):
    """
    Internal function which builds and shows a linear sequence plot 

    build_fun is the linear sequence building function we're using
    """

    plt = build_fun(SeqObj, blobLen)
    plt.show()

