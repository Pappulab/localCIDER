# localCIDER

`v 0.1.1 - August 14th 2014`


# Introduction

**CIDER** [Classification of Intrinsically Disordered Ensemble Regions](http://pappulab.wustl.edu/CIDER) is a [currently under construction] web-server for calculating various properties of a disordered sequence. CIDER's calculations are carried out by a backend app, called localCIDER.

**CIDER** is useful if you want to quickly calculate some parameter on the fly. HOWEVER, for doing a huge range of sequences, we recommend you instead use **localCIDER** and carry out the computation locally.

## localCIDER vs. CIDER?
**localCIDER** lets us provide an implicit distributed computing model to all of CIDER's calculations. Think of it this way - if 100 people have 100 sequences, we (the Pappu lab) have the burden of 10,000 calculations to perform. If everyone submitted these sequences at once, the website would slow for everyone, and maybe crash. This is a lot of work for us, and very annoying for you, the user. HOWEVER, by using **localCIDER** you can take advantage of

* Your own local infrastructure
* The ability to programmatically create analysis pipelines
* The fact you can operate independently of the CIDER web-server

Moreover, if there's a type of analysis you use frequently and you think other users would appreciate, we'd be happy to build this in to **localCIDER**. Adding new functionality is typically much more involved for a web-app.


# Installation

## Introduction
**localCIDER** version 0.1.0 was submitted to `pip` (the Python Package Index) on the 8th of August. This should be considered an alpha release, which is to say we're still ironing out creases, and it would be greatly appreciated if you could report any issues to alex.holehouse@wustl.edu. However, many of the features are ready to go!

## Installing on OSX or linux
We recommend installing using `pip`. `pip` is a command line interface for downloading and installing packages from the Python package index (PyPI). If you don't yet have pip installed [see the documentation here](http://pip.readthedocs.org/en/latest/installing.html).

We also recommend `virtualenv` although this is not at all required, it's just generally a good route to go if you use Python!

Once `pip` is installed **localCIDER** can be installed by running

    [sudo] pip install localcider

Very simple! **localCIDER** depends on both `numpy` and `matplotlib`, both of which will be installed via `pip` if you haven't already installed them.

If you encounter issues installing please let me know. That said, any issues *should* be due to problems with Numpy or `matplotlib` and not **localCIDER**, so investigate the documentation associated with installing those packages first.

### Known installation issues

#### OpenSuse
`pip` will not install system level dependencies, so you may need to install additional packages to get matplotlib working with vanilla `OpenSuse`. Specifically

    [sudo] zypper install freetype2-devel
    [sudo] zypper install gcc-c++
    [sudo] zypper install libpng16-devel
    [sudo] zypper install python-numpy-devel
    
Of course, one or more of these may well already be installed, so your best bet is to try installing **localcider** and if it fails check these guys are installed by using

    zypper search <package-name> # an 'i' next to the package indicates it's installed

An alternative issue may be that installing `python-numpy-devel` (which is required by pip to build `matplotlib`) requires numpy-1.7 and you have numpy-1.8. In *this* case we recommend using

    [sudo] zypper install python-matplotlib-tk
    
To install `matplotlib` (NOTE the `-tk` at the end!)



#### Ubuntu
The issues for OpenSuse are also true for Ubuntu - specifically `Ubuntu` seems to need the `python-dev` package to install `numpy`, and the `libfreetype6-dev` package for to get Freetype working in `matplotlib`

    [sudo] apt-get install python-dev
	[sudo] apt-get install libfreetype6-dev 


### Upgrading via pip
localCIDER does not automatically update, but for the next few months expect fairly regular updates. To see the historic update schedule see the end of this file.


`pip` allows you to upgrade packages to the latest version using

    [sudo] pip install --upgrade localcider
    
However, this command will upgrade all the dependencies as well, which may not be a desirable behaviour. To upgrade *just* **localCIDER** use

    [sudo] pip install -U --no-deps lolcacider 
    
    
## Installation on Windows
   
Tested on Windows 7 64 bit

1) Download and install Anaconda package (http://continuum.io/downloads)

2) Download `get-pip.py` (https://bootstrap.pypa.io/get-pip.py)

3) Open `cmd.exe` from start menu (this is windows command prompt)

4) Change directory to that containing `get-pip.py` (EX: `cd  "D:\downloads" ` )

5) In `cmd.exe` run `python get-pip.py`

6) In `cmd.exe` run `pip install localcider`

7) Open IPython Notebook from the windows start menu (installed with Anaconda package, other editors can be used i.e. PyCharm)

8) Write code to use localCIDER package
 
 
# How to use
Once you have localcider installed via pip, fire up Python and `import` the localcider package;

    import localcider
    
You'll be greeted with a message detailing the current version.

There are two main localcider modules right now;

* `sequenceParameters` which lets you calculate many different parameters, as well as show or save plots relating to a specific sequence - ***this is the main module of interest***
* `plots` which lets you plot certain parameters independent of sequenceParameters

In the next release (V 0.2.0) we'll also introduce the `sequencePermutant` code, which employs a variety of approaches to generate permutants.

#### Quick-start - calculate kappa of my sequence and plot it on the Pappu-Das diagram of states 

The following code will print the kappa value of your sequence to the terminal prompt, and then plot your sequence on a Pappu-Das phase plot, saving it in the current directory as `my_snazzy_plot.png`

    # import the relevant code
    from localcider.sequenceParameters import SequenceParameters
    
    # create a SequenceParameters object with your amino acid sequence
    SeqOb = SequenceParameters("MY AMINO ACID SEQUENCE")
    
    # carry out analysis
    SeqOb.get_kappa()
    SeqOb.save_phaseDiagramPlot('my_snazzy_plot')
    
    
#### Interactive help

**localCIDER** contains highly complete interactive help - for any object or function simply type `help(<item>)` for a complete description. For example;

    # import the relevant code
    from localcider.sequenceParameters import SequenceParameters
    
    # create a SequenceParameters object with your amino acid sequence
    SeqOb = SequenceParameters("MY AMINO ACID SEQUENCE")
    
    help(SeqOb)
    
    help(SeqOb.get_kappa())


# General sequence analysis with SequenceParameters

The approach we recommend for accessing SequenceParameters objects is to use the following Python code;


    from localcider.sequenceParameters import SequenceParameters
    
    
By opening your code with this line you now have direct access to the `SequenceParameters` class, which takes either a string of an amino acid sequence or the filename of a file containing an amino acid sequence, which is then read and parsed. As an example;


    # as before we import the SequenceParameters class directly
    from localcider.sequenceParameters import SequenceParameters
	
	# the sequence below is the first 30 residues from alpha-synuclein
	seqOb = SequenceParameters("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAA")     
    
Or alternatively
    
    # as before we import the SequenceParameters class directly
    from localcider.sequenceParameters import SequenceParameters
	
	# the sequence below is the first 30 residues from alpha-synuclein
	seqOb = SequenceParameters(sequenceFile="syn.fasta")
	
	

Both these code snippets create a `SequenceParameters` object (here called `seqOb`, but obviously this variable could be named anything), on which we can carry out a huge range of operations. The complete function list is shown below for reference.

Most of these functions don't take arguments. Optional arguments are prefixed with a question mark (?). For each function we use the `seqOb.<function>` syntax - e.g.

    seqOb.get_FCR()


### SequenceParameters functions

Function name | Operation 
:---: | :---: 
`get_FCR()`  | Get the fraction of charged residues
`get_NCPR()` | Get the net charge per residue
`get_countNeg()` | Get the number of negatively charged amino acids 
`get_countPos()` | Get the number of positivly charged amino acids
`get_fraction_negative()` | Get the fraction of residues which are negatively charged
`get_fraction_positive()` | Get the fraction of residues which are positivly charged
`get_countNeut()` | Get the number of neutral amino acids
`get_amino_acid_fractions()` | Get a dictionary of the fractions of each amino acid in the sequence
`get_fraction_disorder_promoting()` | Get the fraction of residues predicted to be 'disorder promoting' [REF 1]
`get_kappa()` | Get the kappa value of the sequence [REF 2]
`get_mean_net_charge()` | Get the absolute mean net charge of your sequence
`get_phase_plot_region()` | Get the region on the Pappu-Das diagram of states where your sequence falls [2]
`get_mean_hydropathy()` | Get the mean hydrophobicity as calculated from a skewed Kyte-Doolittle hydrophobicity scale* [REF 3]
`get_uversky_hydropathy()` | Get the mean hydropathy as calculated from a normalized Kyte-Doolittle hydrophobicity scale** [REF 3 REF 4]

\* The skewed hydrophobicity scale shifts the normal KD scale such that the lowest value is 0 (instead of -4.5) and the highest value is 9 (instead of 4.5)

\** The normalized Kyte-Doolittle scale converts all values on the scale to fall between 0 and 1

### Phosphorylation functions
The following functions augment your sequence to consider the impact of phosphorylation on the electrostatic properties. Note this makes the highly simplifying assumption that the phosphorylation of a Ser/Thr/Tyr residue simply adds a negative charge to your protein chain. In reality, many other properties of the chain are impacted by phosphorylation than simply the linear charge patterning.

Function name | Operation 
:---: | :---: 
`get_all_phosphorylatable_sites()` | Get a list of the positions of all Ser/Thr/Tyr residues. Note this function DOES NOT make any inference about accessibility or binding/recognition motifs, it's literally just a list of sites which canonically are able to be phosphorylated.
`set_phosphosites(list_of_positions)` | Sets a residues which may be phosphorylatable (must be T/Y/S). Note this can be called multiple times, and will only set sites which are T/Y/S)
`get_phosphosites()` | Get the list of sites currently designated as phosphosites by set_phosphosites
`clear_phosphosites()` | Clears any previously defined phosphosites
`get_kappa_after_phosphorylation()` | Get the kappa value assuming full phosphorylation (with H<sub>2</sub>PO<sub>4</sub><sup>-1</sup>).
`get_phosphosequence()` | Returns the fully phosphorylated sequence, with "phosphorylated residues" replaced with glutamate (E)
`get_full_phosphostatus_kappa_distribution()`|  ** UNDER DEVELOPMENT ** For each phospho-permutant, given the available phosphosites, calculate the kappa, fraction positive, fraction negative, FCR, NCP and meanHydropathy. This returns a list of tuples. We're still determining if all this information is necessary. Also note that as the number of phosphosites increases the number of calculations here scales as n^2. Be warned!

### Plotting functions
The following functions let you plot parameters from your sequence in a set of defined plots.

Function name | Operation 
:---: | :---: 
`show_phaseDiagramPlot(?label)`| Renders a `matplotlib` Pappu-Das diagram of states plot with your sequence on the diagram<sup>[2]</sup>. `label` can be a string which labels your sequence on the plot.
`show_uverskyDiagramPlot(?label)`| Renders a `matplotlib` Uversky plot with your sequence on the diagram<sup>[4]</sup>. `label` can be a string which labels your sequence on the plot.
`show_linearHydropathy(blobLen)`| Renders a `matplotlib` plot of the moving average hydropathy along the sequence, where the hydropathy is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used (**This made be changed to the default value in future releases**)
`show_linearNCPR(blobLen)` | Renders a `matplotlib` plot of the moving average net charge per residue (NCPR) along the sequence, where the NCPR is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used (**This made be changed to the default value in future releases**)
`show_linearSigma(blobLen)` | Renders a `matplotlib` plot of the moving average sigma parameter along the sequence, where sigma is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used (**This made be changed to the default value in future releases**). Recall that sigma is calculated as the NCPR<sup>2</sup> / FCR. 
`save_phaseDiagram_plot(filename, ?label)` | Renders a `matplotlib` Pappu-Das diagram of states plot with your sequence on the diagram<sup>[2]</sup>. `label` can be a string which labels your sequence on the plot. The plot is saved in the `filename` location (`.png` is appended).
`save_uverskyDiagram_plot(filename, ?label)` | Renders a `matplotlib` Uversky plot with your sequence on the diagram<sup>[4]</sup>. `label` can be a string which labels your sequence on the plot. The plot is saved in the `filename` location (`.png` is appended).
`save_linearHydropathy(blobLen, filename)` | Renders a `matplotlib` plot of the moving average hydropathy along the sequence, where the hydropathy is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used (**This made be changed to the default value in future releases**). The plot is saved in the `filename` location (`.png` is appended).
`save_linearNCPR(blobLen, filename)`| Renders a `matplotlib` plot of the moving average net charge per residue (NCPR) along the sequence, where the NCPR is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used (**This made be changed to the default value in future releases**). The plot is saved in the `filename` location (`.png` is appended).
`save_linearSigma(blobLen, filename)` | Saves a `matplotlib` generated plot of the moving avesequenceParametersrage sigma parameter along the sequence, where sigma is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used (**This made be changed to the default value in future releases**). Recall that sigma is calculated as the NCPR<sup>2</sup> / FCR. The plot is saved in the `filename` location (`.png` is appended).


# Plotting with the plots module

The `plots` module lets you generate or save plots without explicitly invoking sequenceParametersobjects. Such plots can either be generated using the raw values, or by passing in sequenceParametersobjects.

Unlike the `SequenceParameters`, `plots` is a collection of functions, and does not allow for the creation of `Plots` objects.

In the following discussion a **single** plot contains information from just one sequence, whereas a **multiple** plot contains sequence information from many sequences

### Using `Plots` functions

There are two modes with which you can use plots functions;

#### Passing raw parameters to be plotted (single and multiple plots)
This may be convenient if you're plotting previously calculated values, or plotting a theoretical distribution here. Either individual values or lists of values can be passed.

### Passing sequenceParametersobjects (multiple plots only)
If you're analyzing a set of sequenceParametersobjects you may wish to place them on the same plot. Note we do not offer this functionality for single sequences, because you can always plot a single sequence using the associated `SequenceParameters` function (e.g. `.show_uverskyPlot()` or `.save_phaseDiagramPlot()`).




Function name | Input Arguments | Discussion
:---: | :---: | :---:
`show_single_phasePlot(fp, fn, label='', title='Diagram of states', legendOn=True)` | `fp` is the fraction of positive residues in a sequence, `fn` is the fraction of negative residues in a sequence. `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included.  | No programatic output is generated, but the plot is rendered on the screen. 
`show_single_uverskyPlot(hydropathy, mean_net_charge, label='', title='Uversky plot', legendOn=True)` | `hydropathy` is the sequence's 0-1 normalized hydropathy, while `mean_net_charge` is the the absolute value of the sequence's mean net charge. `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included.  | No programatic output is generated, but the plot is rendered on the screen.
`save_single_phasePlot(fp, fn, filename, label='', title='Diagram of states', legendOn=True)` | `fp` is the fraction of positive residues in a sequence, `fn` is the fraction of negative residues in a sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included.  | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_single_uverskyPlot(hydropathy, mean_net_charge, filename, label='', title='Uversky plot', legendOn=True)` | `hydropathy` is the sequence's 0-1 normalized hydropathy, while `mean_net_charge` is the the absolute value of the sequence's mean net charge. `filename` is the name of the file to save the plot as (`.png` is appended). `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included.  | No programatic output is generated, but the plot is constructed and saved at the `filename` location.

If you wish to plot a single sequence's Pappu-Das diagram of states plot or Uversky plot from a `SequenceParameters` object we suggest using the appropriate `SequenceParameters` functions (e.g.  `.show_uverskyDiagram()`, `.show_phaseDiagram()`)

### `*_multiple*_` functions

For the `multiple` plotting functions, we pass in a list of values or objects corresponding to the multiple sequences of interest. 

There are two types of multiple plots

**1)** Those which you pass two lists of the actual values of interest, where the lists represent parallel lists, where the position in each list corresponds to a specific sequence.
**2)** Those which you pass a list of `SequenceParameters` objects to, meaning you only pass it a single list


The `show_*` functions render a plot on the screen in front of you through `matplotlib`. These are useful for directly interacting with sequences and examining how specific sequence might behave.

The `save_*` functions save a plot to the designated file name.

Function name | Input Arguments | Discussion
:---: | :---: | :---:
`show_multiple_phasePlot(fp_list, fn_list, ?label_list=[], ?title="Diagram of states", ?legendOn=True)`  | `fp_list` and `fn_list` are lists of the positive fraction and negative fraction of sequences, where the list position in each of these two 'parallel' lists corresponds to the same sequence. That is to say, if you were looking at 5 sequences, `fp_list` and `fn_list` would both be of length 5, and `fp_list[2]` and `fn_list[2]` would represent the fraction position and fraction negative of your third sequence, respectively. `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. | No programatic output is generated, but the plot is rendered to the screen.  
`show_multiple_phasePlot2(SeqParam_list, label_list=[], title='Diagram of states', legendOn=True)` | All arguments are the same as above, except now instead of explicit fraction of positive and fraction of negative lists, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is rendered to the screen.
`show_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, label_list=[], title='Uversky plot', legendOn=True)` | `hydropathy_list` and `mean_net_charge_list` are lists of the 0-1 normalized hydropathy and the absolute mean net charge, where the list position in each of these two 'parallel' lists corresponds to the same sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. | No programatic output is generated, but the plot is rendered to the screen. 
`show_multiple_uverskyPlot2(SeqParam_list, label_list=[], title='Uversky plot', legendOn=True)` | All arguments are the same as above, except now instead of having an explicit hydropathy list and an explicit mean net charge list, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is rendered to the screen. 
`save_multiple_phasePlot(fp_list, fn_list, filename, ?label_list=[], ?title="Diagram of states", ?legendOn=True)`  | `fp_list` and `fn_list` are lists of the positive fraction and negative fraction of sequences, where the list position in each of these two 'parallel' lists corresponds to the same sequence. That is to say, if you were looking at 5 sequences, `fp_list` and `fn_list` would both be of length 5, and `fp_list[2]` and `fn_list[2]` would represent the fraction position and fraction negative of your third sequence, respectively. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_multiple_phasePlot2(SeqParam_list, filename, label_list=[], title='Diagram of states', legendOn=True)` | All arguments are the same as above, except now instead of explicit fraction of positive and fraction of negative lists, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is constructed and saved at the `filename` location.
`save_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, filename, label_list=[], title='Uversky plot', legendOn=True)` | `hydropathy_list` and `mean_net_charge_list` are lists of the 0-1 normalized hydropathy and the absolute mean net charge, where the list position in each of these two 'parallel' lists corresponds to the same sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_multiple_uverskyPlot2(SeqParam_list, filename, label_list=[], title='Uversky plot', legendOn=True)` | All arguments are the same as above, except now instead of having an explicit hydropathy list and an explicit mean net charge list, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 


# `localcider` exceptions

If you're into exception handling, localcider has a set of defined exceptions (all of which inhert from the `Exception` base class) which can be access by importing localcider - i.e.

    import localcider
    
    localcider.SequenceFileParserException
    localcider.SequenceException
    localcider.ResTableException
    localcider.PlottingException
    
There are a couple more but these relate to functionality not yet available in version 0.1.0

# Tutorials 


## 1) Calculate a bunch of parameters associated with a sequence I found on the internet

    from localcider.sequenceParameters import SequenceParameters
    
    # define my internet sequence and create a new sequenceParametersobject with itxs
    internetSequence="VAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"
    object=SequenceParameters(internetSequencexs)
    
    
    # Let's get some parameters!xs
    object.get_kappa()                       # 0.1960
    
    object.get_fraction_negative()           # 0.1818
    
    object.get_fraction_positive()           # 0.0909
    
    object.get_fraction_disorder_promoting() # 0.75
    
    object.get_amino_acid_fractions()        # {'A': 0.136,
                                             #  'C': 0.0,
                                             #  'D': 0.0682,
                                             #  'E': 0.1136,
                                             #  'F': 0.0227,
                                             #  'G': 0.1364,
                                             #  'H': 0.0,
                                             #  'I': 0.0455,
                                             #  'K': 0.0909,
                                             #  'L': 0.0455,
                                             #  'M': 0.0227,
                                             #  'N': 0.0227,
                                             #  'P': 0.0682,
                                             #  'Q': 0.0682,
                                             #  'R': 0.0,
                                             #  'S': 0.0227,
                                             #  'T': 0.0455,
                                             #  'V': 0.0909,
                                             #  'W': 0.0,
                                             #  'Y': 0.0}


##### 2) Get a list of kappa values for a directory of `fasta` files.

`fasta` files are defined as having the following structure


    > Header, preceded by a '>' character
    SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
    SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
    SEQUENCESEQUENCESEQUENCESEQUENCESEQUENCESEQUENCE
    
Let's imagine we have a directory with 10 fasta files (`seq0.fasta`, `seq1.fasta`, ..., `seq9.fasta`). The following script will calculate the kappa of all such sequences

    # import a few things...
    import os
    from localcider.sequenceParameters import SequenceParameters
    
    # filelist is now a list of each file in the current directory
    filelist = os.listdir(".")
    
    # create an empty list 
    list_of_SeqObjs = []
    
    # populate that list with sequenceParametersobjects, which we construct from the
    # sequence found in each file
    for file in filelist:
       list_of_SeqObjs.append(SequenceParameters(sequenceFile=file))

	# for each 
	for obj in list_of_SeqObjs:
	   print obj.get_kappa()
	   
Uh oh! One of our sequence files has an amino acid "X" in it - this isn't a cannonical amino acid, so triggers a `SequenceFileParserException`. We actually just want to skip any such files, so we change the code above slightly. Specifically, let's change the for-loop that loops over each file;

    # import a few things...
    import os
    from localcider.sequenceParameters import SequenceParameters
    
    ## NOTE - we now import the whole localcider module to give us access
    ## to the exception classes
    import localcider
    
    # get a list of all files in the current directory
    filelist = os.listdir(".")
    
    # create an empty list 
    list_of_SeqObjs = []
    
    # populate that list with sequenceParameters objects, which we construct from the
    # sequence found in each file
    for file in filelist:
       list_of_SeqObjs.append(SequenceParameters(sequenceFile=file))
       
       
    # populate that list with SequenceParameters objects, which we construct from the
    # sequence found in each file
    for file in filelist:
       try:
	      list_of_SeqObjs.append(SequenceParameters(sequenceFile=file))
	   except localcider.SequenceFileParserException:
	      # if we encounter a file parsing error just skip that sequence
	      continue


## Example 3 - plotting two sequence on a single Pappu-Das phase diagram

	# import the plots module
    from localcider import plots
    
    # import the SequenceParameters class from the sequenceParameters class
    from localcider.sequenceParameters import SequenceParameters    
    
    # create two SequenceParameters objects
    SynWT=SequenceParameters("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAA")
    SynMut=SequenceParameters("MDVFMKGLSKAKEGEEKKAEKTKQGVAEAA")
    
    
    # showing with sequence labels
    plots.show_multiple_phasePlot2([SynWT, SynMut],label_list=['WT','MUT'])
    
    # saving without
    plots.save_multiple_uverskyPlot2([SynWT, SynMut], 'example', label_list=['WT','MUT'])


# Epilog

## localCIDER roadmap

We have a fairly well defined roadmap for the next few sets of minor release for localCIDER, with the first major release coming sometime thereafter

#### V 0.2.0
* Turn on the sequence permutant generator and Wang-Landau density of kappa states analysis 

#### V 0.3.0
* Provide stand alone, command line tools for the more computationally intensive operations (to make deploying on HPC infrastructure 


## About
localCIDER was developed in the Pappu lab by Alex Holehouse, with much of the original code written by James Ahad. Importantly, we believe scientific software should be held to same (or higher) software standards as general commercial and open source projects. To this end we use Git for version control, employ standard conventions for versioning, run an ever-expanding set of unit tests, and try to ensure all code is well structured and well documented. The source code is fully open and [can be viewed here](https://github.com/Pappulab/localCIDER).


## References

**[1]** **A. Campen, R.M. Williams, C.J. Brown, J. Meng, V.N. Uversky, A.K. Dunker**. (2008) [*Protein Pept Lett.* 15(9) pp 956 - 963](http://www.ncbi.nlm.nih.gov/pubmed/18991772)

**[2]** **R.K. Das & R.V. Pappu** (2013) [*PNAS* 110, 33, pp 13392 - 13397](http://www.pnas.org/content/110/33/13392.short)

**[3]** **J. Kyte, R. Doolittle** (1982) [*J. Mol. Biol.* 157 pp 105 - 132](http://www.ncbi.nlm.nih.gov/pubmed/7108955)

**[4]** **V. Uversky** (2002) [*Protein Sci. Apr* 11(4) pp 739 - 756](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2373528/)


## Thanks
Many people have been involved in this project. We'll try and include an up-to-date list here

* Paul Nobrega (Matthews lab) for bug reports and the Windows installation instructions

## Update schedule

* version 0.1.0 [August 10th 2014]
* version 0.1.1 [August 14th 2014] (minor bug fixes, more tests, improved sequence parsing)

----------
        <ul>
          <li><a href="https://github.com/Pappulab/localCIDER">View On <strong>GitHub</strong></a></li>
          <li><a href="http://pappulab.wustl.edu">Go to <strong>Pappu lab</strong></a></li>
          <li><a href="http://pappulab.wustl.edu/CIDER/">Go to <strong>CIDER</strong> </a></li>
        </ul>



