# localCIDER

`v 0.1.13 - *** 2017`


# Introduction

**CIDER** [Classification of Intrinsically Disordered Ensemble Regions](http://pappulab.wustl.edu/CIDER) is a web-server for calculating various properties of a disordered sequence. CIDER's calculations are carried out by a backend app, called localCIDER.

**CIDER** is useful if you want to quickly calculate some parameter on the fly. HOWEVER, for doing a large number of sequences, we recommend you instead use **localCIDER** and carry out the computation locally.

## localCIDER vs. CIDER?
**localCIDER** lets us provide an implicit distributed computing model to all of CIDER's calculations. Think of it this way - if 100 people have 100 sequences, we (the Pappu lab) have the burden of 10,000 calculations to perform. If everyone submitted these sequences at once, the website would slow for everyone, and maybe crash. This is a lot of work for us, and very annoying for you, the user. HOWEVER, by using **localCIDER** you can take advantage of

* Your own local infrastructure
* The ability to programmatically create analysis pipelines
* The fact you can operate independently of the CIDER web-server

Moreover, if there's a type of analysis you use frequently and you think other users would appreciate, we'd be happy to build this in to **localCIDER**. Adding new functionality is typically much more involved for a web-app.

# Installation

## Introduction
**localCIDER** version 0.1.10 was submitted to PyPI (the Python Package Index) in November 2016. This should be considered a stable release, however, if you encounter any issues/bugs it would be greatly appreciated if you could report any issues to alex.holehouse@wustl.edu. 

## New Features
1) Added molecular weight 


## Installing on OSX or linux
We recommend installing using `pip`. `pip` is a command line interface for downloading and installing packages from the Python package index (PyPI). If you don't yet have pip installed [see the documentation here](http://pip.readthedocs.org/en/latest/installing.html).

We also recommend `virtualenv` although this is not at all required, it's just generally a good route to go if you use Python! For more information on `virtualenv` check out [this link](https://virtualenv.pypa.io/en/latest/)

Once `pip` is installed, **localCIDER** can be installed by running;

    [sudo] pip install localcider 

Very simple! **localCIDER** depends on the three pillars of scientific Python computing (`numpy`, `scipy`, and  `matplotlib`), all of which will be installed via `pip` if you haven't already installed them.

If you encounter issues installing please let me know. That said, any issues *should* be due to problems with Numpy or `matplotlib` and not **localCIDER**, so investigate the documentation associated with installing those packages first.

## Python 2 and Python 3
As of version 0.1.8 localCIDER is supported under both Python 2.x and 3.x. localCIDER has been extensively tested and used with Python 2.7.x, and we recommend this version of Python. If you encounter issues using localCIDER with Python version 3.x please let us knnow! 

### Known installation issues

#### OpenSuse
`pip` will not install system level dependencies, so you may need to install additional packages to get matplotlib working with vanilla `OpenSuse`. Specifically

    [sudo] zypper install freetype2-devel
    [sudo] zypper install gcc-c++
    [sudo] zypper install libpng16-devel
    [sudo] zypper install python-numpy-devel
    
Of course, one or more of these may already be installed, so your best bet is to try installing **localcider**, and if it fails check these guys are installed by using

    zypper search <package-name> # an 'i' next to the package indicates it's installed

An alternative issue may be that installing `python-numpy-devel` (which is required by pip to build `matplotlib`) requires numpy-1.7 and you have numpy-1.8. In *this* case we recommend using

    [sudo] zypper install python-matplotlib-tk
    
To install `matplotlib` (NOTE the `-tk` at the end!)

Getting Python packaging and OpenSuse to work together can be a pain, so there are likely alternative solutions to those presented here.

#### Ubuntu
The issues for OpenSuse are also true for Ubuntu - specifically `Ubuntu` seems to need the `python-dev` package to install `numpy`, and the `libfreetype6-dev` package for to get Freetype working in `matplotlib`

    [sudo] apt-get install python-dev
    [sudo] apt-get install libfreetype6-dev 


### Upgrading via pip
localCIDER does not automatically update, but can be updated using pip. To see the historic update schedule see the end of this file.

`pip` allows you to upgrade packages to the latest version using

    [sudo] pip install --upgrade localcider
    
However, this command will upgrade all the dependencies as well, which may not be a desirable behaviour. To upgrade *just* **localCIDER** use

    [sudo] pip install -U --no-deps localcider

Alternativly, you can update by uninstalling and re-installing `localcider` using

    pip uninstall localcider
    pip install localcider 
    
    
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
 
 
# Using localCIDER
Once you have localCIDER installed, activate the interactive Python prompt (we recommend using [iPython](https://ipython.org/) and `import` the localcider package;

    import localcider
    
There are two main localcider modules right now;

* `sequenceParameters` which lets you calculate many different parameters, as well as show or save plots relating to a specific sequence - ***this is the main module of interest***
* `plots` which lets you plot certain parameters independent of sequenceParameters


#### Quick-start - calculate kappa of my sequence and plot it on the Das-Pappu diagram of states.

The following code will print the kappa value of your sequence to the terminal prompt, and then plot your sequence on a Das-Pappu diagram of states, saving it in the current directory as `my_snazzy_plot.png`

    # import the relevant code
    from localcider.sequenceParameters import SequenceParameters
    
    # create a SequenceParameters object with your amino acid sequence
    SeqOb = SequenceParameters("DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV")
    
    # carry out analysis
    SeqOb.get_kappa() 
    SeqOb.save_phaseDiagramPlot('my_snazzy_plot')
    
    # For publications we recommend using vector graphics and saving plots as PDF files 
    SeqOb.save_phaseDiagramPlot('my_snazzy_plot', saveFormat='pdf'))
    
    
#### Interactive help

**localCIDER** contains an extensive set of interactive help information - for any object or function simply type `help(<item>)` for a complete description. For example;

    # import the relevant code
    from localcider.sequenceParameters import SequenceParameters
    
    # create a SequenceParameters object with your amino acid sequence
    SeqOb = SequenceParameters("MY AMINO ACID SEQUENCE")
    
    help(SeqOb)
    
    help(SeqOb.get_kappa())


# General sequence analysis with SequenceParameters

The approach we recommend for accessing SequenceParameters objects is to use the following Python code;


    from localcider.sequenceParameters import SequenceParameters
    
    
By opening your code with this line, you now have direct access to the `SequenceParameters` class, which takes either a string of an amino acid sequence or the filename of a file containing an amino acid sequence, which is then read and parsed. As an example;


    # as before we import the SequenceParameters class directly
    from localcider.sequenceParameters import SequenceParameters
	
	# the sequence below is the first 30 residues from alpha-synuclein
	SeqOb = SequenceParameters("MDVFMKGLSKAKEGVVAAAEKTKQGVAEAA")     
    
Or alternatively
    
    # as before we import the SequenceParameters class directly
    from localcider.sequenceParameters import SequenceParameters
	
	# the sequence below is the first 30 residues from alpha-synuclein
	SeqOb = SequenceParameters(sequenceFile="syn.fasta")
	
	

Both these code snippets create a `SequenceParameters` object - here that object is called `SeqOb`, but obviously this variable could be named anything. We can run a huge range of analysis routines on this object. The complete function list is shown below for reference.

Many of these functions don't take arguments. Optional arguments are prefixed with a question mark (?). For each function we use the `seqOb.<function>` syntax - e.g.

    seqOb.get_FCR()


### Single value sequence analysis functions
The functions below perform various analysis over sequences and return a single value. NOTE: Where pH values can be provided, if left blank we assume a neutral pH where only R/K/D/E are charged. If a pH value is provided, then R/K/D/E/C/Y/H are all considered titratable residues using EMBOSS pKa values, listed below:
'C': 8.5, 'Y': 10.1, 'H': 6.5, 'E': 4.1, 'D': 3.9, 'K': 10.0, 'R': 12.5


Function name | Operation 
:---: | :---: 
`get_length()`  | Get the sequence length
`get_FCR(pH=None)`  | Get the fraction of charged residues in the sequence [4] (pH keyword allows for a pH specific value)
`get_NCPR(pH=None)` | Get the net charge per residue of the sequence [5]
`get_isoelectric_point()` | Get the isoelectric point of the sequence
`get_molecular_weight()` | Get the molecular weight of the protein associated with a given amino acid sequence 
`get_countNeg()` | Get the number of negatively charged residues in the sequence (D/E)
`get_countPos()` | Get the number of positively charged residues in the sequence (R/K)
`get_countNeut()` | Get the number of neutral amino acids
`get_fraction_negative()` | Get the fraction of residues which are negatively charged (F-)
`get_fraction_positive()` | Get the fraction of residues which are positively charged (F+)
`get_fraction_expanding(pH=None)` | Get the fraction of residues which are predicted to contribute to chain expansion (E/D/R/K/P)
`get_amino_acid_fractions()` | Get a dictionary of the fractions of each amino acid in the sequence
`get_fraction_disorder_promoting()` | Get the fraction of residues predicted to be 'disorder promoting' [1]. Note this is NOT a disorder prediction!
`get_kappa()` | Get the sequence's kappa value [2]
`get_Omega()` | Get the sequence's Omega value. Omega defines the patterning between charged/proline residues and all other residues 14].
`get_mean_net_charge(pH=None)` | Get the absolute mean net charge of your sequence
`get_phase_plot_region()` | Get the region on the Das-Pappu diagram of states where your sequence falls [2]
`get_mean_hydropathy()` | Get the mean hydropathy as calculated from a skewed Kyte-Doolittle hydrophobicity scale* [3]
`get_uversky_hydropathy()` | Get the mean hydropathy as calculated from a normalized Kyte-Doolittle hydrophobicity scale\*\* [3,4]
`get_PPII_propensity(mode='hilser')` | Get the overall sequence's PPII propensity as defined by one of three PPII propensity scales. By default, the scale by Elam et al.[6] is used, but modes `creamer` and `kallenbach` are also available, which use values from scales by Rucker or Shi, respectively [12,13]. 
`get_delta()` <br> | Returns the delta value of the sequence, as defined when calculating kapp [2] 
`get_deltaMax()` <br> | Returns the maximum possible delta value (delta-max) for a sequence of this composition




\* The skewed hydrophobicity scale shifts the normal KD scale such that the lowest value is 0 (instead of -4.5) and the highest value is 9 (instead of 4.5)

\** The normalized Kyte-Doolittle scale converts all values on the scale to fall between 0 and 1


### Position-specific sequence analysis functions
The following functions generate an array of values which describes some property associated with the sequence as a function of sequence position.

Function name | Operation 
:---: | :---: 
`get_linear_FCR(blobLen=5)`  | Returns a 2D numpy array of the fraction of charged residues (FCR) as defined by a sliding window. The first dimension of the 2D array contains the local FCR values and the second contains the associated residue index values along the sequence. the `blobLen` keyword defines the window size used to calculate the sequence-local FCR. A stepsize of 1 is always used.
`get_linear_NCPR(blobLen=5)`  | Returns a 2D numpy array of the net charge per residue (NCPR) as defined by a sliding window. The first dimension of the 2D array contains the local NCPR values and the second contains the associated residue index values along the sequence. The `blobLen` keyword defines the window size used to calculate the sequence-local NCPR. A stepsize of 1 is always used.
`get_linear_sigma(blobLen=5)`  | Returns a 2D numpy array of the linear sigma value as defined by a sliding window. The first dimension of the 2D array contains the local sigma values and the second contains the associated residue index values along the sequence. the `blobLen` keyword defines the window size used to calculate the sequence-local sigma value. A stepsize of 1 is always used.
`get_linear_hydropathy(blobLen=5)`  | Returns a 2D numpy array of the local hydropathy as defined by a sliding window using the Kyte-Doolite hydropathy scale[3]. The first dimension of the 2D array contains the local hydropathy score and the second contains the associated residue index values along the sequence. The `blobLen` keyword defines the window size used to calculate the sequence-local hydropathy. A stepsize of 1 is always used.
`get_linear_complexity(complexityType='WF', alphabetSize=20, userAlphabet={}, windowSize=10, stepSize=1, wordSize=3)` | Returns the linear sequence complexity as defined by complexityType. Optionally, the sequence complexity of a reduced complexity alphabet can be returned, where that reduced alphabet is defined by either the alphabetSize (which takes advantage of 11 pre-defined simplified alphabets) or via a custom userAlphabet dictionary. <br><br>The `complexityType` Defines the complexity measure being employed. Three different complexity measures are provided by localCIDER, where the measure being used is passed via a string with one of 'WF', 'LC', or 'LZW'. WF is Wooton-Federhen complexity [8], which reports on the sequence's local Shannon entropy, and is the complexity measure used in the SEG algorithm. LC is Linguistic complexity [9], which reports on the number of distinct subsequences over the maximum number of different subsequences given the alphabet size and the word size. Finally, LZW is Lempel-Ziv-Welch [10] complexity, and effectivly asks how efficienctly the sequence can undergo lossless compression using unique subsequences. <br><br>The `alphabetSize` defines the size of the alphabet being used, where pre-defined alphabets are then used based on the specific size. Those pre-defined alphabets are defined below this table for clarity. By default an `alphabetSize` of 20 is used (i.e. no reduction in amino acid complexity). 'userAlphabet' Allows the user to define their own reduced alphabet. The format here is a dictionary where each key-value pair is amino-acid to X. This means you need a dictionary of length 20 where each amino acid is mapped to another amino acid. This is somewhat of tedious, but it helps avoid user-error where specific amino acids are missed. (default=None). <br><br>The `blobLen` is the sliding window size over which complexity is calculated (default=10). The `stepSize` is the size of steps taken as we define a new sliding window. Finally, the `wordSize` keyword is only relevant for the linguistic complexity (LC), ignored for other types, and defines the size of word for the algorithm. Default is 3. We recommend further reading of the associated literature to better understand these complexity measures. 
`get_linear_sequence_composition(blobLen=5, grps=[])`  | Returns an X by n matrix, where n is the number of amino acids in the sequence and X is the number of distinct groups of amino acids provided. This functions allows the local amino acid composition to be explored using a sliding window that computes the local density of one or more groups of amino acids along the sequence. The `blobLen` keyword defines the window size used to calculate the local sequence composition, and a stepsise of 1 is always used. The `grps` variable defines a list of lists, where each sub-list has elements that are amino acids. If no groups are provided a default grouping of amino acids by their physiochemical properties are used, i.e. `grps = [['E','D'], ['R','K'], ['R','K','E','D'], ['Q','N','S','T','G','H', 'C'], ['A','L','M','I','V'], ['F','Y','W'], ['P']]`, where groupsings are negativly charged, positivly chagred, charged, polar, aliphatic, aromatic, proline.
`get_Omega_sequence()` | Returns the reduced-alphabet sequence used to calculate the Omega parameter, where by R/K/D/E/P residues are defined as 'X' and all other residues are 'O'. This provides a clear visual description of the charge/proline patterning.


#### Reduced alphabets
Predefined alphabets shown below - all except eleven are based on alphabets defined in the reference below [11].

    two      - [(LVIMCAGSTPFYW), (EDNQKRH)]
    three    - [(LVIMCAGSTP), (FYW), (EDNQKRH)]
    four     - [(LVIMC), (AGSTP), (FYW), (EDNQKRH)]
    five     - [(LVIMC), (ASGTP), (FYW), (EDNQ), (KRH)]
    six      - [(LVIM), (ASGT), (PHC), (FYW), (EDNQ), (KR)]
    eight    - [(LVIMC), (AG), (ST), (P), (FYW), (EDNQ), (KR), (H)]
    ten      - [(LVIM), (C), (A), (G), (ST), (P), (FYW), (EDNQ), (KR), (H)]
    eleven   - [(LVIM), (C), (A), (G), (ST), (P), (FYW), (ED), (NQ), (KR), (H)]
    twelve   - [(LVIM), (C), (A), (G), (ST), (P), (FY), (W), (EQ), (DN), (KR), (H)]
    fifteen  - [(LVIM), (C), (A), (G), (S), (T), (P), (FY), (W), (E), (Q), (D), (N), (KR), (H)]
    eighteen - [(LM), (VI), (C), (A), (G), (S), (T), (P), (F), (Y), (W), (E), (D), (N), (Q), (K), (R), (H)]
    twenty   - all twenty!



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
`get_full_phosphostatus_kappa_distribution()`|  For each possible phospho-permutant, given the available phosphosites, calculate the kappa, fraction positive, fraction negative, FCR, NCPR, mean hydropathy and the phosphostatus. Phosphostatus is itself a tuple, where each position defines the phosphorylation status of each consecutive phosphosite in the sequence; 0 indicates unphosphorylated and 1 indicates phosphorylated. <br> <br> The function call returns a list of tuples, where each tuple contains the information described for each unique phospho-permutant. The assumption currently being made is that phosphorylation introduces a -1 charge only. While this may not be fully accurate, it provides a good and simple first approximation for the effect of phosphorylation.  Also note that as the number of phosphosites increases the number of calculations here scales as n^2. Be warned!


### Miscellaneous functions
The functions below represent a variety of miscellaneous functions.

Function name | Operation 
:---: | :---: 
`get_HTMLColorString()` <br> | Returns a fully formated HTML string which can be used to represent your sequence. The coloring used has a default, but can be defined using the `set_HTMLColorResiduePalette` function
`set_HTMLColorResiduePalette(colorDictionary)` <br> | Allows you to custom define a colour pallete. The colorDictionary must be a dictionary object that maps each of the 20 amino acids to a color. Currently 17 possible colors can be assigned to the 20 amino acids. These are;<br> aqua, black, blue, fuchsia, gray, green, lime, maroon, navy, olive, orange, purple, red, silver, teal, white, and yellow. This set of 17 colors represents the HTML browser compatible set of colors.

### Sequence permutation functions
These functions perform some operation on the sequence, returning a permuted SequenceParameter object populated with a different sequence. The underlying sequence object the function is called on is not altered. 

Function name | Operation 
:---: | :---: 
`get_shuffle(frozen=None)` <br> | Returns a SequenceParameter object with the primary amino acid sequence shuffled. If residue index positions are past as a list to 'frozen' those residues are considered imutable and are not shuffled. 

### Plotting functions (on-screen 'show' functions)
The following functions let you plot parameters from your sequence and display the results immediately on screen.

Function name | Operation 
:---: | :---: 
`show_phaseDiagramPlot(label="", title='Diagram of states', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)`| Renders a `matplotlib` Das-Pappu diagram of states plot with your sequence on the diagram<sup>[2]</sup>. If a `label` is provided this is a string which annotates your sequence on the plot. If a `title` is provided this sets the plot title. `legendOn` defines if the region labels are included as a legend. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. `getFig` defines if a matplotlib object is returned instead of being rendered on screen. 
`show_uverskyDiagramPlot(label="", title='Uversky Plot', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)`| Renders a `matplotlib` Uversky plot with your sequence on the diagram<sup>[4]</sup>. `label` can be a string which labels your sequence on the plot. If a `title` is provided this sets the plot title. `legendOn` defines if the regions labels are included as a legend. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. `getFig` defines if a matplotlib object is returned instead of being rendered.
`show_linearHydropathy(blobLen=5, getFig=False)`| Renders a `matplotlib` plot of the moving average hydropathy along the sequence, where the hydropathy is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used. `getFig` defines if a matplotlib object is returned instead of being rendered.
`show_linearNCPR(blobLen=5, getFig=False)` | Renders a `matplotlib` plot of the moving average net charge per residue (NCPR) along the sequence, where the NCPR is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used. `getFig` defines if a matplotlib object is returned instead of being rendered.
`show_linearFCR(blobLen=5, getFig=False)` | Renders a `matplotlib` plot of the moving average fraction of charged residues (FCR) along the sequence, where the FCR is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used. `getFig` defines if a matplotlib object is returned instead of being rendered.
`show_linearSigma(blobLen=5, getFig=False)` | Renders a `matplotlib` plot of the moving average sigma parameter along the sequence, where sigma is calculated in overlapping windows of size `blobLen`. Typically a blob length of 5-7 is used. Recall that sigma is calculated as the NCPR<sup>2</sup> / FCR. `getFig` defines if a matplotlib object is returned instead of being rendered.
`show_linearComplexity(complexityType='WF', alphabetSize=20, userAlphabet={}, windowSize=10, stepSize=1, wordSize=3, getFig=False)` | Renders a `matplotlib` plot of the linear sequence complexity. For a discussion of the various options see the `get_linear_complexity` description under the SequenceParameters functions table. `getFig` defines if a matplotlib object is returned instead of being rendered.


### Plotting functions (file-creating 'save' functions)
The following functions let you plot parameters from your sequence and save those plots to file for future use.

Function name | Operation 
:---: | :---: 
`save_phaseDiagramPlot(filename, label='', title='Diagram of states', legendOn=True, xLim=1, yLim=1, fontSize=10, saveFormat='png')` | Generates a `matplotlib` Das-Pappu diagram of states plot which is then saved to disk. `filename` is required and defines the file to be saved. Adding extensions is recommended but not required. All options are the same as in `show_phaseDiagramPlot`, with the addition of the `saveFormat` keyword, which defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz. (DEFAULT = png) 
`save_uverskyPlot(filename, label='', title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10, saveFormat='png')` | Generates a `matplotlib` Uversky plot with your sequence on the diagram<sup>[4]</sup> which is then saved to disk. `filename` is required and defines the file to be saved. Adding extensions is recommended but not required. All options are the same as in `show_uverskyPlot`, with the addition of the `saveFormat` keyword, which defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz. (DEFAULT = png) 
`save_linearHydropathy(filename, blobLen=5, saveFormat='png')` | Renders a `matplotlib` plot of the moving average hydropathy along the sequence, where the hydropathy is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
`save_linearNCPR(filename, blobLen=5, saveFormat='png')`| Renders a `matplotlib` plot of the moving average net charge per residue (NCPR) along the sequence, where the NCPR is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
`save_linearFCR(filename, blobLen=5, saveFormat='png')`| Renders a `matplotlib` plot of the moving average fraction of charged residues (FCR) along the sequence, where the FCR is calculated in overlapping windows of size `blobLen`. Typically 5-7 is used. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
`save_linearSigma(filename, blobLen=5, saveFormat='png')`  | Renders a `matplotlib` plot of the moving sigma value, where sigma defines the local charge assymetry and is used in the calculation of kappa. Sigma is calculated over blobs of `blobLen` size, typically with blobs of 5-7 residues. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
`save_linearComplexity(filename, complexityType='WF', alphabetSize=20, userAlphabet={}, windowSize=10, stepSize=1, wordSize=3, saveFormat='png')` | Renders a `matplotlib` plot of the linear sequence complexity. For a discussion of the various options see the `get_linear_complexity` description under the SequenceParameters functions table. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs savefig command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.
`save_linearCoposition(filename, blobLen=5, saveFormat='png', title='', plot_data=False))` | Renders a `matplotlib` plot of the local, position-specific amino acid composition. In version 0.1.9 the only grouping is the standard physiochemical grouping of amino acids, but in future versions we plan to add customizable groups to the plotting functions (customizable groups are available in the equivalent analysis function `get_linear_sequence_composition`). The local density is initially calculated and the fit to a univariate spline to remove noise and make the local sequence features more easily identifiable. If the `plot_data` variable is set to true the raw data is plotted alongside this spline fit, to ensure the fitting procedure is capturing the relevant sequence features. The plot is saved in the `filename` location. Adding extensions is recommended but not required. The `saveFormat` keyword defines the output format - this parameter is passed to matplotlibs `savefig` command which supports the following filetypes: emf, eps, pdf, png, ps, raw, rgba, svg, svgz.



# Plotting with the plots module

The `plots` module lets you generate or save plots without explicitly invoking sequenceParameters objects. Such plots can either be generated using the raw values, or by passing in sequenceParameters objects.

Unlike the `SequenceParameters`, `plots` is a collection of functions, and does not allow for the creation of `Plots` objects.

In the following discussion a **single** plot contains information from just one sequence, whereas a **multiple** plot contains sequence information from many sequences

### Using `Plots` functions

There are two modes with which you can use plots functions;

#### Passing raw parameters to be plotted (single and multiple plots)
This may be convenient if you're plotting previously calculated values, or plotting a theoretical distribution here. Either individual values or lists of values can be passed.

### Passing SequenceParameters objects (multiple plots only)
If you're analyzing a set of SequenceParameters objects you may wish to place them on the same plot. Note we do not offer this functionality for single sequences, because you can always plot a single sequence using the associated `SequenceParameters` function (e.g. `.show_uverskyPlot()` or `.save_phaseDiagramPlot()`).


Function name | Input Arguments | Discussion
:---: | :---: | :---:
`show_single_phasePlot(fp, fn, label='', title='Diagram of states', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)` | `fp` is the fraction of positive residues in a sequence, `fn` is the fraction of negative residues in a sequence. `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. `getFig` defines if a matplotlib object is returned instead of being rendered.   | No programatic output is generated, but the plot is rendered on the screen. If `getFig` is set to `True` then the matplotlib object is returned, rather than rendered on the screen.
`show_single_uverskyPlot(hydropathy, mean_net_charge, label='', title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)` | `hydropathy` is the sequence's 0-1 normalized hydropathy, while `mean_net_charge` is the the absolute value of the sequence's mean net charge. `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. `getFig` defines if a matplotlib object is returned instead of being rendered. | No programatic output is generated, but the plot is rendered on the screen. If `getFig` is set to `True` then the matplotlib object is returned, rather than rendered on the screen.
`save_single_phasePlot(fp, fn, filename, label='', title='Diagram of states', legendOn=True, xLim=1, yLim=1, fontSize=10)` | `fp` is the fraction of positive residues in a sequence, `fn` is the fraction of negative residues in a sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_single_uverskyPlot(hydropathy, mean_net_charge, filename, label='', title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10)` | `hydropathy` is the sequence's 0-1 normalized hydropathy, while `mean_net_charge` is the the absolute value of the sequence's mean net charge. `filename` is the name of the file to save the plot as (`.png` is appended). `label`, if included, should be a string with which the sequence on the plot is labelled. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. | No programatic output is generated, but the plot is constructed and saved at the `filename` location.

If you wish to plot a single sequence's Das-Pappu diagram of states plot or Uversky plot from a `SequenceParameters` object we suggest using the appropriate `SequenceParameters` functions (e.g.  `.show_uverskyDiagram()`, `.show_phaseDiagram()`)

### `*_multiple*_` functions

For the `multiple` plotting functions, we pass in a list of values or objects corresponding to the multiple sequences of interest. 

There are two types of multiple plots

**1)** Those which you pass two lists of the actual values of interest, where the lists represent parallel lists, where the position in each list corresponds to a specific sequence.
**2)** Those which you pass a list of `SequenceParameters` objects to, meaning you only pass it a single list


The `show_*` functions render a plot on the screen in front of you through `matplotlib`. These are useful for directly interacting with sequences and examining how specific sequence might behave. You can also directly access the matplotlib object through the `getFig` argument.

The `save_*` functions save a plot to the designated file name.

Function name | Input Arguments | Discussion
:---: | :---: | :---:
`show_multiple_phasePlot(fp_list, fn_list, label_list=[""], title="Diagram of states", legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)`  | `fp_list` and `fn_list` are lists of the positive fraction and negative fraction of sequences, where the list position in each of these two 'parallel' lists corresponds to the same sequence. That is to say, if you were looking at 5 sequences, `fp_list` and `fn_list` would both be of length 5, and `fp_list[2]` and `fn_list[2]` would represent the fraction position and fraction negative of your third sequence, respectively. `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. If `getFig` is set to `True` then the matplotlib object is returned, rather than rendered on the screen. | No programatic output is generated, but the plot is rendered to the screen.  
`show_multiple_phasePlot2(SeqParam_list, label_list=[""], title='Diagram of states', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)` | All arguments are the same as above, except now instead of explicit fraction of positive and fraction of negative lists, we have a single list containing `SequenceParameters` objects, from which these properties are derived.  | No programatic output is generated, but the plot is rendered to the screen.
`show_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, label_list=[""], title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)` | `hydropathy_list` and `mean_net_charge_list` are lists of the 0-1 normalized hydropathy and the absolute mean net charge, where the list position in each of these two 'parallel' lists corresponds to the same sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the max values for the X and Y axes. `fontSize` defines the size of the label font. If `getFig` is set to `True` then the matplotlib object is returned, rather than rendered on the screen. | No programatic output is generated, but the plot is rendered to the screen. 
`show_multiple_uverskyPlot2(SeqParam_list, label_list=[""], title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10, getFig=False)` | All arguments are the same as above, except now instead of having an explicit hydropathy list and an explicit mean net charge list, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is rendered to the screen. 
`save_multiple_phasePlot(fp_list, fn_list, filename, ?label_list=[""], ?title="Diagram of states", ?legendOn=True, ?xLim=1, ?yLim=1,?fontSize=10)`  | `fp_list` and `fn_list` are lists of the positive fraction and negative fraction of sequences, where the list position in each of these two 'parallel' lists corresponds to the same sequence. That is to say, if you were looking at 5 sequences, `fp_list` and `fn_list` would both be of length 5, and `fp_list[2]` and `fn_list[2]` would represent the fraction position and fraction negative of your third sequence, respectively. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the maximum values for the fraction of negative and positive residues. `fontSize` defines the size of the label font. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_multiple_phasePlot2(SeqParam_list, filename, label_list=[""], title='Diagram of states', legendOn=True, ?xLim=1, ?yLim=1,?fontSize=10)` | All arguments are the same as above, except now instead of explicit fraction of positive and fraction of negative lists, we have a single list containing `SequenceParameters` objects, from which these properties are derived.| No programatic output is generated, but the plot is constructed and saved at the `filename` location.
`save_multiple_uverskyPlot(hydropathy_list, mean_net_charge_list, filename, label_list=[""], title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10`) | `hydropathy_list` and `mean_net_charge_list` are lists of the 0-1 normalized hydropathy and the absolute mean net charge, where the list position in each of these two 'parallel' lists corresponds to the same sequence. `filename` is the name of the file to save the plot as (`.png` is appended). `label_list`, if included, should be a list of strings corresponding to the length of `fn_list` and `fp_list`, and includes names for each sequence on the plot. `title`, if provided, sets the plot title. `legendOn` defines if the 5-region legend for the phase diagram plot is included. `xLim` and `yLim` define the maximum values for the fraction of negative and positive residues. `fontSize` defines the size of the label font.| No programatic output is generated, but the plot is constructed and saved at the `filename` location. 
`save_multiple_uverskyPlot2(SeqParam_list, filename, label_list=[""], title='Uversky plot', legendOn=True, xLim=1, yLim=1, fontSize=10)` | All arguments are the same as above, except now instead of having an explicit hydropathy list and an explicit mean net charge list, we have a single list containing `SequenceParameters` objects, from which these properties are derived. | No programatic output is generated, but the plot is constructed and saved at the `filename` location. 


# `localcider` exceptions

If you're into exception handling, localcider has a set of defined exceptions (all of which inherits from the `Exception` base class) which can be access by importing localcider - i.e.

    import localcider
    
    localcider.SequenceFileParserException
    localcider.SequenceException
    localcider.ResTableException
    localcider.PlottingException
    
There are a couple more but these relate to functionality not yet available in version 0.1.9

# Tutorials 


## 1) Calculate a bunch of parameters associated with a sequence I found on the internet

    from localcider.sequenceParameters import SequenceParameters
    
    # define my internet sequence and create a new SequenceParameters object with it
    internetSequence="VAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP"
    object=SequenceParameters(internetSequencexs)
    
    
    # Let's get some parameters!
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
    
    # populate that list with SequenceParameters objects, which we construct from the
    # sequence found in each file
    for file in filelist:
       list_of_SeqObjs.append(SequenceParameters(sequenceFile=file))

	# for each 
	for obj in list_of_SeqObjs:
	   print obj.get_kappa()
	   
Uh oh! One of our sequence files has an amino acid "X" in it - this isn't a canonical amino acid, so triggers a `SequenceFileParserException`. We actually just want to skip any such files, so we change the code above slightly. Specifically, let's change the for-loop that loops over each file;

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


## Example 3 - plotting two sequence on a single Das-Pappu diagram of states

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


## FAQ

1. **I get a kappa value greater than 1?!**

kappa is a ratio of your sequence's delta over the maximium possible value for a sequence of that composition. The original mechanism for computing what the maximum delta value is was using Monte Carlo simulations, but this is *extremely* computationally expensive. To avoid this we have developed a new deterministic algorithm that uses a number of rule-based heuristics. During early testing of these heuristics there were examples of kappa values being greater than 1, which indicates that these heuristics have failed. As of Fall 2015 we have run many millions of disordered sequences without encountering sequences with a kappa > 1. However, should you encounter a sequence with a kappa value greater than 1 PLEASE email the sequence to me, as it represents a very rare edge case that we haven't encountered in our testing procedures.

2. **I want to run localCIDER on my cluster but it has not graphical frontend and I can't get matplotlib to work**

There is a plan to release a higher performance version of localCIDER which carries out the core functionality using C++. However, that's a way off, so for now a quick hack is to do the following;

    * Download the sourcecode and unpack - this will be a local package, and means you don't have to instal localCIDER globally
    * Move your local package to where you want to carry out the analysis
    * Edit localcider/backend/plotting.py and comment out the two import lines which reference matplotlib
    * Run your analysis

If you call any functions which depend on `matplotlib` this will obviously fail, but if not it'll work just fine.

# Epilog

## Feature request
We'd like localCIDER to be the single-software solution for sequence analysis. To this end, if you have analysis routines you would like to be incorporated used please contact us about adding these into localCIDER. If you can provide us with Python code for performing the analysis in question, we can almost guarentee incorporation into the next version of localCIDER. 

## About
localCIDER was developed in the Pappu lab by Alex Holehouse, with additional code by James Ahad and Mary Richardson. Many of the early ideas were pioneered by [Dr. Rahul Das](https://sites.google.com/site/rahulkdas82/home). 

We believe scientific software should be held to same (or higher) software standards as general commercial and open source projects. To this end we use Git for version control, employ standard conventions for versioning, run an ever-expanding set of unit tests, and try to ensure all code is well structured and well documented. The source code is fully open and [can be viewed here](https://github.com/Pappulab/localCIDER).


## References

**[1]** Campen, A. et al. TOP-IDP-scale: a new amino acid scale measuring propensity for intrinsic disorder. Protein Pept. Lett. 15, 956–963 (2008).

**[2]** Das, R. K. & Pappu, R. V. Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues. Proc. Natl. Acad. Sci. U. S. A. 110, 13392–13397 (2013).

**[3]** Kyte, J. & Doolittle, R. F. A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157, 105–132 (1982).

**[4]** Uversky, V. N. Natively unfolded proteins: a point where biology waits for physics. Protein Sci. 11, 739–756 (2002).

**[5]** Mao, A. H., Crick, S. L., Vitalis, A., Chicoine, C. L. & Pappu, R. V. Net charge per residue modulates conformational ensembles of intrinsically disordered proteins. Proc. Natl. Acad. Sci. U. S. A. 107, 8183–8188 (2010).

**[6]** Elam WA, Schrank TP, Campagnolo AJ, Hilser VJ. Evolutionary conservation of the polyproline II conformation surrounding intrinsically disordered phosphorylation sites. Protein Sci. 2013; 22: 405- 417. doi: 10.1002/pro.2217 PMID: 23341186

**[7]** Tomasso, M. E., Tarver, M. J., Devarajan, D. & Whitten, S. T. Hydrodynamic Radii of Intrinsically Disordered Proteins Determined from Experimental Polyproline II Propensities. PLoS Comput. Biol. 12, e1004686 (2016).

**[8]** Wootton, J. C., & Federhen, S. (1993). Statistics of local complexity in amino acid sequences and sequence databases. Computers & Chemistry, 17(2), 149-163.

**[9]** Troyanskaya, O. G., Arbell, O., Koren, Y., Landau, G. M., & Bolshoy, A. (2002). Sequence complexity profiles of prokaryotic genomic sequences: a fast algorithm for calculating linguistic complexity. Bioinformatics , 18(5), 679 - 688.

**[10]** Lempel, A.& Ziv, J. (1976). On the complexity of finite sequence. IEEE Trans. Inf. Theory, vol. IT-22, no. 1, 75-81.

**[11]** Murphy, L. R., Wallqvist, A., & Levy, R. M. (2000). Simplified amino acid alphabets for protein fold recognition and implications for folding. Protein Engineering, 13(3), 149-152.

**[12]** Rucker, A.L., Pager, C.T., Campbell, M.N., Qualls, J.E., and Creamer, T.P. (2003). Host-guest scale of left-handed polyproline II helix formation. Proteins 53, 68-75.

**[13]** Shi, Z., Chen, K., Liu, Z., Ng, A., Bracken, W.C., and Kallenbach, N.R. (2005). Polyproline II propensities from GGXGG peptides reveal an anticorrelation with beta-sheet scales. Proc. Natl. Acad. Sci. U. S. A. 102, 17964-17968.

**[14]** Martin, E.W., Holehouse, A.S., Grace, C.R., Hughes, A., Pappu, R.V., and Mittag, T. (2016). Sequence determinants of the conformational properties of an intrinsically disordered protein prior to and upon multisite phosphorylation. J. Am. Chem. Soc. 138, 15323-15335.
        
        
## Thanks
Many people have been involved in this project. We'll try and include an up-to-date list here;

* Paul Nobrega (University of Massachusetts) for bug reports and the Windows installation instructions
* Davide Mercadante (Heidelberger Institut fur Theoretische Studien) for bug reports and code development to improve plotting customization
* Katra Kolsek (Heidelberger Institut fur Theoretische Studien) for bug reports and code development to improve plotting customization
* Alex Chin (Johns Hopkins University) for a crucial bug report
* Thomas Pranzatelli (Washington University in St. Louis) for testing and bug reports
* Carlos Hernández (Stanford University) for Python 3 support and PEP8 compliance 
* Luke Wheeler (University of Oregon) for Python 3 testing
* Xiaohan Li (Yale University) for corrections to text
* Sean Cascarina (Colorado State University) for finding a bug where stop-codons are not dealt with correctly in FASTA files
* David Sanders and Anastasia Repouliou (Princeton University) for the suggestion of introducing pH dependent charge analysis
 

## Update schedule

* **version 0.1.0** - August 10th 2014: Initial release 
* **version 0.1.1** - August 14th 2014: Minor bug fixes, more tests, improved sequence parsing
* **version 0.1.2** - October 1st 2014: More tests, fixed heuristic bugs in extreme kappa sequences, significant improvement on plotting customization, cleaned up docs and code documentation
* **version 0.1.3** - October 13th 2014: More tests, corrected a bug introduced in 0.1.2 for sequences with no neutral residues, corrected a bug in how titles/labels were assigned by default, improved how classification into the 5 regions is done, and explicitly allow for regions 4 (positive electrolytes) and 5 (negative electrolytes) rather than general polyelectrolytes
* **version 0.1.4** - January 16th 2015: Corrected bug in the way KD hydrophobicity was calculated - alanine incorrectly set to 0
* **version 0.1.7** - February 25th 2015: Bunch of additional updates and code readability improvements, as well as changing how NCPR/FCR plots are made (bar instead of line)
* **version 0.1.8** - March 28th 2016: Major update with a number of additional features. All figures can now be generated as PDFs, sequence complexity analysis has been added, PPII propensity added, kappa-P added. Pleasantly, no bugs needed fixing, however!
* **version 0.1.9** - September 31st 2016: Added local amino acid composition, general patterning parameter, and improved plot formatting 
* **version 0.1.10** - Explicit sequence shuffling with residue freezing, Omega code update, additional PPII scales, updated references, updated plotting functions, font sizes fixed, tests updated, (dynamic rescaling of font size and line-widths)
* **version 0.1.11** - Update to colors for amino acid output string, added the get_linear_sigma function, moved the import location for the scipy dependency so it can 
* **version 0.1.12** - Deal with stop codons when reading fasta files. Introduced pH dependence analysis. 




----------
        <ul>
          <li><a href="https://github.com/Pappulab/localCIDER">View On <strong>GitHub</strong></a></li>
          <li><a href="http://pappulab.wustl.edu">Go to <strong>Pappu lab</strong></a></li>
          <li><a href="http://pappulab.wustl.edu/CIDER/">Go to <strong>CIDER</strong> </a></li>
        </ul>



