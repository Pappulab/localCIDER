localCIDER
==========

localCIDER is the Python backend for CIDER developed by the Pappu lab at Washington University in St. Louis.

For more information please
[see the documentation](http://pappulab.github.io/localCIDER/).

localCIDER was written by Alex Holehouse and James Ahad in the [Pappu Lab](http://pappulab.wustl.edu/).

While Alex has since left the Pappu lab, he continues to maintain localCIDER - for now, please address all questions [to him](http://www.holehouse.wustl.edu).


Nardini
=======

### NB: Nardini installation required
To use nardini analysis described below, one must install nardini 'manually' - i.e. as of version 0.1.21 localCIDER does not enforce a hard dependency of nardini during installation. Nardini can be installed by running:

	pip install nardini 
	
### Nardini analysis

To perform Nardini-related analysis via a localCIDER SequenceParameter object, the following example is illustrative:

```
from localcider.sequenceParameters import SequenceParameters

# create a SequenceParameters object with your amino acid sequence
SeqObj = SequenceParameters('YOURSEQUENCEHERE')

# `num_scrambles` and `random_seed` are optional arguments. Shown
# below are their default values.
SeqObj.save_zscoresAndPlots(num_scrambles=100000, random_seed=None)
```

An interface for setting a random seed is provided for reproducibility. If no random seed is provided, one will be chosen and reported to the user. As analysis proceeds, the progress will be printed to the user in the command-line or Jupyter notebook. Upon completion, the generated plots (PNGs) as well as a zip file containing TSVs of the parameters as well as matrix representations of the plots will be exported.

The function `save_zscoresAndPlots` encapsulates several methods from the [Nardini package](https://github.com/mshinn23/nardini) into a user-friendly entry point for subsequent analysis. A complementary shell-interface (`nardini`) is also made available when the Nardini package is installed which can be paired with this method. To learn more, run `nardini -h` from the command-line or visit the Github repository for [Nardini](https://github.com/mshinn23/nardini#command-line-usage) for more details.

In addition to `save_zscoresAndPlots`, LocalCIDER also integrates a couple other [Nardini-based](https://github.com/mshinn23/nardini) methods. All three methods are detailed below:

| LocalCIDER Nardini Function | Parameters | Description |
|:--:|:--:|--|
| `save_zscoresAndPlots` | `num_scrambles=100000, random_seed=None` | Calculates the Nardini z-score of the input sequences and generates a ZIP file containing the Nardini analysis (TSV and CSV) as well as accompanying plots. |
| `calculate_zscore` | `num_scrambles=100000, random_seed=None` | This function performs the Nardini z-score analysis and returns a 5-item dictionary containing: 1) the original sequence; 2) the scrambled sequence; 3) the sequence number (used for book-keeping for many sequences); 4) the `reshaped_zvecdb` corresponding to the original sequence; and, 5) the `reshaped_zvecdbscr` corresponding to the scrambled sequences. |
| `plot_nardini_zscores` | `seq_name`, `zvec_db`, `typeall`, `index`, `savename`, `is_scrambled` | This function generates a ZIP file containing plots of the Nardini analysis. Parameters: `seq_name` (the name of the sequence); `zvec_db` (the numpy array corresponding to the zscore-vector used for calculations); `typeall` (the amino acid types used for the analysis contained within `zvec_db`); `index` (the index of `seq_name` in the `zvec_db`); `savename` (the base filename under which to save the plots); and, `is_scrambled` (a boolean that indicates whether or not the matrix corresponds to the scrambled sequence analysis). |

It should be noted that:

1. `save_zscoreAndPlots` is a convenience function which performs and exports the Nardini analysis in a single step. This method call produces a ZIP file that will contain the file `sequences.tsv` (the sequences analyzed). For every sequence contained within the input file, 4 additional files will also be generated and included in the ZIP file: `regular-<seq_name>.png` (the plot of the Nardini matrix of the original sequence); `scrambled-<seq_name>.png` (the plot of the Nardini matrix of the closest matching scrambled sequence); `zscore-original-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the original sequence); and, `zscore-scrambled-sequence-<seq_name>.tsv` (the text file corresponding to the z-score matrix of the scrambled sequence). Contextually, this function encapsulates the functions of items 2 and 3 while also exporting them into a ZIP file.
2. `calculate_zscore` is meant for performing the Nardini analysis and using the output in another method. A Python dictionary is returned containing the relevant matrices and objects. No plots are generated.
3. `plot_nardini_zscores` is meant for plotting and saving images of the Nardini analysis (i.e. matrices) that were previously generated.


## Running tests
Before committing changes please run tests as

	cd localcider/tests
	python runtests.py

Note that subsets of tests can be run by editing `runtests.py` as needed.

Alternatively one can run

	pytest -vv

To run all tests
