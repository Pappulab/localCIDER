localCIDER
==========



localCIDER is the Python backend for CIDER developed by the Pappu lab at Washington University in St. Louis.

For more information please
[see the documentation](http://pappulab.github.io/localCIDER/).

localCIDER was written by Alex Holehouse and James Ahad in the [Pappu Lab](http://pappulab.wustl.edu/).

While Alex has since left the Pappu lab, he continues to maintain localCIDER - for now, please address all questions [to him](http://www.holehouse.wustl.edu).


Nardini
=======

To perform Nardini-related analysis via a localCIDER SequenceParameter object, the following example is illustrative:

```
from localcider.sequenceParameters import SequenceParameters

# create a SequenceParameters object with your amino acid sequence
SeqObj = SequenceParameters('YOURSEQUENCEHERE')

# `num_scrambles` and `random_seed` are optional arguments. Shown
# below are their default values.
SeqObj.save_zscoresAndPlots(num_scrambles=100000, random_seed=None)
```

An interface to setting a random seed is provided for reproducibility. As analysis proceeds, a report will be printed to the user in the command-line or Jupyter notebook. Upon completion, the generated plots (PNGs) as well as a zip file containing TSVs of the parameters as well as matrix representations of the plots will be exported. The function `save_zscoresAndPlots` encapsulates several methods from the [Nardini package](https://github.com/mshinn23/nardini) into a user-friendly entry point for subsequent analysis. A complementary shell-interface (`nardini`) is also made available when the Nardini package is installed which can be paired with this method. To learn more, run `nardini -h` from the command-line or visit the Github repository for Nardini for more details.


## Running tests
Before committing changes please run tests as

	cd localcider/tests
	python runtests.py

Note that subsets of tests can be run by editing `runtests.py` as needed.

Alternatively one can run

	pytest -vv

To run all tests
