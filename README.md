localCIDER
==========

localCIDER is the Python backend for CIDER developed by the Pappu lab at Washington University in St. Louis.

For more information please 
[see the documentation](http://pappulab.github.io/localCIDER/).

localCIDER was written by Alex Holehouse and James Ahad in the [Pappu Lab](http://pappulab.wustl.edu/).

While Alex has since left the Pappu lab, he continues to maintain localCIDER - for now, please address all questions [to him](http://www.holehouse.wustl.edu). 

## Running tests
Before committing changes please run tests as

	cd localcider/tests
	python runtests.py
	
Note that subsets of tests can be run by editing `runtests.py` as needed.

Alternatively one can run

	pytest -vv
	
To run all tests