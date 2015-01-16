==========
localCIDER
==========

`version 0.1.4 - January 2015`

**Introduction**

**localCIDER** is a Python package developed by the lab of `Rohit Pappu at Washington University in St. Louis <http://pappulab.wustl.edu>`_ for calculating and plotting parameters associated with intrinsically disordered proteins (IDPs) and disodered regions (IDRs). **localCIDER** is the Python backend for `CIDER <http://pappulab.wustl.edu/CIDER.html>`_, ( **C**\lassification of **I**\ntrinsically **D**\isordered **E**\nsemble **R**\egions) - a webserver currently in the final stages of development, also by the Pappu lab. Essentially, localCIDER lets you run CIDER's calculations locally, allowing you to create custom analysis pipelines which do not rely on the webserver. It also allows you to take advantage of your own local infrastructure, rather than competing with everyone else for a common set of hardware provided by the Pappu lab.

This project was motivated by the need to rapidly and easily calculate the |kgr| (kappa) parameter, as defined in the 2013 Das & Pappu PNAS paper [1], as well as provide a tool to easily plot a sequence on the Pappu-Das diagram of states;

.. image :: http://pappulab.wustl.edu/img/phase_diagram.png
   :width: 500 px

The original method for calculating |kgr| involved an incredibly computationally expensive Monte Carlo step (used to determine the maximally segregated sequence) , whereas localCIDER makes use of a number of heuristics to reduce the Monte Carlo search space to an optimized and highly manageable size (most |kgr| calculations are carried out in less than a second).

Beyond |kgr|, **localCIDER** lets you calculate a range of additional parameters, as well as create the Pappu-Das phase diagram (shown below) and a number of other plots.

For more information please `see the full documentation <http://pappulab.github.io/localCIDER/>`_ .

------------------
Installation
------------------

To install run

    [sudo] pip install localcider

Note that localcider requires `numpy` and `matplotlib` to run.

------------------------------
Usage, bugs, and questions
------------------------------

Please see the `see the full documentation <http://pappulab.github.io/localCIDER/>`_ for usage guidelines. Please address all questions and bug reports to `Alex <http://pappulab.wustl.edu/people.html#grads>`_ and he'll do his best to get back to you. 

----------
About
----------

**localCIDER** was written by Alex Holehouse and James Ahad in the `Pappu Lab <http://pappulab.wustl.edu>`_ . A manuscript is currently under preparation for citation, but until that time please cite localCIDER as;

    A.S. Holehouse, J. Ahad, R.K. Das, R.V. Pappu (2014) localCIDER, http://pappulab.wustl.edu/localCIDER/

--------------
References 
--------------

[1] `Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues <http://www.pnas.org/content/110/33/13392.short>`_ R.K. Das & R.V. Pappu (2013) PNAS **110**\, 33, pp13392–13397.


.. |kgr|  unicode:: U+003BA .. GREEK SMALL LETTER KAPPA
