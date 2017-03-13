==========
localCIDER
==========

`version 0.1.12 - February 2017`

**Introduction**

**localCIDER** is a Python package developed by the lab of `Rohit Pappu at Washington University in St. Louis <http://pappulab.wustl.edu>`_ for calculating and plotting parameters associated with intrinsically disordered proteins (IDPs) and disodered regions (IDRs). **localCIDER** is the Python backend for `CIDER <http://pappulab.wustl.edu/CIDER.html>`_, ( **C**\lassification of **I**\ntrinsically **D**\isordered **E**\nsemble **R**\egions) - a webserver for the calculation of many of those same properties. Essentially, localCIDER lets you run CIDER's calculations locally, allowing you to create custom analysis pipelines which do not rely on the webserver. It also allows you to take advantage of your own local computing hardware, rather than competing with everyone else for a common set of hardware provided by the Pappu lab.

This project was motivated by the need to rapidly and easily calculate the |kgr| (kappa) parameter, as defined in the 2013 Das & Pappu PNAS paper [1], as well as provide a tool to easily plot a sequence on the Das-Pappu diagram-of-states;

.. image :: http://pappulab.wustl.edu/img/phase_diagram_small.png
   :width: 400 px

The original method for calculating |kgr| involved an incredibly computationally expensive Monte Carlo step (used to determine the maximally segregated sequence), whereas localCIDER makes use of a newly developed algorithm that computes |kgr| in O(1) time. More generally, localCIDER provides a wide range of sequence analysis routines for understanding, classifying, and designing disordered proteins. Beyond charge-based parameters, a wide range of additional parameters can also be calculated. **localCIDER** represents an ideal tool for the construction of a high-throughput sequence analysis pipelines, and we believe provides a direct route to explore the sequence-to-ensemble relationship of IDPs.

For more information please `see the full documentation <http://pappulab.github.io/localCIDER/>`_ .

------------------
Installation
------------------

To install run

    [sudo] pip install localcider

We strongly recommend using `virtualenv`, in which case `sudo` should not be required. Note that localcider requires `numpy`, `scipy`, and `matplotlib` to run.

------------------------------
Usage, bugs, and questions
------------------------------

Please see the `see the full documentation <http://pappulab.github.io/localCIDER/>`_ for usage guidelines. Please address all questions and bug reports to `Alex <http://pappulab.wustl.edu/people.html#grads>`_ and he'll do his best to get back to you!

----------
About
----------

**localCIDER** was written by Alex Holehouse and James Ahad in the `Pappu Lab <http://pappulab.wustl.edu>`_ . Please cite localCIDER (and CIDER) as:

A.S. Holehouse, R.K. Das, J.N. Ahad, M.O.G. Richardson, R.V. Pappu (2017) CIDER: Resources to analyze sequence-ensemble relationships of intrinsically disordered proteins. Biophysical Journal, 112: 16-21

Also, make sure to cite the various studies used if you use analyses that make use of other datasets (e.g. PPII propensity, charge patterning, proline patterning, hydrophobicity etc.).

--------------
References 
--------------

[1] `Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues <http://www.pnas.org/content/110/33/13392.short>`_ R.K. Das & R.V. Pappu (2013) PNAS **110**\, 33, pp13392–13397.


.. |kgr|  unicode:: U+003BA .. GREEK SMALL LETTER KAPPA
