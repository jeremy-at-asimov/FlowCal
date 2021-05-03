===========
``FlowCal``
===========
``FlowCal`` is a library for processing and analyzing flow cytometry data in Python.
It features:

* Extraction of Flow Cytometry Standard (FCS) files into numpy array-like structures
* Traditional and non-standard gating, including automatic density-based two-dimensional gating.
* Traditional transformation functions, such as exponentiation.
* Analysis of calibration beads data, standard curve generation, and transformation to absolute units (Molecules of Equivalent Fluorophore, MEF).
* Plotting, including generation of histograms, density plots and scatter plots.
* A user-fiendly Excel UI to gate, transform, plot, and generate statistics from a list of flow cytometry samples in a simple fashion.

What we've added in the Asimov Fork:
* Saving and applying gating contours
* Binning by fluorescence

Documentation
=============
The official documentation for the Tabor lab repo can be found in https://flowcal.readthedocs.io. Info on additions made for the Asimov fork is currently in the source code, but may be added to a separate location eventually.

Report Bugs
===========
The official way to report a bug for the Tabor lab is through the issue tracker on github (https://github.com/taborlab/FlowCal/issues). Try to be as explicit as possible when describing your issue. Ideally, a set of instructions to reproduce the error should be provided, together with the version of all the relevant packages you are using.

Request Features
================
Features can also be requested through the issue tracker on github. Try to be as descriptive as possible about the desired feature.
