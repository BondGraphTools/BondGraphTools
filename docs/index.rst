.. BondGraphTools documentation master file, created by
   sphinx-quickstart on Tue Sep 25 14:32:01 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:
   :maxdepth: 1

   tutorials
   guides
   discussion
   api

BondGraphTools
==============
Overview
^^^^^^^^
BondGraphTools is a python library for systems modelling based on the bond graph
methodology [Gaw1996]_, [Gaw2007]_.

BondGraphTools is intended to be used by:
   * software developers; as a framework to build modelling interfaces on top of
   * engineers; to quickly build and simulate physical systems
   * mathematicians; to perform model reduction and analysis

BondGraphTools is built upon the scientific python stack; numpy, scipy, sympy
and matplotlib, and relies on Julia and DifferentialEquations.jl to handle
numnerical simulation.

How to read the docs
^^^^^^^^^^^^^^^^^^^^
The documentation for this library is organised into five sections.
   * This page show how to install BondGraphTools
   * :doc:`tutorials` contains a list of step-by-step tutorials demonstrating how
     to use BondGraphTools.
   * :doc:`guides` contains recipes for common tasks.
   * :doc:`discussion` contains high level discussion about to library.
   * :doc:`api` is the library reference documentation.

Getting Started
^^^^^^^^^^^^^^^
Requirements
------------

* Python 3.6 or above.
* Julia 0.6. *Note: Julia 0.7 and above are not yet supported*

Installation
------------

1. Install Julia 0.6.4 which can be downloaded from here_.
2. Make sure Julia 0.6.4 is in your path variables. You can test this by typing
   ``julia -v`` at the command prompt. If this returns ``julia version 0.6.4``
   then julia is in your path variables.
3. Install BondGraphTools via pip: ``pip install BondGraphTools``

BondGraphTools is now installed and will download and install the required
julia packages the first time you import it.


Usage
-----

BondGraphTools can be loaded inside a jupyter-notebook, python interpreter, or
script by using::

    import BondGraphTools

Users need an internet connection during the first load as the it will download
the required packages from the julia package repository.

Api reference for the package and it's contents can be accessed using the help
command::

    help(BondGraphTools)

Contributing
------------
The best way to contribute right now is to use it!
If you wish to help out, please visit the project github_.

Bibliography
^^^^^^^^^^^^
.. _github: https://github.com/BondGraphTools/BondGraphTools
.. _here: https://julialang.org/downloads/oldreleases.html
.. [Gaw2007] https://ieeexplore.ieee.org/document/4140745
.. [Gaw1996] http://www.gawthrop.net/Books/GawSmi96.pdf
