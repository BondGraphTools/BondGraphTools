Getting Started
===============
Overview
--------
BondGraphTools is a physical systems modelling library based on the Bond Graph
technique. This library seeks to provide a set of user-friendly tools for
quickly building, analysing and simulating Bond Graphs.

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


Using BondGraphTools
--------------------

BondGraphTools can be loaded inside a jupyter-notebook, python interpreter, or
script by using::

    import BondGraphTools

Users need an internet connection during the first load as the it will download
the required packages from the julia package repository.

Api reference for the package and it's contents can be accessed using the help
command::

    help(BondGraphTools)

.. _here: https://julialang.org/downloads/oldreleases.html
