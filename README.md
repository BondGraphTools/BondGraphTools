[![PyPI version](https://badge.fury.io/py/BondGraphTools.svg)](https://badge.fury.io/py/BondGraphTools)
[![Build Status](https://travis-ci.com/BondGraphTools/BondGraphTools.svg?branch=master)](https://travis-ci.com/BondGraphTools/BondGraphTools)
[![Test Coverage](https://api.codeclimate.com/v1/badges/4735c13a87b24d3a1899/test_coverage)](https://codeclimate.com/github/BondGraphTools/BondGraphTools/test_coverage)
# BondGraphTools - A Toolkit for modelling multi-physics systems.
## Summary

This toolkit is for rapid modelling and design of networked phsyical systems.
It is conceptually based upon the Bond Graph modelling methodology.

## Documentation

https://bondgraphtools.readthedocs.io/

## Installation

### Dependencies
BondGraphTools requires:
- python 3.7 or greater,
- LaTeX,
- gFortran, BLAS, LAPACK and SUNDIALS (for the numerical solvers).

Python dependencies:
- sympy
- numpy
- scipy
- matplotlib
- scikits.odes 

## Instructions:

#### Recommended installation using Conda
1. Install [Anaconda](https://anaconda.org/) (or miniconda).
2. Install dependencies from conda-forge: liblas, liblapack, sundials using
   `conda install -c conda-forge libblas liblapack sundials fortran-compiler`
3. Install LaTeX via `conda install -c conda-forge miktex` (windows) or `conda install -c conda-forge texlive-core` (osx/linux)
4. Install the package using `pip install bondgraphtools`
