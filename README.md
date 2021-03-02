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

BondGraph requires:
- python 3.6 
- sundials 5.*

Python dependencies:
- sympy
- numpy>=1.14
- scipy>=1.0.1
- matplotlib>=2.2.2
 
 LaTEX Dependencies (Used by matplotlib):
 - case for Windows: You can install MiKTex (here https://miktex.org/2.9/setup)
 when you star the draw() MikTex will automaticly try to install all required, if unsuccessfully try to use MiKTex console-->Packages-> and install required packages by searching.
 - for other OS should be similar.

### Instructions:
1. Install python > 3.6 for your operating system. 
2. Install sundials (either from https://computing.llnl.gov/projects/sundials/, or via Anaconda https://anaconda.org/conda-forge/sundials)
3. Install BondGraphTools using PyPI; `pip install BondGraphTools`

