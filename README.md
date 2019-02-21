[![PyPI version](https://badge.fury.io/py/BondGraphTools.svg)](https://badge.fury.io/py/BondGraphTools)
[![Build Status](https://travis-ci.org/BondGraphTools/BondGraphTools.svg?branch=master)](https://travis-ci.org/BondGraphTools/BondGraphTools)
[![Test Coverage](https://api.codeclimate.com/v1/badges/4735c13a87b24d3a1899/test_coverage)](https://codeclimate.com/github/BondGraphTools/BondGraphTools/test_coverage)
# BondGraph - A Bond graph toolkit
## Summary

This toolkit is for rapid modelling and design of networked thermodynamic systems.
It is conceptually based upon the Bond Graph modelling methodology.

## Documentation

https://bondgraphtools.readthedocs.io/

## Installation

### Dependencies

BondGraph requires:
- python 3.6
- julia 0.6.4

Python dependencies:
- sympy>=1.1.1
- numpy>=1.14
- scipy>=1.0.1
- matplotlib>=2.2.2
- julia>=0.1.5
- diffeqpy>=0.4
- requests>=2.19

Julia dependencies:
 - PyCall
 - DifferentialEquations.jl

### Instructions:
1. Install python > 3.6 for your operating system.
2. Install Julia 0.6.4 (https://julialang.org/downloads/) for your operating
 system. _Julia 0.7 and 1.0 are not yet supported_
3. Make sure Julia 0.6.4 is in your os path. (test this by running `julia -v`)
4. Install using PyPI; `pip install BondGraphTools`