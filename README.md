# BondGraph - A Bond graph toolkit
## Summary

This toolkit is for rapid modelling and design of networked thermodynamic systems.
It is conceptually based upon the Bond Graph modelling methodology.

Currently you can:
- Build power system models using basic linear components.
- Build networked models of biochemical systems
- Automatically derive constitutive relations.
- Run basic simulations
- Build new components

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
1. Install python 3.6 for your operating system.
2. Install Julia 0.6.4 (https://julialang.org/downloads/) for your operating
 system. _Julia 0.7 and 1.0 are not yet supported_
3. Make sure Julia 0.6.4 is in your os path. (test this by running `julia -v`)
4. checkout or download this repo
5. Set up Julia dependencies by executing the julia setup file `julia setup.jl`
6. Install the python 3.6 dependencies; `python3 -m pip install requirements.txt`


## Examples
