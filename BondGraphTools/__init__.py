"""
BondGraphTools
==============
BondGraphTools is a python library for symbolic modelling
of energetic systems.

Package Documentation::

    https://bondgraphtools.readthedocs.io/

Source::

    https://github.com/BondGraphTools/BondGraphTools

Bug reports:

    https://github.com/BondGraphTools/BondGraphTools/issues


Simple Example
--------------

Build and simulate a RLC driven RLC circuit::

    import BondGraphTools as bgt

    # Create a new model
    model = bgt.new(name="RLC")

    # Create components
    # 1 Ohm Resistor
    resistor = bgt.new("R", name="R1", value=1.0)

    # 1 Henry Inductor
    inductor = bgt.new("L", name="L1", value=1.0)
    # 1 Farad Capacitor
    capacitor = bgt.new("C", name="C1", value=1.0)
    # Conservation Law
    law = bgt.new("0") # Common voltage conservation law

    # Connect the components
    connect(law, resistor)
    connect(law, capacitor)
    connect(law, inductor)

    # produce timeseries data
    t, x = simulate(model, x0=[1,1], timespan=[0, 10])

Bugs
----

Please report any bugs `here <https://github.com/BondGraphTools/BondGraphTools/issues>`_,
or fork the repository and submit a pull request.

License
-------
Released under the Apache 2.0 License::

    Copyright (C) 2018
    Peter Cudmore <peter.cudmore@uqconnect.edu.au>
"""

from .actions import *
from .fileio import save, load
from .compound import BondGraph
from .sim_tools import simulate
from .view import draw
from .version import __version__ as version
