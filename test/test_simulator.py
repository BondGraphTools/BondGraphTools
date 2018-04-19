import numpy

import BondGraphTools as bgt
from BondGraphTools.sim_tools import simulate


def test_c_sim():

    c = bgt.new("C", value=0.001)
    r = bgt.new("R", value=1000)
    se = bgt.new("Se")
    kcl = bgt.new("0")

    bg = c + r + se + kcl

    bg.connect([(c, kcl),
                (r, kcl),
                (se, kcl)])

    observables = se.ports[0] + c.state_vars

    t, X = simulate(
        bg,
        intial=[0],
        params=["exp(-t)"],
        observable=observables
    )

    t_test = numpy.exp(-t)

    assert abs(t - X[0]) < 0.0001