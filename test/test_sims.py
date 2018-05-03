import pytest
import numpy as np

import BondGraphTools as bgt

from BondGraphTools.base import ModelException
from BondGraphTools.sim_tools import simulate, _build_ode


def test_c_sim_fail():

    c = bgt.new("C")
    with pytest.raises(ModelException):

        t, x = simulate(c, timespan=[0, 1], initial_state=[1])

@pytest.mark.xfail
def test_c_se_build_ode():

    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = c + se + kcl + r

    bg.connect(c,kcl)
    bg.connect(r, kcl)
    bg.connect(se, kcl)

    # "dx_0 - u_0 + x_0"
    # so f(x,t) = exp(-t) - x

    func, jac = _build_ode(bg, input=['exp(-t)'])

    assert func(0, 0, 0, 0) == 1

@pytest.mark.skip
def test_c_se_sim():

    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = c + se + kcl + r

    bg.connect(c,kcl)
    bg.connect(r, kcl)
    bg.connect(se, kcl)

    with pytest.raises(ModelException) as ex:
        t, x = simulate(
            c, timespan=[0, 10], initial_state=[1]
        )
        assert "Control variable not specified" in ex.args
    t, x = simulate(
        bg, timespan=[0,10], initial_state=[1],
        input=lambda XT: [np.exp(-XT(-1))]
    )

    assert t[0] == 0
    assert t[-1] == 10
    t_cols, = t.shape
    assert t_cols > 1
    assert (1, t_cols) == x.shape
    assert x[0,0] == 1

    solution = (1+t)*np.exp(-t)
    solution = solution.reshape(1, t_cols)

    assert ((x - solution) < 0.001).all()

