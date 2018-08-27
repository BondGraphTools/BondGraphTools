import pytest
import numpy as np
import sympy as sp
import BondGraphTools as bgt

from BondGraphTools.exceptions import ModelException
from BondGraphTools.sim_tools import simulate, _build_dae
from BondGraphTools.algebra import inverse_coord_maps

def test_c_sim_fail():

    c = bgt.new("C")
    with pytest.raises(ModelException):

        t, x = simulate(c, timespan=[0, 1], x0=[1],dx0=[1])

@pytest.mark.slow
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

    func, diff_vars = _build_dae(bg, control_vars=['-exp(-t)'])

    assert func(0, 0, 0, 0) == -1
    assert func(2, 0, 0, 0) == 1
    assert func(0, 2, 0, 0) == 1

@pytest.mark.slow
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
            c, timespan=[0, 10], x0=[0], dx0=[1]
        )
        assert "Control variable not specified" in ex.args
    t, x = simulate(
        bg, timespan=[0, 10], x0=[0], dx0=[1], control_vars=['-exp(-t)']
    )

    assert t[0] == 0
    assert t[-1] == 10
    t_cols, _ = t.shape
    assert t_cols > 1
    assert (t_cols, 1) == x.shape
    assert x[0, 0] == 0

    solution = (1+t)*np.exp(-t)
    solution = solution.reshape(t_cols, 1)

    assert ((x - solution) < 0.001).all()

@pytest.mark.slow
def test_c_se_sum_switch():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = c + se + kcl + r

    bg.connect(c, kcl)
    bg.connect(r, kcl)
    bg.connect(se, kcl)

    bang_bang = ["x_0 >= 1 ?  1.5: -2 "]

    t, x = simulate(
        bg, timespan=[0, 10], x0=[0],dx0=[1], control_vars=bang_bang
    )

    assert (x[0, -1] - 1) < 0.001


@pytest.mark.slow
def test_rlc():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    kvl = bgt.new("0")
    bg = c + se + kvl + r + l

    bg.connect(c, kvl)
    bg.connect(r, kvl)
    bg.connect(se, kvl)
    bg.connect(l, kvl)

    t, x = simulate(
        bg, timespan=[0, 10], x0=[0,0], control_vars=[1]
    )

