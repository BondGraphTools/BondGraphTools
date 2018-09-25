import pytest
import numpy as np
import sympy as sp
from BondGraphTools import *
from BondGraphTools.config import config
from BondGraphTools.exceptions import ModelException
from BondGraphTools.sim_tools import simulate, to_julia_function_string
from BondGraphTools.algebra import inverse_coord_maps

@pytest.mark.slow
def test_c_sim_fail():

    c = new("C")
    with pytest.raises(ModelException):
        t, x = simulate(c, timespan=[0, 1], x0=[1],dx0=[1])

@pytest.mark.slow
def test_c_se_build_ode():

    c = new("C", value=1)
    se = new("Se")
    r = new("R", value=1)
    kcl = new("1")
    bg = new()
    bg.add([c, se, kcl, r])

    connect(c,(kcl, kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    # "dx_0 - u_0 + x_0"
    # so f(x,t) = exp(-t) - x

    j = config.julia
    func_str, diff_vars = to_julia_function_string(bg,
                                                   control_vars=['-exp(-t)'],
                                                   in_place=False)
    func = j.eval(func_str)
    assert func(0, 0, 0, 0) == -1
    assert func(2, 0, 0, 0) == 1
    assert func(0, 2, 0, 0) == 1

@pytest.mark.slow
def test_c_se_sim():

    c = new("C", value=1)
    se = new("Se")
    r = new("R", value=1)
    kcl = new("1")
    bg = new()
    bg.add([c, se, kcl, r])

    connect(c, (kcl, kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

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
    c = new("C", value=1)
    se = new("Se")
    r = new("R", value=1)
    kcl = new("1")
    bg = new()
    bg.add([c, se, kcl, r])

    connect(c, (kcl, kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    bang_bang = ["x_0 >= 1 ?  1.5: -2 "]

    t, x = simulate(
        bg, timespan=[0, 10], x0=[0],dx0=[1], control_vars=bang_bang
    )

    assert (x[0, -1] - 1) < 0.001

@pytest.mark.skip
@pytest.mark.slow
def test_rlc():
    c = new("C", value=1)
    se = new("Se")
    r = new("R", value=1)
    l = new("I", value=1)
    kvl = new("0")

    bg = new()
    bg.add([c, se, kvl, r, l])


    connect(c, kvl)
    connect(r, kvl)
    connect(se, kvl)
    connect(l, kvl)

    t, x = simulate(
        bg, timespan=[0, 10], x0=[1,0], control_vars=[1]
    )

