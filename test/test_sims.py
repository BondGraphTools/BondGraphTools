import pytest
import numpy as np
from BondGraphTools import connect, new
from BondGraphTools.exceptions import ModelException
from BondGraphTools.sim_tools import simulate, _bondgraph_to_residuals

from .test_biochem import biochemical_cycle, open_cycle


def test_c_sim_fail():

    c = new("C")
    with pytest.raises(ModelException):
        t, x = simulate(c, timespan=[0, 1], x0=[1], dx0=[1])


def test_c_se_build_ode():

    c = new("C", value=1)
    se = new("Se")
    r = new("R", value=1)
    kcl = new("1")
    bg = new()
    bg.add([c, se, kcl, r])

    connect(c, (kcl, kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    # "dx_0 - u_0 + x_0"
    # so f(x,t) = exp(-t) - x

    def u(t, x, dx):
        return -np.exp(-t)

    residual_func, diff_vars = _bondgraph_to_residuals(bg, control_vars=[u])

    r = [0]
    residual_func(0, [0], [0], r)
    assert r == [-1]
    residual_func(0, [0], [2], r)
    assert r == [1]
    residual_func(0, [2], [0], r)
    assert r == [1]

    residual_func(0, [0], [1], r)
    assert r == [0]

    t_test = np.log(4)
    residual_func(t_test, [1 / 8], [1 / 8], r)
    assert r == [0]


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

    def u(t, x, dx):
        return -np.exp(-t)
    assert str(bg.constitutive_relations) == '[dx_0 + u_0 + x_0]'
    with pytest.raises(ModelException) as ex:
        _ = simulate(
            c, timespan=[0, 10], x0=[0]
        )
        assert "Control variable u_0 must be specified" in ex.args
    t, x = simulate(
        bg, timespan=[0, 10], x0=[0], control_vars=[u]
    )

    assert t[0] == 0
    assert t[-1] == 10

    assert (len(t), 1) == x.shape
    assert x[0, 0] == 0
    solution = t * np.exp(-t)

    res = abs(x - solution)
    assert (res < 0.001).all()


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

    def bang_bang(t, x, dx):
        return 1.5 if x >= 1 else -2.0

    t, x = simulate(
        bg, timespan=[0, 10], x0=[0], dx0=[1], control_vars=[bang_bang])

    assert (x[0, -1] - 1) < 0.001


@ pytest.mark.skip
@ pytest.mark.slow
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

    _ = simulate(
        bg, timespan=[0, 10], x0=[1, 0], control_vars=[1]
    )


@ pytest.mark.slow
def test_closed_cycle():
    model = biochemical_cycle()

    x0 = [2.0, 2.0, 2.0]
    tspan = (0, 1.0)
    K_X_vals = [1.0, 2.0, 3.0, 4.0]

    r1 = (model / "R:r1").params['r']['value']
    K_Y = (model / "C:Y").params['k']['value']
    def v1(r, K_X, K_Y, x_X, x_Y): return r * K_X * x_X - r * K_Y * x_Y

    for K_X in K_X_vals:
        (model / "C:X").set_param('k', K_X)
        t, x = simulate(model, tspan, x0, dt=0.01)
        # Check initial reaction rate
        assert v1(r1, K_X, K_Y, x[0][0], x[0][1]) == -4 + 2 * K_X
        # Check that the final reaction rate is small
        assert v1(r1, K_X, K_Y, x[-1][0], x[-1][1]) < 1e-4


@ pytest.mark.slow
def test_open_cycle():
    model = open_cycle()

    x0 = [2.0, 2.0, 2.0]
    tspan = (0, 10.0)
    K_A = 1
    x_A_vals = [0.02, 2.0, 200]
    expected_fluxes = [-3.0906200316489003, 0.0, 10.860335195530652]

    r1 = (model / "R:r1").params['r']['value']
    K_X = (model / "C:X").params['k']['value']
    K_Y = (model / "C:Y").params['k']['value']

    def v1(r, K_X, K_Y, K_A, x_X, x_Y, x_A): return r * \
        K_X * x_X * K_A * x_A - r * K_Y * x_Y

    for x_A, v_true in zip(x_A_vals, expected_fluxes):
        mu_A = np.log(K_A * x_A)
        (model / "SS:A").set_param('e', mu_A)
        t, x = simulate(model, tspan, x0)

        x_end = x[-1]
        flux = v1(r1, K_X, K_Y, K_A, x_end[0], x_end[1], x_A)
        assert flux == pytest.approx(v_true, abs=1e-6)
