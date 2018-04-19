import BondGraphTools as bgt
from BondGraphTools.model import AtomicComponent, CompositeBondGraph


def test_new():
    c = bgt.new("C")

    assert isinstance(c, AtomicComponent)
    assert len(c.ports) == 1
    assert len(c.state_vars) == 1
    assert len(c.control_vars) == 1


def test_clone():

    c = bgt.new("C")
    cprime = bgt.new(c)

    assert id(c) == id(cprime)
    assert c.__dict__ == cprime.__dict__


def test_add():

    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    assert isinstance(bg, CompositeBondGraph)

    assert len(bg.ports) == 2
    assert len(bg.state_vars) == 1
    assert len(bg.control_vars) == 2

    assert bg is not c
    assert bg is not se

    assert c in bg
    assert se in bg

    assert c.parent is bg
    assert se.parent is bg
    assert not bg.parent


def test_set_param_float():

    c = bgt.new("C")

    capacitance = next(c.control_vars)
    assert c.params[capacitance] == capacitance

    val = 1.0
    c.params[capacitance] = val

    assert not c.contol_vars
    assert c.params[capacitance] == val

#
# def test_cv_merge():
#
#     c1 = bgt.new("C")
#     c2 = bgt.new("C")
#
#     se = bgt.new("Se")
#
#




