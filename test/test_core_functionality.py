import BondGraphTools as bgt
from BondGraphTools.model import AtomicComponent, BondGraph


def test_new():
    c = bgt.new("C")

    assert isinstance(c, AtomicComponent)
    assert len(c.ports) == 1
    assert len(c.state_vars) == 1
    assert len(c.params) == 1


def test_new_parameterised():

    c = bgt.new("C", value=1.0)
    v = next(iter(c.params.values()))
    assert v["value"] == 1.0

def test_clone():

    c = bgt.new("C")
    cprime = bgt.new(c)

    assert id(c) != id(cprime)
    assert c == cprime

def test_add():

    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    assert len(bg.ports) == 2
    assert len(bg.state_vars) == 1

    assert bg is not c
    assert bg is not se

    assert c in bg
    assert se in bg

    # assert c.parent is bg
    # assert se.parent is bg
    # assert not bg.parent


def test_equal():

    c_1 = bgt.new("C", value=1)
    c_2 = bgt.new("C", value=1)

    assert c_1 == c_2
    assert c_1 is not c_2


def test_long_sum():
    c_1 = bgt.new("C", value=1)
    c_2 = bgt.new("C", value=1)
    c_3 = bgt.new("C", value=1)

    bg = c_1 + c_2 + c_3

    assert c_1 in bg
    assert c_2 in bg
    assert c_3 in bg


def test_find_port():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    comp_id, port_id = bg._find_port(c)

    assert bg[comp_id] is c
    assert port_id in c.ports


def test_connect_components():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    bg.connect(c, se)


def test_connect():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    k1, k2 = tuple(bg.components)
    print(k1)
    print(k2)

    bg.connect(k1, k2)


def test_connect_ports():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    k1, k2 = tuple(bg.ports)

    bg.connect(k1, k2)




# def test_set_param_float():
#
#     c = bgt.new("C")
#
#     capacitance = next(c.param)
#     assert c.params[capacitance] == capacitance
#
#     val = 1.0
#     c.params[capacitance] = val
#
#     assert not c.contol_vars
#     assert c.params[capacitance] == val

#
# def test_cv_merge():
#
#     c1 = bgt.new("C")
#     c2 = bgt.new("C")
#
#     se = bgt.new("Se")
#
#




