import pytest

import BondGraphTools as bgt
from BondGraphTools import *
from BondGraphTools.exceptions import InvalidComponentException, InvalidPortException
from BondGraphTools.atomic import *


def test_new():
    c = bgt.new("C")

    assert isinstance(c, BaseComponent)
    assert len(c.ports) == 1
    assert len(c.state_vars) == 1
    assert len(c.params) == 1


def test_new_zero():
    j = bgt.new("0")

    assert isinstance(j, EqualEffort)
    assert len(j.ports) == 0


def test_new_parameterised():

    c = bgt.new("C", value=1.0)
    v = next(iter(c.params.values()))
    assert v == 1.0


def test_add():

    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se])


    assert len(bg.internal_ports) == 2
    assert len(bg.state_vars) == 1

    assert bg is not c
    assert bg is not se

    assert c in bg.components
    assert se in bg.components


def test_equal():

    c_1 = bgt.new("C", value=1)
    c_2 = bgt.new("C", value=1)

    assert c_1 == c_2
    assert c_1 is c_1
    assert c_1 is not c_2

def test_connect_components():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se])

    connect(c, se)
    bond = list(bg.bonds)[0]
    assert isinstance(bond.tail.component, BaseComponent)
    assert isinstance(bond.head.component, BaseComponent)
    assert bond.head.component in (c, se)
    assert bond.tail.component in (c, se)


def test_connect():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se])

    k1, k2 = tuple(bg.components)

    connect(k1, k2)

    (c1, p1), (c2, p2) = list(bg.bonds)[0]
    assert isinstance(c1, BaseComponent)
    assert isinstance(c2, BaseComponent)
    assert c1 in (c, se)
    assert c2 in (c, se)

def test_one_port():
    c = new('C')
    se = new('Se')
    r = new('R')
    one = new('1')
    bg = new()
    bg.add(c,se,r,one)

    connect(c, one.non_inverting)
    connect(se, one)
    connect(one, r)

    p0,p1,p2 = tuple(one.ports)

    assert p0.weight == 1
    assert p1.weight == 1
    assert p2.weight == -1


def test_connect_ports():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se])

    k1 = bg.internal_ports[0]
    k2 = bg.internal_ports[1]


    connect(k1, k2)

    (c1, p1), (c2, p2) = list(bg.bonds)[0]
    assert isinstance(c1, BaseComponent)
    assert isinstance(c2, BaseComponent)
    assert c1 in (c, se)
    assert c2 in (c, se)


def test_disconnect_ports():

    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se])

    connect(c, se)

    with pytest.raises(InvalidPortException):
        connect(c, se)

    (c1, p1), (c2, p2) = list(bg.bonds)[0]
    assert c1 in (c, se)
    assert c2 in (c, se)

    disconnect(c, se)
    assert not bg.bonds

def test_many_port():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se, j])

    connect(c, j)
    connect(se, j)

    for (c1, p1), (c2, p2) in bg.bonds:
        assert c1 in (c, se)
        assert c2 is j

    assert len(j.ports) == 2


def test_delete_one_from_many_port():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se, j])

    connect(c, j)
    connect(se, j)

    (p,q), (r,s) = bg.bonds[1]

    assert p in (se, j)
    assert r in (se, j)
    assert len(bg.bonds) == 2
    disconnect(c, j)
    assert len(bg.bonds) == 1
    (pp, qq), (rr, ss) = bg.bonds[0]
    assert pp in (se, j)
    assert rr in (se, j)


def test_set_param():
    c = bgt.new("C")

    c.params["c"] = 1

    assert c.params["c"] == 1


def test_set_compound_param():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, se, j])


    assert len(bg.params) == 2

    bg.params[0] = 1


def test_connect_unadded():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = bgt.new()
    bg.add([c, j])

    connect(c, j)

    assert j in bg.components
    assert c in bg.components
    assert j is not se

    with pytest.raises(InvalidComponentException):
        connect(se, j)


def test_disconnect_multiport():
    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")
    one = bgt.new("1")
    bg = bgt.new()
    bg.add([zero, r, c, one])

    connect(zero, one.non_inverting)
    connect(r,zero)
    connect(c, one.inverting)

    assert len(bg.bonds) == 3
    disconnect(zero, one)
    assert len(bg.bonds) == 2
    assert ((zero, 0), (one,0)) not in bg.bonds



def test_disconnect_component():
    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")
    one = bgt.new("1")

    bg = bgt.new()
    bg.add([zero, r, c, one])

    connect(c, (one, one.non_inverting))

    with pytest.raises(TypeError):
        disconnect(c)


def test_remove_component():
    zero = bgt.new("0")
    r = bgt.new("R", value=1)
    c = bgt.new("C", value=1)

    bg = bgt.new()
    bg.add([zero, r, c])


    assert c in bg.components
    bg.remove(c)
    assert c not in bg.components

    bg.add(c)

    connect(r, zero)
    connect(c, zero)
    r_p, = r.ports
    c_p, = c.ports
    z0,z1, = zero.ports
    assert c in bg.components
    assert set(bg.bonds) == {
        (r_p, z0),
        (c_p, z1)
    }

def test_remove_components_fail_states():
    zero = bgt.new("0")
    r = bgt.new("R", value=1)
    c = bgt.new("C", value=1)

    bg = bgt.new()
    bg.add([zero, r])

    with pytest.raises(InvalidComponentException):
        bg.remove(c)


def test_swap_component_1():
    ###
    # We must be able to swap a 1-port for a 1-port
    #
    #

    zero = bgt.new("0")
    r = bgt.new("R", value=1)
    c = bgt.new("C", value=1)

    bg = bgt.new()
    bg.add([zero, r, c ])

    connect(r,zero)
    connect(c, zero)
    assert c in bg.components
    r_p, = r.ports
    c_p, = c.ports
    z0,z1, = zero.ports

    assert set(bg.bonds) == {
        (r_p, z0),
        (c_p, z1)
    }

    assert len(bg.state_vars) == 1
    Sf = bgt.new('Sf')
    swap(c, Sf)
    sf_port, = Sf.ports
    assert len(bg.state_vars) == 0
    assert len(bg.control_vars) == 1

    assert set(bg.bonds) == {
        (r_p, z0),
        (sf_port, z1)
    }

    assert c not in bg.components
    assert Sf in bg.components

    swap(Sf, c)

    assert set(bg.bonds) == {
        (r_p, z0),
        (c_p, z1)
    }

    assert c in bg.components
    assert Sf not in bg.components


def test_swap_components_2():

    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")

    bg = bgt.new()
    bg.add([zero, r, c ])

    connect(r,zero)
    connect(c, zero)

    p1, = r.ports
    p2, = c.ports
    p3, p4, = zero.ports
    assert set(bg.bonds) == {
        (p1, p3),
        (p2, p4)
    }
    one = bgt.new('1')
    with pytest.raises(InvalidComponentException):
        swap(zero, one)

def test_swap_failstate():
    zero = bgt.new("0")
    r = bgt.new("R", value=1)
    c = bgt.new("C", value=1)
    l = bgt.new("I")

    bg = bgt.new()
    bg.add([zero, r, c ])
    connect(r, zero)
    connect(c, zero)

    with pytest.raises(InvalidComponentException):
        swap(zero, l)
