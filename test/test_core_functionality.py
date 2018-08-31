import pytest

import BondGraphTools as bgt
from BondGraphTools.exceptions import InvalidComponentException, InvalidPortException
from BondGraphTools.atomic import BaseComponent


def test_new():
    c = bgt.new("C")

    assert isinstance(c, BaseComponent)
    assert len(c.ports) == 1
    assert len(c.state_vars) == 1
    assert len(c.params) == 1


def test_new_zero():
    j = bgt.new("0")

    assert isinstance(j, BaseComponent)
    assert len(j.ports) == 0


def test_new_parameterised():

    c = bgt.new("C", value=1.0)
    v = next(iter(c.params.values()))
    assert v == 1.0


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


def test_equal():

    c_1 = bgt.new("C", value=1)
    c_2 = bgt.new("C", value=1)

    assert c_1 == c_2
    assert c_1 is c_1
    assert c_1 is not c_2

@pytest.mark.skip
def test_iadd():
    c_1 = bgt.new("C", value=1)
    c_2 = bgt.new("C", value=1)
    c_3 = bgt.new("C", value=1)

    bg = bgt.new(name="bg")
    bg += c_1 + c_2 + c_3

    assert c_1 in bg
    assert c_2 in bg
    assert c_3 in bg


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

    comp, port_id = bg._find_port(c)

    assert comp is c
    assert port_id in c.ports


def test_find_port_by_name():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    c_str = list(bg.components.keys())[0]

    comp, port_id = bg._find_port(c_str)

    assert comp is c
    assert port_id in c.ports


def test_port_connection():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    assert len(bg.ports) ==2
    bg.connect(c, se)
    assert len(bg.ports) == 0


def test_connect_components():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    bg.connect(c, se)
    (c1, p1), (c2, p2) = bg.bonds[0]
    assert isinstance(c1, BaseComponent)
    assert isinstance(c2, BaseComponent)
    assert c1 in (c, se)
    assert c2 in (c, se)


def test_connect():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    k1, k2 = tuple(bg.components)

    bg.connect(k1, k2)

    (c1, p1), (c2, p2) = bg.bonds[0]
    assert isinstance(c1, BaseComponent)
    assert isinstance(c2, BaseComponent)
    assert c1 in (c, se)
    assert c2 in (c, se)


def test_connect_ports():
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    k1 = bg.ports[0]
    k2 = bg.ports[1]

    bg.connect(k1, k2)

    (c1, p1), (c2, p2) = bg.bonds[0]
    assert isinstance(c1, BaseComponent)
    assert isinstance(c2, BaseComponent)
    assert c1 in (c, se)
    assert c2 in (c, se)


def test_disconnect_ports():

    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se

    bg.connect(c, se)

    with pytest.raises(InvalidPortException):
        bg.connect(c, se)

    _ = bg.bonds.pop()
    assert not bg.bonds

    bg.connect(c, se)
    (c1, p1), (c2, p2) = bg.bonds[0]
    assert c1 in (c, se)
    assert c2 in (c, se)

    bg.disconnect(c, se)
    assert not bg.bonds

def test_many_port():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se + j

    bg.connect(c, j)
    bg.connect(se, j)

    for (c1, p1), (c2, p2) in bg.bonds:
        assert c1 in (c, se)
        assert c2 is j

    assert len(j.ports) == 2


def test_delete_one_from_many_port():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se + j

    bg.connect(c, j)
    bg.connect(se, j)
    (p,q), (r,s) = bg.bonds[1]
    assert p in (se, j)
    assert r in (se, j)

    bg.disconnect(c, j)
    assert len(bg.bonds) == 1
    (pp, qq), (rr, ss) = bg.bonds[0]

    assert pp in (se, j)

    assert pp is p
    assert rr is r
    assert qq is q
    assert ss is s


def test_set_param():
    c = bgt.new("C")

    c.params["c"] = 1

    assert c.params["c"] == 1


def test_set_compound_param():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se + j

    assert len(bg.params) == 2

    bg.params[0] = 1


def test_connect_unadded():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + j

    bg.connect(c, j)

    assert j in bg
    assert c in bg

    assert j is not se
    with pytest.raises(InvalidComponentException):
        bg.connect(se, j)


def test_disconnect_multiport():
    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")
    one = bgt.new("1")

    bg = zero+r+c+one

    bg.connect((zero,0), (one,0))
    bg.connect(r,zero)
    bg.connect(c, one)

    assert len(bg.bonds) == 3
    bg.disconnect(zero, one)
    assert len(bg.bonds) == 2
    assert ((zero, 0), (one,0)) not in bg.bonds


def test_disconnect_component():
    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")
    one = bgt.new("1")

    bg = zero+r+c+one

    bg.connect((zero,0), (one,0))
    bg.connect(r,zero)
    bg.connect(c, one)

    assert len(bg.bonds) == 3
    bg.disconnect(c)
    assert len(bg.bonds) == 2


def test_swap_component_1():
    ###
    # We must be able to swap a 1-port for a 1-port
    #
    #

    zero = bgt.new("0")
    r = bgt.new("R", value=1)
    c = bgt.new("C", value=1)

    bg = zero+r+c

    bg.connect(r,zero)
    bg.connect(c, zero)
    assert c in bg.components

    assert bg.bonds == [
        ((r, 0), (zero, 0)),
        ((c, 0), (zero, 1))
    ]
    assert len(bg.state_vars) == 1
    Sf = bgt.new('Sf')
    bg.replace(c, Sf)
    assert len(bg.state_vars) == 0
    assert len(bg.control_vars) == 1

    assert bg.bonds == [
        ((r,0), (zero, 0)),
        ((Sf, 0), (zero, 1))
    ]

    assert c not in bg.components
    assert Sf in bg.components

    bg.replace(Sf, c)

    assert bg.bonds == [
        ((r, 0), (zero, 0)),
        ((c, 0), (zero, 1))
    ]

    assert c in bg.components
    assert Sf not in bg.components


def test_swap_components_2():
    ###
    # We should also be able to swap zero and one juncitons.
    # However, due to the nature of one ports, we assume that simga = 1.

    zero = bgt.new("0")
    r = bgt.new("R")
    c = bgt.new("C")

    bg = zero+r+c

    bg.connect(r,zero)
    bg.connect(c, zero)

    assert bg.bonds == [
        ((r, 0), (zero, 0)),
        ((c, 0), (zero, 1))
    ]
    one = bgt.new('1')
    bg.replace(zero, one)

    assert bg.bonds == [
        ((r, 0), (one, 0)),
        ((c, 0), (one, 1))
    ]

    assert one in bg.components
    assert zero not in bg.components