import pytest

import BondGraphTools as bgt
from BondGraphTools.base import InvalidComponentException
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

    k1, k2 = tuple(bg.ports)

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

    with pytest.raises(InvalidComponentException):
        bg.connect(c, se)

    _ = bg.bonds.pop()
    assert not bg.bonds

    bg.connect(c, se)
    (c1, p1), (c2, p2) = bg.bonds[0]
    assert c1 in (c, se)
    assert c2 in (c, se)


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

def test_set_compuond_param():

    j = bgt.new("0")
    c = bgt.new("C")
    se = bgt.new("Se")

    bg = c + se + j

    assert len(bg.params) == 2
