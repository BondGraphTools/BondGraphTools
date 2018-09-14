import pytest
import pathlib
import logging

import sympy as sp

from BondGraphTools import new, connect, load, expose
from BondGraphTools.exceptions import *

file_path = pathlib.Path(__file__).parent / 'files'


def test_source_sensor():

    ss = new("SS", name="Plug 1")
    model = new()

    with pytest.raises(InvalidComponentException):
        expose(ss, "A")

    model.add(ss)

    assert len(model.ports) == 0
    expose(ss, "A")
    assert len(model.ports) == 1

def test_ss_exposure():
    model = new()
    ss = new("SS")
    sf = new("Sf", name="Sf")
    zero = new("0")

    model.add(ss,sf,zero)

    connect(sf, zero)
    connect(ss, zero)
    assert not model.ports
    assert set(model.control_vars.values()) == {
        (sf, 'f'), (ss, 'e'), (ss, 'f')
    }
    expose(ss, 'pin')
    assert model.ports
    p, = list(model.ports)
    assert p.name is "pin"
    assert set(model.control_vars.values()) == {(sf, 'f')}

def test_load_modular():
    model_1 = load(file_path / "modular.bg")

    assert model_1.name == "system"

    assert {c.name for c in model_1.components}  == {"Vs", "Z", "kvl"}
    for c in model_1.components:
        assert c
        assert c.parent is model_1

    Vs = model_1 / "Vs"
    Z = model_1 / "Z"
    kvl = model_1 / "kvl"

    assert Vs in model_1.components

    assert len(model_1.bonds) == 2
    assert len(Vs.ports) == 1
    assert list(Vs.ports)[0].name == 'A'
    assert Vs.state_vars == {}
    assert len(Vs.control_vars) == 1



def test_modularity():
    model_1 = load(file_path / "modular.bg")

    Vs = model_1 / "Vs"
    Z = model_1 / "Z"

    Vs_cr = Vs.constitutive_relations
    assert Vs_cr ==[sp.sympify('f_0 + u_0')]
    
    assert (set(model_1.constitutive_relations) == {
        sp.sympify("dx_0 - x_1"), sp.sympify("dx_1 - u_0 +x_0 + x_1")
    }) or ( set(model_1.constitutive_relations) == {
        sp.sympify("dx_1 - x_0"), sp.sympify("dx_0 - u_0 +x_0 + x_1")}
           )




