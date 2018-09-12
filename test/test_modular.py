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


def test_build():
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

    model_2 = new()

    Se = new("Sf", value=1)
    kvl = new("0")
    r = new('R', value=1)
    l = new('I', value=1)
    c = new("C", value=1)

    model_2.add(Se, r,l, c, kvl)
    for comp in model_2.components:
        if comp is not kvl:
            connect(comp, kvl)

    Vs = model_1 / "Vs"
    Z = model_1 / "Z"

    Vs_cr = Vs.constitutive_relations
    assert Vs_cr ==[sp.sympify('e_0 - u_0')]


    rel_1 = set(model_1.constitutive_relations)
    rel_2 = set(model_2.constitutive_relations)

    assert rel_1 == rel_2
