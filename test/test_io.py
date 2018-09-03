import pytest
import pathlib
import sympy as sp

from BondGraphTools import load


file_path = pathlib.Path(__file__).parent / 'files'

def test_load_rlc():

    path = str(file_path / 'rlc.bg')

    model = load(path)

    c = model.components['C1']
    r = model.components['R1']
    l = model.components['L1']
    kcl = model.components['kcl']
    sf = model.components['Sf']

    assert len(model.control_vars) == 1
    assert 'u_0' in model.control_vars

    assert c.params['C']['value'] == 10
    assert r.params['r']['value']  == 100
    assert l.params['L']['value']  == 10

    assert model.bonds == [
        ((r,0), (kcl, 0)),
        ((c, 0), (kcl, 1)),
        ((l, 0), (kcl, 2)),
        ((sf, 0), (kcl, 3))
    ]

    eqns = {
        sp.sympify("dx_0 - u_0 + x_0 / 1000 + x_1 / 10"),
        sp.sympify("dx_1 - x_0 / 10")
    }

    assert set(model.constitutive_relations) == eqns