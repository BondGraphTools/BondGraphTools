import pytest
import pathlib
import sympy as sp

from BondGraphTools import load


file_path = pathlib.Path(__file__).parent / 'files'

def test_load_rlc():

    path = str(file_path / 'rlc.bg')

    model = load(path)
    uris = ["/C1",
            "/R1",
            "/L1",
            "/kcl",
            "/Sf"]

    c, r, l, kcl, sf =  (comp for uri in uris for comp in model.components if
                      comp.uri == uri)

    assert len(model.control_vars) == 1
    assert 'u_0' in model.control_vars

    assert c.params['C']['value'] == 10
    assert r.params['r']['value']  == 100
    assert l.params['L']['value']  == 10

    r_0, = r.ports
    c_0, = c.ports
    l_0, = l.ports
    sf_0,  = sf.ports
    kcl_0, kcl_1,kcl_2,kcl_3 = kcl.ports

    assert set(model.bonds) == {
        (r_0, kcl_0),
        (c_0, kcl_1),
        (l_0, kcl_2),
        (sf_0, kcl_3)
    }

    eqns = {
        sp.sympify("dx_0 - u_0 + x_0 / 1000 + x_1 / 10"),
        sp.sympify("dx_1 - x_0 / 10")
    }

    eqns_2 = {
        sp.sympify("dx_1 - u_0 + x_1 / 1000 + x_0 / 10"),
        sp.sympify("dx_0 - x_1 / 10")
    }

    assert (set(model.constitutive_relations) == eqns) or \
                (set(model.constitutive_relations) == eqns_2)


def test_load_rlc_parallel():
    path = str(file_path / 'rlc_parallel.bg')

    model = load(path)

    one, = (comp for comp in model.components if comp.metaclass == "1")

    assert one

    rel = model.constitutive_relations

    _, v = model.state_vars['x_0']

    if str(v) != 'p_0':
        eq1 = sp.sympify("dx_0 - x_1")
        eq2 = sp.sympify("dx_1 - u_0 + x_0 + x_1")
    else:
        eq1 = sp.sympify("dx_1 - x_0")
        eq2 = sp.sympify("dx_0 - u_0 + x_0 + x_1")

    for r in rel:
        assert r in (eq1, eq2)

def test_load_modular():

    model_1 = load(file_path / "modular.bg")

    assert model_1.uri == "/"
    assert model_1.name == "system"

    tree = set()
    def uri_tree(bg):
        tree.add(bg.uri)
        try:
            for c in bg.components:
                uri_tree(c)

        except AttributeError:
            pass
        return

    uri_tree(model_1)

    assert tree == {
        "/",
        "/Vs",
        "/Z",
        "/kvl",
        "/Vs/Sf",
        "/Vs/pout",
        "/Vs/kvl",
        "/Z/R1",
        "/Z/pin",
        "/Z/kvl",
        "/Z/L1",
        "/Z/C1"
    }




