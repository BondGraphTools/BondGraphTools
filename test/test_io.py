import pytest
import pathlib
import sympy as sp
import os
import yaml
from BondGraphTools import load, save, new, connect
import BondGraphTools.datamodel as dm
import logging

file_path = pathlib.Path(__file__).parent / 'files'
logging.basicConfig(level=logging.DEBUG)

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


def test_load_model_from_modular():

    Vs = load(file_path / "modular.bg", model='source', as_name="Source 1")
    Vs2 = load(file_path / "modular.bg", model='source', as_name="Source 2")
    assert Vs is not Vs2

    assert len(Vs.components) == 3
    assert len(Vs2.components) == 3

@pytest.mark.usefixture("rlc")
def test_save_build_component(rlc):

    r, = (c for c in rlc.components if c.metaclass == "R")
    c, = (c for c in rlc.components if c.metaclass == "C")
    r.params["r"] = 10
    c.params["C"] = None

    r_str = dm._build_component_string(r)
    c_str = dm._build_component_string(c)
    assert r_str == f"{r.name} base/R r=10"
    assert c_str == f"{c.name} base/C"

@pytest.mark.usefixture("rlc")
def test_save_build_model(rlc):
    r, = (c for c in rlc.components if c.metaclass == "R")
    c, = (c for c in rlc.components if c.metaclass == "C")
    l, = (c for c in rlc.components if c.metaclass == "I")
    kvl, =  (c for c in rlc.components if c.metaclass == "0")
    model_dict = dm._build_model_data(rlc, {})

    test_strings = {f"{r.name} base/R r=1",
                    f"{c.name} base/C C=1",
                    f"{l.name} base/I L=1",
                    f"{kvl.name} base/0"}

    assert set(model_dict["components"]) == test_strings

    assert set(model_dict["netlist"]) == {
        f"{r.name} {kvl.name}", f"{c.name} {kvl.name}", f"{l.name} {kvl.name}"
    }

def test_build_templated_model():
    root = new(name="System")
    model = new(name="Vs")
    root.add(model)
    ss = new("SS", name="pout")
    Se = new("Se", name="Vs")
    zero = new("0", name="0_0")
    model.add(ss,Se,zero)
    connect(ss,zero)
    connect(Se, zero)

    file = file_path / "temp.bg"
    assert root.name == "System"

    with TempFile(file):

        save(root, file)

        with open(file, 'r') as fs:
            temp_data = yaml.load(fs)

    assert temp_data["root"] == root.name
    assert temp_data["models"].keys() == {
        "/", "/Vs"
    }

    assert temp_data["models"]["/"]["components"] == ["Vs /Vs"]
    assert set(temp_data["models"]["/Vs"]["components"]) == {
        "pout base/SS", "Vs base/Se", "0_0 base/0"
    }
    assert set(temp_data["models"]["/Vs"]["netlist"]) == {
        "pout 0_0", "Vs 0_0"
    }


@pytest.mark.usefixture("rlc")
def test_rlc_save(rlc):
    filename = str(file_path / "test_rlc.bg")
    rlc.name = "RLC"
    r, = (c for c in rlc.components if c.metaclass == "R")
    c, = (c for c in rlc.components if c.metaclass == "C")
    l, = (c for c in rlc.components if c.metaclass == "I")
    kvl, =  (c for c in rlc.components if c.metaclass == "0")

    r.name= "R1"
    c.name= "C1"
    l.name= "L1"

    with TempFile(filename):
        save(rlc, filename)

        with open(filename, 'r') as fs:
            test_data = yaml.load(fs)

    assert test_data["version"] == dm.FILE_VERSION
    assert test_data["root"] == "RLC"
    assert set(test_data["models"].keys()) == {"/"}
    model_data = test_data["models"]["/"]

    assert set(model_data.keys()) == {"components", "netlist"}

    assert set(model_data["components"]) == {
        "R1 base/R r=1", "C1 base/C C=1", "L1 base/I L=1","kvl base/0"}

    assert set(model_data["netlist"]) == {
        "R1 kvl", "C1 kvl", "L1 kvl"
    }


def test_save_idempotent():
    model_1 = load(file_path / "modular.bg")
    new_file = str(file_path / "modular_2.bg")

    with TempFile(new_file):
        save(model_1, new_file)
        model_2 = load(new_file)

    assert model_1 is not model_2
    assert model_1.name == model_2.name
    assert_components_are_equal(model_1, model_2)

def assert_components_are_equal(bg_1, bg_2):
    names = {c.name for c in bg_1.components} & {c.name for c in bg_2.components}
    if len(names) != len(bg_1.components):
        assert False, names

    for name in names:
        c1, = (c for c in bg_1.components if c.name == name)
        c2, = (c for c in bg_2.components if c.name == name)
        if hasattr(c1, 'components'):
            assert_components_are_equal(c1, c2)
        else:
            assert_atomics_are_equal(c1,c2)

def assert_atomics_are_equal(c1,c2):
    ports_1 = {repr(p) for p in c1.ports}
    ports_2 = {repr(p) for p in c2.ports}
    assert ports_1 == ports_2

    for p in ["__library__", "__component__"]:
        assert c1.__dict__[p] == c2.__dict__[p]

    assert (hasattr(c1, 'params') == hasattr(c2, 'params'))
    if hasattr(c1, 'params'):
        assert c1.params == c2.params


class TempFile():
    def __init__(self, filename):
        self.file = filename

    def __enter__(self):
        try:
            os.remove(self.file)
        except FileNotFoundError:
            pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        with open(self.file, 'r' ) as file:
            print("".join(line for line in file.readlines()))
        try:
            os.remove(self.file)
        except FileNotFoundError:
            pass