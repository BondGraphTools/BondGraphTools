import pytest
import sympy

import BondGraphTools as bgt
from BondGraphTools.algebra import extract_coefficients, smith_normal_form

def test_build_relations():
    c = bgt.new("C")

    eqns = c._build_relations()
    assert len(eqns) == 1

    eqn = eqns.pop()

    test_eqn = sympy.sympify("q_0 - c*e_0")

    assert eqn == test_eqn


def test_zero_junction_relations():
    r = bgt.new("R", value=sympy.symbols('r'))
    l = bgt.new("I", value=sympy.symbols('l'))
    c = bgt.new("C", value=sympy.symbols('c'))
    kvl = bgt.new("0", name="kvl")

    rlc = r + l + c + kvl

    rlc.connect(r, kvl)
    rlc.connect(l, kvl)
    rlc.connect(c, kvl)

    rels = kvl._build_relations()

    assert sympy.sympify("e_1 - e_0") in rels
    assert sympy.sympify("e_2 - e_0") in rels
    assert sympy.sympify("f_0 + f_1 + f_2") in rels


@pytest.mark.usefixture('rlc')
def test_basis_vectors(rlc):

    model_basis_vects = set()

    for component in rlc.components.values():
        for vects in component.basis_vectors:
            basis_vects = set(vects.values())

            assert not basis_vects & model_basis_vects

            model_basis_vects |= basis_vects


def test_build_junction_dict():
    c = bgt.new("C")
    kvl = bgt.new("0")

    bg = kvl+c
    bg.connect(kvl, c)
    index_map = {(c,"0"):0, (kvl,"0"):1}
    M = bg._build_junction_dict(index_map, offset=1)
    assert M[(0, 1)] == 1
    assert M[(0, 3)] == -1
    assert M[(1, 2)] == 1
    assert M[(1, 4)] == 1


def test_build_model_fixed_cap():
    c = bgt.new("C", value=0.001)

    eqns = c.constitutive_relations
    assert len(eqns) == 2

    test_eqn1 = sympy.sympify("q_0 - 0.001*e_0")
    test_eqn2 = sympy.sympify("dq_0-f_0")

    assert test_eqn1 in eqns
    assert test_eqn2 in eqns


@pytest.mark.usefixture("rlc")
def test_rlc_basis_vectors(rlc):

    tangent, ports, cv = rlc.basis_vectors

    assert len(tangent) == 2
    assert len(cv) == 0
    assert len(ports) == 6


def test_extract_coeffs_lin():
    eqn = sympy.sympify("y -2*x -3")
    local_map = {
        sympy.symbols("y"): 0,
        sympy.symbols("x"): 1
    }
    coords = [sympy.symbols("r_1"), sympy.symbols("r_0"), sympy.S(1)]

    lin, nlin = extract_coefficients(eqn, local_map, coords)
    assert lin[1] == -2
    assert lin[2] == -3
    assert lin[0] == 1
    assert not nlin


def test_extract_coeffs_nlin():
    eqn = sympy.sympify("y -2*x -3 + exp(x)")
    local_map = {
        sympy.symbols("y"): 0,
        sympy.symbols("x"): 1
    }
    coords = [sympy.symbols("r_1"), sympy.symbols("r_0"), sympy.S(1)]

    lin, nlin = extract_coefficients(eqn, local_map, coords)
    assert lin[1] == -2
    assert lin[2] == -3
    assert lin[0] == 1
    assert nlin == sympy.sympify("exp(r_0)")


def test_smith_normal_form():

    m = sympy.SparseMatrix(2,3,{(0,2):2, (1,1):1})
    mp = smith_normal_form(m)
    assert mp.shape == (3,3)
    assert mp[2, 2] != 0


def test_relations_iter():
    c = bgt.new("C", value=1)

    mappings = ({(c, 'q_0'): 0}, {(c,'0'): 0}, {})
    coords = list(sympy.sympify("dx_0,e_0,f_0,x_0, 1"))
    relations = c.get_relations_iterator(mappings, coords)

    d1 = {0:1, 2:-1} #dx/dt - f = 0
    d2 = {1:1, 3:-1} #e - x = 0$
    d1m = {0:-1, 2:1} #dx/dt - f = 0
    d2m = {1:-1, 3:1} #e - x = 0$

    for lin, nlin in relations:
        assert not nlin
        assert lin in (d1, d2, d1m, d2m)

