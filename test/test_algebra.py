import pytest
import sympy

import BondGraphTools as bgt
from BondGraphTools.algebra import extract_coefficients

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


def test_extract_coeffs():

    eqn = sympy.sympify("y -2*x -3 + exp(x)")
    local_map = {
        sympy.symbols("y"): 0,
        sympy.symbols("x"): 1
    }
    coords = [sympy.symbols("r_1"), sympy.symbols("r_0"), sympy.S(1)]

    lin, nlin = extract_coefficients(eqn, local_map, coords)

    print(lin)

    assert nlin == sympy.sympify("exp(r_0)")
