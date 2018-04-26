import pytest
import sympy

import BondGraphTools as bgt


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
        vects = component.basis

        basis_vects = set(vects.values())

        assert not basis_vects & model_basis_vects

        model_basis_vects |= basis_vects

