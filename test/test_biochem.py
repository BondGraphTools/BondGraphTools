import pytest

import sympy
import numpy as np

import BondGraphTools as bgt
from BondGraphTools import connect
from BondGraphTools.exceptions import InvalidPortException
from BondGraphTools.algebra import extract_coefficients, inverse_coord_maps, get_relations_iterator
from BondGraphTools.reaction_builder import Reaction_Network

import logging


def test_make_a_to_b():

    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])

    Re = bgt.new("Re", library="BioChem", value={"R": 1, "T": 1})
    a_to_b = bgt.new()
    a_to_b.add(A, Re, B)

    connect(A, (Re, 0))
    connect(B, (Re, 1))

    with pytest.raises(InvalidPortException):
        connect(A, Re)

    assert Re.control_vars == ['r']
    assert list(a_to_b.state_vars.keys()) == ['x_0', 'x_1']
    assert list(a_to_b.control_vars.keys()) == ['u_0']
    assert not a_to_b.ports

    state_basis, _, control_basis = a_to_b.basis_vectors
    (x_A, dx_A), = (k for k, (v, _) in state_basis.items() if v is A)
    (x_B, dx_B), = (k for k, (v, _) in state_basis.items() if v is B)
    (r,) = (k for k, (v, _) in control_basis.items() if v is Re)

    solutions = {dx_A + r * x_A - r * x_B, dx_B - r * x_A + r * x_B}
    relations = set(a_to_b.constitutive_relations)
    assert not solutions ^ relations


def test_make_a_to_b_inplicit():

    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])

    Re = bgt.new("Re", library="BioChem", value={"R": 1, "T": 1})
    a_to_b = bgt.new()
    a_to_b.add(A, Re, B)

    connect(A, Re)
    connect(B, Re)

    assert Re.control_vars == ['r']
    assert list(a_to_b.state_vars.keys()) == ['x_0', 'x_1']
    assert list(a_to_b.control_vars.keys()) == ['u_0']
    assert not a_to_b.ports


def test_re_con_rel():
    Re = bgt.new("Re", library="BioChem", value={"R": 1, "T": 1})

    coords = list(sympy.sympify("e_0,f_0,e_1,f_1,r"))

    mappings, coords = inverse_coord_maps(*Re.basis_vectors)
    assert mappings
    assert coords

    for r in get_relations_iterator(Re, mappings, coords):
        assert r in [
            ({1: 1, 3: 1}, 0), ({1: 1}, sympy.sympify("-r*exp(e_0) + r*exp(e_1)"))
        ]


def test_a_to_b_model():
    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])

    Re = bgt.new("Re", library="BioChem", value={'r': 1, "R": 1, "T": 1})

    Y_A = bgt.new('1')
    Y_B = bgt.new('1')

    a_to_b = bgt.new()
    a_to_b.add(A, Re, B, Y_A, Y_B)

    connect(A, Y_A.non_inverting)
    connect(B, Y_B.non_inverting)
    connect((Re, 0), Y_A.inverting)
    connect((Re, 1), Y_B.inverting)

    eqns = {
        sympy.sympify("dx_0 + x_0 -x_1"), sympy.sympify("dx_1 + x_1 -x_0")
    }
    for relation in a_to_b.constitutive_relations:
        assert relation in eqns


def test_ab_to_c_model():

    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    C = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    Re = bgt.new("Re", library="BioChem", value={'r': 1, "R": 1, "T": 1})
    Y_AB = bgt.new('1')
    Y_C = bgt.new('1')

    bg = bgt.new()
    bg.add(A, B, Re, Y_AB, C, Y_C)
    connect(A, Y_AB)
    connect(B, Y_AB)
    connect(Y_AB, Re)
    connect(Re, Y_C)
    connect(Y_C, C)

    state_basis, _, _ = bg.basis_vectors
    (x_0, dx_0), = (k for k, (v, _) in state_basis.items() if v is A)
    (x_1, dx_1), = (k for k, (v, _) in state_basis.items() if v is B)
    (x_2, dx_2), = (k for k, (v, _) in state_basis.items() if v is C)

    eqns = {dx_0 + x_0 * x_1 - x_2,
            dx_1 + x_0 * x_1 - x_2,
            dx_2 - x_0 * x_1 + x_2}

    relations = set(bg.constitutive_relations)
    assert not relations ^ eqns


def test_new_reaction_network():

    rn = Reaction_Network(name="A+B to C")
    rn.add_reaction(
        "A+B=C"
    )
    assert rn.species == ["A", "B", "C"]
    assert rn.forward_stoichiometry == sympy.Matrix([[1], [1], [0]])
    assert rn.reverse_stoichiometry == sympy.Matrix([[0], [0], [1]])
    assert rn.stoichiometry == sympy.Matrix([[-1], [-1], [1]])


def test_rn_to_bond_graph():
    rn = Reaction_Network(name="A+B to C", reactions="A+B=C")

    system = rn.as_network_model(normalised=True)
    assert len(system.state_vars) == 3
    assert len(system.control_vars) == 4


def test_cat_rn():
    reactions = [
        "E + S = ES", "ES = E+P"
    ]
    rn = Reaction_Network(reactions=reactions, name="Catalysed Reaction")
    assert rn._species["ES"] == 2
    assert len(rn._reactions) == 2


# TODO: fix when we rework parameters
def test_nlin_se():
    rn = Reaction_Network(name="A+B to C", reactions="A+B=C")
    system = rn.as_network_model(normalised=True)

    for param in system.params:
        system.set_param(param, 1)

    Ce_A = system / "A"
    Ce_B = system / "B"
    Ce_C = system / "C"

    assert Ce_A is not Ce_B

    Y = system / "AB"

    bgt.disconnect(Ce_A, Y)

    J_A = bgt.new("0", name="Ce_A")
    Se_A = bgt.new('Se', value=1, name='e=1')

    system.add(J_A),
    system.add(Se_A)
    connect(J_A, Y)
    connect(Ce_A, J_A)
    connect(Se_A, J_A)

    assert len(system.state_vars) == 3
    state_basis, _, _ = system.basis_vectors

    (x0, dx0), = (k for k, (v, _) in state_basis.items() if v is Ce_A)
    (x1, dx1), = (k for k, (v, _) in state_basis.items() if v is Ce_B)
    (x2, dx2), = (k for k, (v, _) in state_basis.items() if v is Ce_C)
    E = sympy.S("E")
    solutions = {dx1 + E * x1 - x2, dx2 - E * x1 + x2, x0 - E, dx0}

    relations = set(system.constitutive_relations)

    assert not solutions ^ relations


def biochemical_cycle(name="Cycle"):
    # Reactions: X = Y = Z = X
    R = 1
    T = 1

    model = bgt.new(name=name)

    X = bgt.new("Ce", name="X", library="BioChem",
                value={'k': 1, 'R': R, 'T': T})
    Y = bgt.new("Ce", name="Y", library="BioChem",
                value={'k': 2, 'R': R, 'T': T})
    Z = bgt.new("Ce", name="Z", library="BioChem",
                value={'k': 3, 'R': R, 'T': T})
    common_X = bgt.new("0", name="X")
    common_Y = bgt.new("0", name="Y")
    common_Z = bgt.new("0", name="Z")
    r1 = bgt.new("Re", name="r1", library="BioChem",
                 value={'r': 1, 'R': R, 'T': T})
    r2 = bgt.new("Re", name="r2", library="BioChem",
                 value={'r': 2, 'R': R, 'T': T})
    r3 = bgt.new("Re", name="r3", library="BioChem",
                 value={'r': 3, 'R': R, 'T': T})
    bgt.add(model, X, Y, Z, common_X, common_Y, common_Z, r1, r2, r3)

    bgt.connect(common_X, X)
    bgt.connect(common_X, r1)
    bgt.connect(r1, common_Y)
    bgt.connect(common_Y, Y)
    bgt.connect(common_Y, r2)
    bgt.connect(r2, common_Z)
    bgt.connect(common_Z, Z)
    bgt.connect(common_Z, r3)
    bgt.connect(r3, common_X)

    return model


def test_closed_cycle():
    model = biochemical_cycle("Closed cycle")

    assert len(model.state_vars) == 3
    assert len(model.control_vars) == 0

    X = model / "C:X"
    Y = model / "C:Y"
    Z = model / "C:Z"
    state_basis, _, _ = model.basis_vectors
    (x0, dx0), = (k for k, (v, _) in state_basis.items() if v is X)
    (x1, dx1), = (k for k, (v, _) in state_basis.items() if v is Y)
    (x2, dx2), = (k for k, (v, _) in state_basis.items() if v is Z)

    solutions = {
        dx0 + 4 * x0 - 2 * x1 - 9 * x2,
        dx1 - x0 + 6 * x1 - 6 * x2,
        dx2 - 3 * x0 - 4 * x1 + 15 * x2
    }
    relations = set(model.constitutive_relations)
    assert not solutions ^ relations


def open_cycle():
    model = biochemical_cycle("Open cycle")
    r1 = model / ("R:r1")
    r3 = model / ("R:r3")
    common_X = model / ("0:X")
    bgt.disconnect(common_X, r1)
    bgt.disconnect(r3, common_X)

    R = 1
    T = 1
    K_A = 1
    x_A = 1
    K_B = 2
    x_B = 1

    A = bgt.new("Se", name="A", value={'e': R * T * np.log(K_A * x_A)})
    B = bgt.new("Se", name="B", value={'e': R * T * np.log(K_B * x_B)})
    XA = bgt.new("1", name="XA")
    XB = bgt.new("1", name="XB")
    bgt.add(model, A, B, XA, XB)

    bgt.connect(common_X, XA)
    bgt.connect(A, XA)
    bgt.connect(XA, r1)
    bgt.connect(r3, XB)
    bgt.connect(XB, common_X)
    bgt.connect(XB, B)

    return model


def test_open_cycle():
    model = open_cycle()

    assert len(model.state_vars) == 3
    assert len(model.control_vars) == 0

    X = model / "C:X"
    Y = model / "C:Y"
    Z = model / "C:Z"
    state_basis, _, _ = model.basis_vectors
    (x0, dx0), = (k for k, (v, _) in state_basis.items() if v is X)
    (x1, dx1), = (k for k, (v, _) in state_basis.items() if v is Y)
    (x2, dx2), = (k for k, (v, _) in state_basis.items() if v is Z)

    solutions = {
        dx0 + 7 * x0 - 2 * x1 - 9 * x2,
        dx1 - x0 + 6 * x1 - 6 * x2,
        dx2 - 6 * x0 - 4 * x1 + 15 * x2
    }
    relations = set(model.constitutive_relations)
    assert not solutions ^ relations


class TestReactionNames(object):
    def test_naming(self):
        # see issue 79
        from BondGraphTools.reaction_builder import Reaction_Network
        from BondGraphTools.atomic import SymmetricComponent
        rn = Reaction_Network(name="Reaction 1")
        rn.add_reaction("A+B=C", name="E1")
        model = rn.as_network_model()

        assert model.name == "Reaction 1"

        re = model / "E1"

        assert isinstance(re, SymmetricComponent)
        assert re.metamodel == 'R'
