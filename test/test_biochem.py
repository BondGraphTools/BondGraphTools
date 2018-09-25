import pytest

import sympy

import BondGraphTools as bgt
from BondGraphTools import connect
from BondGraphTools.exceptions import InvalidPortException
from BondGraphTools.algebra import extract_coefficients, inverse_coord_maps,get_relations_iterator
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

    r = Re.constitutive_relations[1]

    coords = list(sympy.sympify("e_0,f_0,e_1,f_1,r"))
    local_map = {c:i for i,c in enumerate(coords)}

    coeff_dict, nlin = extract_coefficients(r, local_map, coords)

    assert coeff_dict == {1: 1}
    assert nlin

    mappings, coords = inverse_coord_maps(*Re.basis_vectors)
    assert mappings
    assert coords

    for r in get_relations_iterator(Re, mappings, coords):
        assert r in [
            ({1:1, 3:1}, 0), ({1:1}, sympy.sympify("-r*exp(e_0) + r*exp(e_1)"))
        ]

def test_a_to_b_model():
    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])

    Re = bgt.new("Re", library="BioChem", value={'r': 1, "R": 1, "T": 1})

    Y_A = bgt.new('1')
    Y_B = bgt.new('1')

    a_to_b = bgt.new()
    a_to_b.add(A , Re, B,Y_A, Y_B)

    connect(A, Y_A.non_inverting)
    connect(B, Y_B.non_inverting)
    connect((Re,0), Y_A.inverting)
    connect((Re,1), Y_B.inverting)

    eqns ={
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

    state_basis ,_,_ = bg.basis_vectors
    (x_0, dx_0), = (k for k, (v, _) in state_basis.items() if v is A)
    (x_1, dx_1), = (k for k, (v, _) in state_basis.items() if v is B)
    (x_2, dx_2), = (k for k, (v, _) in state_basis.items() if v is C)

    eqns = {dx_0 + x_0*x_1 - x_2,
            dx_1 + x_0*x_1 - x_2,
            dx_2 - x_0*x_1 + x_2}

    relations = set(bg.constitutive_relations)
    assert not relations ^ eqns


def test_new_reaction_network():

    rn = Reaction_Network(name="A+B to C")
    rn.add_reaction(
        "A+B=C"
    )
    assert rn.species ==["A","B","C"]
    assert rn.forward_stoichiometry == sympy.Matrix([[1],[1],[0]])
    assert rn.reverse_stoichiometry == sympy.Matrix([[0], [0], [1]])
    assert rn.stoichiometry == sympy.Matrix([[-1],[-1],[1]])


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


#TODO: fix when we rework parameters
def test_nlin_se():
    rn = Reaction_Network(name="A+B to C", reactions="A+B=C")
    system = rn.as_network_model(normalised=True)

    for param in system.params:
        system.set_param(param, 1)

    Ce_A = system/ "A"
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
    state_basis , _, _ = system.basis_vectors

    (x0, dx0), = (k for k, (v,_) in  state_basis.items() if v is Ce_A)
    (x1, dx1), = (k for k, (v,_) in  state_basis.items() if v is Ce_B)
    (x2, dx2), = (k for k, (v,_) in  state_basis.items() if v is Ce_C)
    E = sympy.S("E")
    solutions = {dx1 + E * x1 - x2, dx2 - E * x1 + x2, x0 - E, dx0}

    relations = set(system.constitutive_relations)

    assert not solutions ^ relations


