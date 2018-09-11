import pytest

import sympy

import BondGraphTools as bgt
from BondGraphTools import connect
from BondGraphTools.exceptions import InvalidPortException
from BondGraphTools.algebra import extract_coefficients, inverse_coord_maps,get_relations_iterator
from BondGraphTools.reaction_builder import Reaction_Network

@pytest.mark.skip
def test_make_a_to_b():

    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])

    Re = bgt.new("Re", library="BioChem", value={"R": 1, "T": 1})
    a_to_b = bgt.new()
    a_to_b.add(A, Re, B)
    with pytest.raises(InvalidPortException):
        connect(A, Re)

    assert (Re,0) in Re.ports
    connect(A, (Re, 0))
    connect(B, (Re, 1))

    assert Re.control_vars == ['r']
    assert list(a_to_b.state_vars.keys()) == ['x_0', 'x_1']
    assert list(a_to_b.control_vars.keys()) == ['u_0']
    assert not a_to_b.ports

@pytest.mark.skip
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

    connect(A, Y_A.input)
    connect(B, Y_B.input)
    connect((Re,0), Y_A.output)
    connect((Re,1), Y_B.output)

    eqns ={
        sympy.sympify("dx_0 + x_0 -x_1"), sympy.sympify("dx_1 + x_1 -x_0")
    }
    for relation in a_to_b.constitutive_relations:
        assert relation in eqns

# TODO: we need to find a better way to compare equations...
@pytest.mark.skip
def test_ab_to_c_model():

    A = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    B = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    C = bgt.new("Ce", library="BioChem", value=[1, 1, 1])
    Re = bgt.new("Re", library="BioChem", value={'r': 1, "R": 1, "T": 1})
    Y_AB = bgt.new('Y', library="BioChem")
    Y_C = bgt.new('Y', library="BioChem")

    assert Y_AB.ports[0]
    assert Y_AB._fixed_ports == {0}
    bg = bgt.new()
    bg.add(A, B, Re, Y_AB, C, Y_C)
    connect(A, (Y_AB, 1))

    assert Y_AB.ports[1] == 1
    connect(B, (Y_AB, 2))
    connect((Re, 0), (Y_AB, 0))
    connect((Re, 1), (Y_C, 0))
    connect((Y_C, 1), C)

    assert 0 in Y_AB._fixed_ports
    assert 0 in Y_AB.ports

    eqns = {
        sympy.sympify("dx_0 + x_0*x_1 - x_2"),
        sympy.sympify("dx_1 + x_0*x_1 - x_2"),
        sympy.sympify("dx_2 - x_0*x_1 + x_2")
    }

    for relation in bg.constitutive_relations:
        assert relation in eqns

@pytest.mark.skip
def test_new_reaction_network():

    rn = Reaction_Network(name="A+B to C")
    rn.add_reaction(
        "A+B=C"
    )
    assert rn.species ==["A","B","C"]
    assert rn.forward_stoichiometry == sympy.Matrix([[1],[1],[0]])
    assert rn.reverse_stoichiometry == sympy.Matrix([[0], [0], [1]])
    assert rn.stoichiometry == sympy.Matrix([[-1],[-1],[1]])

@pytest.mark.skip
def test_rn_to_bond_graph():
    rn = Reaction_Network(name="A+B to C", reactions="A+B=C")

    system = rn.as_network_model(normalised=True)
    for comp in system.components:
        if comp.metaclass == 'Y':
            for port in comp.ports:
                if port != 0:
                    assert comp.ports[port] == 1

    assert len(system.state_vars) == 3
    assert len(system.control_vars) == 4

@pytest.mark.skip
def test_cat_rn():
    reactions = [
        "E + S = ES", "ES = E+P"
    ]
    rn = Reaction_Network(reactions=reactions, name="Catalysed Reaction")
    assert rn._species["ES"] == 2
    assert len(rn._reactions) == 2

#TODO: fix when we rework parameters
@pytest.mark.skip
def test_nlin_se():
    rn = Reaction_Network(name="A+B to C", reactions="A+B=C")
    system = rn.as_network_model(normalised=True)

    for param in system.params:
        system.set_param(param, 1)

    Ce_A = system.find(name="A", component="Ce")
    Ce_B = system.find(name="B", component="Ce")

    assert Ce_A is not Ce_B

    Y = system.find(component="Y", name="AB")

    system.disconnect(Ce_A, Y)

    J_A = bgt.new("0", name="Ce_A")
    Se_A = bgt.new('Se', value=1, name='1')

    system.add(J_A),
    system.add(Se_A)
    connect(J_A, Y)
    connect(Ce_A, J_A)
    connect(Se_A, J_A)

    eqns= "dx_0, dx_1 + E * x_1 - x_2,dx_2 - E * x_1 + x_2, x_0 - E"
    relations = set(system.constitutive_relations)
    solutions = set(sympy.sympify(eqns))

    assert not solutions ^ relations


