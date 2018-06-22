import pytest
import sympy

import BondGraphTools as bgt
from BondGraphTools.algebra import extract_coefficients, smith_normal_form, \
    adjacency_to_dict, augmented_rref,_generate_substitutions


def test_extract_coeffs_lin():
    eqn = sympy.sympify("y -2*x -3")
    local_map = {
        sympy.symbols("y"): 0,
        sympy.symbols("x"): 1
    }
    coords = [sympy.symbols("r_1"), sympy.symbols("r_0"), sympy.S(1)]

    lin, nlin = extract_coefficients(eqn, local_map, coords)
    assert lin[1] == -2
    assert lin[0] == 1
    assert nlin == sympy.sympify("-3")


def test_extract_coeffs_nlin():
    eqn = sympy.sympify("y -2*x -3 + exp(x)")
    local_map = {
        sympy.symbols("y"): 0,
        sympy.symbols("x"): 1
    }
    coords = [sympy.symbols("r_1"), sympy.symbols("r_0"), sympy.S(1)]

    lin, nlin = extract_coefficients(eqn, local_map, coords)
    assert lin[1] == -2
    assert lin[0] == 1
    assert nlin == sympy.sympify("exp(r_0) - 3")


def test_smith_normal_form():

    m = sympy.SparseMatrix(2,3,{(0,2):2, (1,1):1})
    mp,_,_ = smith_normal_form(m)
    assert mp.shape == (3, 3)
    assert mp[2, 2] != 0




def test_smith_normal_form_2():
    matrix = sympy.eye(3)
    matrix.row_del(1)

    m,_,_ = smith_normal_form(matrix)

    diff = sympy.Matrix([[0,0,0],[0,1,0], [0,0,0]])

    assert (sympy.eye(3) - m) == diff


def test_smith_normal_form_3():
    m = sympy.SparseMatrix(5,3,{(0,1):1, (1,0):1,
                                (4,2):1})
    mp,_,_ = smith_normal_form(m)
    assert mp.shape == (3, 3)
    assert (mp - sympy.eye(3)).is_zero


def test_aug_rref():
    matrix = sympy.Matrix([
        [0, 0, 0, 1],
        [1, 0, 1, 0],
        [0, 2, 0, 0],
        [0, 0, 4, 0]
    ])

    adj = sympy.eye(4)
    M = matrix.row_join(adj)
    assert M.shape == (4, 8)
    MA = augmented_rref(M, 4)

    m = MA[:, :4]
    a = MA[:, 4:]
    assert m is not matrix
    assert a is not adj

    assert m != matrix
    assert a != adj

    assert (m - sympy.eye(4)).is_zero

    assert (a * matrix - sympy.eye(4)).is_zero


def test_augmented_rref():
    M = sympy.Matrix(
        [[1, 0, 1, 0],
         [1, 0, 1, 0],
         [0, 0, 1 ,0]])

    A = sympy.Matrix(3,1, list(sympy.symbols('a,a,c')))
    target_A = sympy.Matrix(3, 1,
                            [sympy.sympify('a - c'), sympy.symbols('c'),0])
    Mr, _ = M.rref()

    assert Mr == sympy.Matrix([[1,0,0,0],
                               [0,0,1,0],
                               [0,0,0,0]])

    MA = M.row_join(A)
    MA = augmented_rref(MA, 1)
    assert MA[:,:-1]== Mr
    assert target_A == MA[:,-1:]


def test_augmented_rref_2():
    M = sympy.Matrix(
        [[1, 0, 1, 0],
         [1, 0, 1, 0],
         [0, 0, 1 ,0]])

    A = sympy.Matrix(3,1, [sympy.symbols('a'),sympy.symbols('a'), 1])
    target_A = sympy.Matrix(3, 1,
                            [sympy.sympify('a - 1'), 1, 0])
    Mr, _ = M.rref()

    assert Mr == sympy.Matrix([[1,0,0,0],
                               [0,0,1,0],
                               [0,0,0,0]])
    assert M.shape == (3,4)
    assert A.shape == (3, 1)
    MA = M.row_join(A)
    MA = augmented_rref(MA, 1)
    assert MA[:, :-1] == Mr
    assert target_A == MA[:, -1:]


def test_build_relations():
    c = bgt.new("C")

    eqns = c._build_relations()
    assert len(eqns) == 1

    eqn = eqns.pop()

    test_eqn = sympy.sympify("q_0 - C*e_0")

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
    index_map = {(c,0):0, (kvl,0):1}
    M = adjacency_to_dict(index_map, bg.bonds, offset=1)
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
    assert len(ports) == 0


def test_relations_iter():
    c = bgt.new("C", value=1)

    mappings = ({(c, 'q_0'): 0}, {(c,0): 0}, {})
    coords = list(sympy.sympify("dx_0,e_0,f_0,x_0, 1"))
    relations = c.get_relations_iterator(mappings, coords)

    d1 = {0:1, 2:-1} #dx/dt - f = 0
    d2 = {1:1, 3:-1} #e - x = 0$
    d1m = {0:-1, 2:1} #dx/dt - f = 0
    d2m = {1:-1, 3:1} #e - x = 0$

    for lin, nlin in relations:
        assert not nlin
        assert lin in (d1, d2, d1m, d2m)


def test_relation_iter_twoport():
    tf = bgt.new("TF", value=1)
    mappings = ({},{(tf,0):0, (tf,1):1},{})
    coords = list(sympy.sympify("e_0,f_0,e_1,f_1"))
    # d1 = {0: 1, 2: -1}
    # d2 = {1: 1, 2: 1}
    relations = tf.get_relations_iterator(mappings, coords)
    for lin, nlin in relations:
        assert not nlin


@pytest.mark.usefixture("rlc")
def test_interal_basis_vectors(rlc):
    tangent, ports, cv = rlc._build_internal_basis_vectors()

    assert len(tangent) == 2
    assert len(cv) == 0
    assert len(ports) == 6


def test_cv_relations():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = c + se + kcl + r

    bg.connect(c,kcl)
    bg.connect(se, kcl)
    bg.connect(r, kcl)

    assert bg.constitutive_relations == [sympy.sympify("dx_0 + u_0 + x_0")]


def test_parallel_crv_relations():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("0")
    bg = c + se + kcl + r

    bg.connect(c, kcl)
    bg.connect(se, kcl)
    bg.connect(r, kcl)

    assert bg.constitutive_relations == [sympy.sympify("dx_0 - du_0"),
                                         sympy.sympify("x_0 - u_0")]


def test_generate_subs():

    w, x, y, z = sympy.sympify("w,x,y,z")
    size_tuple  =(0, 2, 0,4 )
    coords = [w, x, y, z]
    #  w + w^2 + x^2
    #  x + 1 + y^2   < should appear in subs as x = -y^2  - 1
    #  y + 1         <                          y = - 1
    #  0 + z^2 + w^2

    L = sympy.SparseMatrix(4, 4, {(0,0): 1,
                                  (1,1): 1,
                                  (2,2): 1})

    N = sympy.SparseMatrix(4, 1, {(0, 0): w**2 + x**2,
                                  (1, 0): 1 + y**2,
                                  (2, 0): 1})

    constraint = [z**2 + w**2]
    subs = _generate_substitutions(L, N,constraint, coords, size_tuple)
    target_subs = [(y,-1), (x, -1-y**2)]

    assert subs == target_subs
