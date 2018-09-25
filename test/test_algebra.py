import pytest
import sympy
import logging

import BondGraphTools as bgt
from BondGraphTools import connect, new, expose
from BondGraphTools.algebra import extract_coefficients, smith_normal_form, \
    adjacency_to_dict, augmented_rref,_generate_substitutions,\
    inverse_coord_maps, _generate_cv_substitutions, get_relations_iterator

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

    test_eqn = {sympy.sympify("q_0 - C*e_0"),
                sympy.sympify("dq_0 - f_0")}

    assert set(eqns) == test_eqn

def test_zero_junction_relations():
    r = bgt.new("R", value=sympy.symbols('r'))
    l = bgt.new("I", value=sympy.symbols('l'))
    c = bgt.new("C", value=sympy.symbols('c'))
    kvl = bgt.new("0", name="kvl")

    rlc = bgt.new()
    rlc.add([c, l, kvl, r])


    connect(r, kvl)
    connect(l, kvl)
    connect(c, kvl)

    rels = kvl.constitutive_relations

    assert sympy.sympify("e_1 - e_2") in rels
    assert sympy.sympify("e_0 - e_2") in rels
    assert sympy.sympify("f_0 + f_1 + f_2") in rels


@pytest.mark.usefixture('rlc')
def test_basis_vectors(rlc):

    model_basis_vects = set()

    for component in rlc.components:
        for vects in component.basis_vectors:
            basis_vects = set(vects.values())

            assert not basis_vects & model_basis_vects

            model_basis_vects |= basis_vects


def test_build_junction_dict():
    c = bgt.new("C")
    kvl = bgt.new("0")
    bg = bgt.new()
    bg.add([c, kvl])
    connect(kvl, c)
    cp,kp = list(c.ports) + list(kvl.ports)
    index_map = {cp:0, kp:1}
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
    cp, = list(c.ports)
    mappings = ({(c, 'q_0'): 0}, {cp: 0}, {})
    coords = list(sympy.sympify("dx_0,e_0,f_0,x_0, 1"))
    relations = get_relations_iterator(c, mappings, coords)

    d1 = {0:1, 2:-1} #dx/dt - f = 0
    d2 = {1:1, 3:-1} #e - x = 0$
    d1m = {0:-1, 2:1} #dx/dt - f = 0
    d2m = {1:-1, 3:1} #e - x = 0$

    for lin, nlin in relations:
        assert not nlin
        assert lin in (d1, d2, d1m, d2m)


def test_relation_iter_twoport():
    tf = bgt.new("TF", value=1)
    tf_port_0, tf_port_1 = list(tf.ports)
    mappings = ({},{tf_port_0:0, tf_port_1:1},{})
    coords = list(sympy.sympify("e_0,f_0,e_1,f_1"))
    # d1 = {0: 1, 2: -1}
    # d2 = {1: 1, 2: 1}
    relations = get_relations_iterator(tf, mappings, coords)
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
    bg = bgt.new()
    bg.add([c, se, kcl, r])

    connect(c,(kcl,kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    assert bg.constitutive_relations == [sympy.sympify("dx_0 + u_0 + x_0")]


def test_parallel_crv_relations():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("0")
    bg = bgt.new()
    bg.add([c, se, kcl, r])

    connect(c, kcl)
    connect(se, kcl)
    connect(r, kcl)

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


def test_cv_subs_func():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = bgt.new()
    bg.add([c, se, kcl, r])

    connect(c,(kcl,kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    cv_s = {'u_0': ' -exp(-t)'}

    subs = [(sympy.Symbol('u_0'), sympy.sympify('-exp(-t)'))]

    mappings, coords = inverse_coord_maps(*bg.basis_vectors)
    assert _generate_cv_substitutions(cv_s, mappings,coords) == subs


def test_cv_subs_const():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = bgt.new()
    bg.add([c, se, kcl, r])

    connect(c,(kcl,kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    cv_s = {'u_0': ' 2'}

    subs = [(sympy.Symbol('u_0'), sympy.S(2))]

    mappings, coords = inverse_coord_maps(*bg.basis_vectors)
    assert _generate_cv_substitutions(cv_s, mappings,coords) == subs


def test_cv_subs_state_func():
    c = bgt.new("C", value=1)
    se = bgt.new("Se")
    r = bgt.new("R", value=1)
    kcl = bgt.new("1")
    bg = bgt.new()
    bg.add([c, se, kcl, r])

    connect(c,(kcl,kcl.non_inverting))
    connect(r, (kcl, kcl.non_inverting))
    connect(se, (kcl, kcl.non_inverting))

    cv_s = {'u_0': ' -exp(-x_0)'}

    subs = [(sympy.Symbol('u_0'), sympy.sympify('-exp(-x_0)'))]

    mappings, coords = inverse_coord_maps(*bg.basis_vectors)
    assert _generate_cv_substitutions(cv_s, mappings,coords) == subs



def test_ported_cr():
    model = bgt.new()
    Sf = bgt.new('Sf', name="Sf")
    R = bgt.new("R", value=2)
    zero = bgt.new("0")
    ss = bgt.new("SS")

    model.add(Sf, R, zero, ss)
    connect(Sf, zero)
    connect(R, zero)
    connect(ss, zero)

    bgt.expose(ss, 'A')
    assert len(model.control_vars) == 1

    ts, ps, cs = model._build_internal_basis_vectors()
    assert len(cs) == 1
    assert len(ps) == 7
    assert len(ts) == 0

    mapping, coords = inverse_coord_maps(ts, ps, cs)
    assert len(coords) == 15

    coords, mappings, lin_op, nl_op, conttr = model.system_model()
    assert nl_op.is_zero
    assert not conttr

    assert model.constitutive_relations == [
        sympy.sympify('e_0 - 2*f_0 - 2*u_0')
    ]

def test_ported_series_resistor():

    Se = new("Se")
    r1 = new("R", value=1)
    r2 = new("R", value=2)
    kvl = new('1')
    ss = new("SS")
    model = new()
    model.add(
        Se,r1,r2,kvl, ss
    )
    expose(ss)
    connect(Se, kvl.non_inverting)
    connect(kvl.inverting, r1)
    connect(kvl.inverting, r2)
    connect(kvl.inverting, ss)

    assert len(model.ports) == 1

    assert model.constitutive_relations == [
        sympy.sympify("e_0 - 3*f_0 - u_0")
    ]

def test_ported_cap():
    model = new()
    c = new("C", value=3)
    zero = new("0")
    ss = new("SS")
    model.add(
        c, zero, ss
    )

    connect(c, zero)
    connect(ss, zero)

    expose(ss)
    assert len(model.ports) == 1


    assert model.constitutive_relations == [
        sympy.sympify("dx_0 - f_0"),sympy.sympify("e_0 - x_0/3")
    ]
def test_ported_parallel_rc():

    model = new()
    r = new("R", value=2)
    c = new("C", value=3)
    zero = new("0")
    ss = new("SS")
    model.add(
        r,c,zero, ss
    )

    connect(r,zero)
    connect(c,zero)
    connect(ss, zero)

    expose(ss)
    assert len(model.ports) == 1

    assert model.constitutive_relations == [
        sympy.sympify("dx_0 + x_0/6 - f_0"),
        sympy.sympify("e_0 - x_0/3")
    ]
