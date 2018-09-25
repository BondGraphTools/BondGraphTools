import pytest
import sympy
import BondGraphTools as bgt
import BondGraphTools.sim_tools as sim


@pytest.mark.use_fixture("rlc")
def test_build(rlc):
    assert len(rlc.state_vars) == 2
    assert len(rlc.ports) == 0

def test_build_rlc():
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    kvl = bgt.new("0", name="kvl")
    rlc = bgt.new()
    rlc.add([r, l, c, kvl])

    bgt.connect(r, kvl)
    bgt.connect(l, kvl)
    bgt.connect(c, kvl)
    assert len(kvl.ports) == 3


@pytest.mark.use_fixture("rlc")
def test_build_and_drive(rlc):
    se = bgt.new("Se")
    assert len(se.control_vars) == 1
    rlc.add(se)

    for comp in rlc.components:
        if comp.metaclass == "0":
            bgt.connect(se, comp)
            break

    assert len(rlc.bonds) == 4
    assert len(rlc.control_vars) == 1


# TODO: we should fix this when we rethink how parameters are set up.
@pytest.mark.skip
def test_symbolic_params():
    r = bgt.new("R", value=sympy.symbols('r'))
    l = bgt.new("I", value=sympy.symbols('l'))
    c = bgt.new("C", value=sympy.symbols('c'))
    kvl = bgt.new("0", name="kvl")
    rlc = bgt.new()
    rlc.add([r,l, c , kvl])

    bgt.connect(r, kvl)
    bgt.connect(l, kvl)
    bgt.connect(c, kvl)

    assert len(rlc.params) == 3

    params = {k for _, k in rlc.params.values()}

    assert params & set(sympy.symbols('r, l, c'))


@pytest.mark.use_fixture("rlc")
def test_rlc_con_rel(rlc):

    rel = rlc.constitutive_relations

    _, v = rlc.state_vars['x_0']


    if str(v) != 'q_0':
        eq1 = sympy.sympify("dx_0 - x_1")
        eq2 = sympy.sympify("dx_1 + x_0 + x_1")
    else:
        eq1 = sympy.sympify("dx_1 - x_0")
        eq2 = sympy.sympify("dx_0+ x_0 + x_1")

    for r in rel:
        assert r in (eq1, eq2)

    assert "x_0" in rlc.state_vars
    assert "x_1" in rlc.state_vars

def test_tf():
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    tf = bgt.new("TF", value=0.5)
    tflc = bgt.new()
    tflc.add([tf,l, c])



    bgt.connect(l, (tf, 1))
    bgt.connect(c, (tf, 0))

    c,m,lp,nlp,const = tflc.system_model()
    assert nlp.is_zero
    assert const ==[]


def test_se():

    Se = bgt.new('Se', value=1)
    c = bgt.new('C', value=1)
    vc = bgt.new()
    vc.add([Se, c])
    assert Se.constitutive_relations == [sympy.sympify("e_0 - 1")]
    bgt.connect(Se, c)


    assert vc.constitutive_relations == [sympy.sympify("dx_0"),
                                         sympy.sympify("x_0 - 1")]

def test_one():
    loop_law = bgt.new('1')
    Se = bgt.new('Se', value=1)
    c = bgt.new('C', value=1)
    r = bgt.new('R', value=1)
    vc = bgt.new()
    vc.add([Se, c, loop_law, r])

    bgt.connect(Se, (loop_law, loop_law.non_inverting))
    bgt.connect(c, (loop_law, loop_law.inverting))
    bgt.connect(r, (loop_law, loop_law.inverting))

    assert vc.constitutive_relations == [
        sympy.sympify("dx_0 + x_0 - 1")
    ]


