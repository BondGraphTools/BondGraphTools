
import sympy
import BondGraphTools as bgt
import BondGraphTools.sim_tools as sim


def build_rlc():
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    kvl = bgt.new("0", name="kvl")

    rlc = r + l + c + kvl

    rlc.connect(r, kvl)
    rlc.connect(l, kvl)
    rlc.connect(c, kvl)

    return rlc


def test_build():
    rlc = build_rlc()

    assert len(rlc.state_vars) == 2
    assert len(rlc.ports) == 0


def test_build_and_drive():
    rlc = build_rlc()
    se = bgt.new("Se")
    assert len(se.control_vars) == 1
    rlc += se

    for comp in rlc.components.values():
        if comp.type == "0":
            rlc.connect(se, comp)
            break

    assert len(rlc.bonds) == 4
    assert len(rlc.control_vars) == 1


def test_symbolic_params():
    r = bgt.new("R", value=sympy.symbols('r'))
    l = bgt.new("I", value=sympy.symbols('l'))
    c = bgt.new("C", value=sympy.symbols('c'))
    kvl = bgt.new("0", name="kvl")

    rlc = r + l + c + kvl

    rlc.connect(r, kvl)
    rlc.connect(l, kvl)
    rlc.connect(c, kvl)

    assert len(rlc.params) == 3

    assert set(rlc.params.values()) & set(sympy.symbols('r, l, c'))


def test_simulate():

    rlc = build_rlc()

    t, X = sim.simulate(
        rlc,
        timespan=[0, 10],
        initial_state={k: 1 for k in rlc.state_vars}
    )

