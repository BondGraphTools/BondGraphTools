import BondGraphTools as bgt


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


