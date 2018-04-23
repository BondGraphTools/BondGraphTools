import BondGraphTools as bgt


def build_rlc():
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    kvl = bgt.new("0")

    rlc = r + l + c + kvl

    rlc.connect(r, kvl)
    rlc.connect(l, kvl)
    rlc.connect(c, kvl)

    return rlc


def test_build():

    rlc = build_rlc()

    assert len(rlc.state_vars) == 2
    assert len(rlc.ports) == 0