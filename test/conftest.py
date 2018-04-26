import pytest
import BondGraphTools as bgt

@pytest.fixture(scope='class')
def rlc():
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    kvl = bgt.new("0", name="kvl")

    rlc = r + l + c + kvl

    rlc.connect(r, kvl)
    rlc.connect(l, kvl)
    rlc.connect(c, kvl)

    return rlc