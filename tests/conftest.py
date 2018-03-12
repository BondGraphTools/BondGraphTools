import pytest

import BondGraphTools


@pytest.fixture(scope='class')
def BondGraph():
    return BondGraphTools.BondGraph("Test")

@pytest.fixture(scope='class')
def RLC():
    bg = BondGraphTools.BondGraph("RLC")

    Se = bg.add_component("Se", name="Vs")
    c = bg.add_component("C", name="Capacitor")
    i = bg.add_component("I", name="Inductor")
    r = bg.add_component("R", name="Resistor")
    kvl = bg.add_component("0", pos=(0, 0))

    bg.add_bond(kvl, c)
    bg.add_bond(kvl, i)
    bg.add_bond(kvl, Se)
    bg.add_bond(kvl, r)
