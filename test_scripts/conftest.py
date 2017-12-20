import pytest
import sys
import BondGraphTools


@pytest.fixture(scope="class")
def BondGraph():
    return BondGraphTools.BondGraph()
