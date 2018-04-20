import os
import pytest

import BondGraphTools.component_manager as cm


def test_load():

    dir, _ = os.path.split(__file__ )
    file = os.path.join(dir, "lib/test_lib.json")

    assert cm.load_library(file)
    assert "test" in cm.__libraries


def test_duplicate_load():

    del cm.__libraries["test"]
    dir, _ = os.path.split(__file__ )
    file = os.path.join(dir, "lib/test_lib.json")

    assert cm.load_library(file)
    assert "test" in cm.__libraries

    assert not cm.load_library(file)


def test_bad_libary():
    dir, _ = os.path.split(__file__)
    file = os.path.join(dir, "lib/test_lib2.json")

    assert not cm.load_library(file)


