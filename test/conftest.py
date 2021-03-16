import pytest
import BondGraphTools as bgt


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",
                     default=False, help="run slow tests")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark tests as slow")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture(scope='class')
def rlc():
    r = bgt.new("R", value=1)
    l = bgt.new("I", value=1)
    c = bgt.new("C", value=1)
    kvl = bgt.new("0", name="kvl")

    rlc = bgt.new()
    rlc.add([r, l, c, kvl])

    bgt.connect(r, kvl)
    bgt.connect(l, kvl)
    bgt.connect(c, kvl)

    return rlc
