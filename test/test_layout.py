import pytest

import BondGraphTools as bgt
from BondGraphTools.view import _build_graph, _metro_layout


@pytest.mark.usefixture("rlc")
def test_edge_list(rlc):

    graph = _build_graph(rlc)

    assert graph.nnz == 6

@pytest.mark.usefixture("rlc")
def test_metro_layout(rlc):

    graph = _build_graph(rlc)
    nodes = _metro_layout(graph)

