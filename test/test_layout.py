import pytest
from BondGraphTools.view import _build_graph


@pytest.mark.usefixture("rlc")
def test_edge_list(rlc):

    graph = _build_graph(rlc)

    assert graph.nnz == 6


