import pytest

from BondGraphTools.core.component_manager import base_id, get_components_list
from BondGraphTools.core.bondgraph import DeletionError

@pytest.mark.usefixture("BondGraph")
class TestBasicFunctionality(object):
    def test_add_components(self, BondGraph):

        component_list = get_components_list(base_id)
        assert component_list

        for comp_id,_ in component_list:
            BondGraph.add_component(comp_id)

        assert len(component_list) == len(BondGraph.nodes)

    def test_add_bonds(self, BondGraph):

        BondGraph.add_bond(0, 1)
        assert BondGraph.bonds == {(0, 1, 0, 0): 1}

    def test_remove_connected_component(self, BondGraph):
        with pytest.raises(DeletionError):
            BondGraph.remove_component(0)

    def test_remove_bond(self, BondGraph):
        for bond in BondGraph.find_bonds(1):
            BondGraph.remove_bond(bond)
        assert not BondGraph.bonds

    def test_remove_free_component(self, BondGraph):

        node = BondGraph.nodes[0]
        BondGraph.remove_component(0)

        with pytest.raises(KeyError):
            node_id = BondGraph.nodes[0]

        node_id = BondGraph.find_component(node=node)

        assert node_id == None


