

def from_file(filename):
    pass

class BondGraph(object):
    _n = 0

    def __init__(self):

        self._nodes = {}
        self.ports = {}
        self.name = "Untitled Bond Graph"
        self.local_params = {}
        self.global_params = {}
        self.bonds = {}

    def add_node(self, node, node_id=None):

        if not node_id or node_id not in self._nodes:
            node_id = self._n
            self._n += 1

        self._nodes[node_id] = node

        return node_id


class NodeBase(object):
    def __init__(self, name, node_type, ports=None,
                 pos=None, local_params=None, global_params=None):
        self.ports = {}
        """ 
        port_id: {domain, [domain_restrictions]}
        This should be handled by inherited classes.
        """
        self.node_type = node_type
        """(str) Type classifier for this node"""
        self.name = name
        """(str) Name of this particular node"""
        self.local_params = self.build_local_params(local_params)
        """Dictionary of local parameters; of the form `param_str:param_dict`
        """
        self.global_params = global_params
        """List of strings, where each string is a model parameter that is 
        specified outside this class"""

        self.pos = pos
        """The `(x,y)` position of the center of this object"""

    def build_local_params(self, local_params):
        return None

    def __str__(self):
        return "{node_type}: {name}".format(
            node_type=self.node_type, name=self.name
        )


class AtomicNode(NodeBase):
    def __init__(self, *args, constitutive_relations=None, ports=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.constitutive_relations = constitutive_relations

        if ports:
            self.build_port_list(ports)

    def build_port_list(self, ports):
        pass

class CompositeNode(NodeBase):
    def __init__(self, *args, bond_graph, port_map, **kwargs):
        super().__init__(*args, **kwargs)

        self.bond_graph = bond_graph
        """Internal bond graph representation of the component"""
        self.port_map = port_map
        """Mapping between the ports exposed to the outer model, and the ports
        contained inside the internal bond graph representation"""