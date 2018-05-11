
import logging
import sympy as sp

from . import component_manager as cm

logger = logging.getLogger(__name__)

def from_file(filename):
    pass


class BondGraph(object):
    _n = 0

    def __init__(self, name=None):

        self.nodes = {}
        self.ports = {} # rename me
        self.name = name if name else "Untitled Bond Graph"
        # self.local_params = {}
        self.global_params = {}
        self.bonds = dict()

    @property
    def local_params(self):
        out = dict()
        for node in self.nodes.values():
            if node.local_params:
                for param, vals in node.local_params.items():
                    p_name = "{}_{{{}}}".format(param, node.name)
                    out[p_name] = vals
        return out

    def add_component(self, component,
                      pos=None, name=None, library=None, node_id=None, **kwargs):

        names = {node.name for node in self.nodes.values()
                 if node.node_type == component}
        if name and name in names:
            raise ValueError("{} already exists".format(name))
        elif not name:
            i = 1
            name = "{}_{}".format(component, i)
            while name in names:
                i += 1
                name = "{}_{}".format(component, i)

        if not library:
            library = cm.find(component, ensure_unique=True)
        build_args = cm.get_component_data(library, component)
        build_args.update(kwargs)

        if not build_args:
            raise NotImplementedError

        node = NodeBase.factory(**build_args, name=name,
                                pos=pos, node_type=component)

        return self._add_node(node, node_id)

    def add_bond(self, from_component, to_component,
                 from_port=None, to_port=None):

        if isinstance(from_component, str):
            fr = self.find_component(name=from_component)
        elif isinstance(from_component, int):
            fr = from_component
        elif isinstance(from_component, NodeBase):
            fr = self.find_component(node=from_component)
        else:
            raise NotImplementedError("Could not find base component")

        if isinstance(to_component, str):
            to = self.find_component(name=to_component)
        elif isinstance(to_component, int):
            to = to_component
        elif isinstance(to_component, NodeBase):
            to = self.find_component(node=to_component)
        else:
            raise NotImplementedError("Could not find destination component")

        if from_port != 0 and not from_port:
            from_port = self.nodes[fr].next_port()
        if to_port != 0 and not to_port:
            to_port = self.nodes[to].next_port()

        self.nodes[fr].reserve_port(from_port)

        try:
            self.nodes[to].reserve_port(to_port)
        except Exception as e:
            self.nodes[fr].release_port(from_port)
            raise e

        self.bonds[(from_component, to_component, from_port, to_port)] = 1

    def find_component(self, name=None, node_type=None, node=None):
        if node:
            for node_id, test_node in self.nodes.items():
                if test_node is node:
                    return node_id
        elif name:
            for node_id, test_node in self.nodes.items():
                if test_node.name == name:
                    if not node_type or node_type == test_node.node_type:
                        return node_id
        else:
            raise ValueError("Must specify search conditions")

        return None

    def _add_node(self, node, node_id=None):

        if not node_id or node_id not in self.nodes:
            node_id = self._n
            self._n += 1

        self.nodes[node_id] = node

        if isinstance(node, IOPort):
            p = 0
            while p < len(self.ports):
                p += 1
            self.ports[p] = node_id

        return node_id

    def find_bonds(self, node):
        if isinstance(node, NodeBase):
            node_id = self.find_component(node=node)
        elif isinstance(node, int):
            node_id = node

        out = []

        for bond in self.bonds:
            n1, n2, _, _ = bond

            if n1 == node_id or n2 == node_id:
                out.append(bond)

        return out

    def remove_component(self, node):
        if isinstance(node, NodeBase):
            node_id = self.find_component(node=node)
        elif isinstance(node, int):
            node_id = node
        else:
            raise NotImplementedError("Could not find node %s of type %s",
                                      node,
                                      node.__class__)

        if not self.find_bonds(node_id):

            del self.nodes[node_id]
            search = self.ports.items()
            for k, v in search:
                if v == node_id:
                    del self.ports[k]
        else:
            logger.warning("Could not delete; attached bond")
            raise DeletionError

    def remove_bond(self, bond):
        del self.bonds[bond]


class DeletionError(Exception):
    pass


class NodeBase(object):
    factories = {}

    def __init__(self,
                 name,
                 node_type,
                 ports=None,
                 pos=None,
                 local_params=None,
                 global_params=None,
                 control_variables=None,
                 string=None,
                 **kwargs):

        self.ports = None if not ports else {
            int(port): data for port, data in ports.items()
        }
        """ 
        port_id: {domain, [domain_restrictions]}
        This should be handled by inherited classes.
        """

        self.node_type = node_type
        """(str) Type classifier for this node"""

        self.name = name
        """(str) Name of this particular node"""

        self.local_params = local_params
        """Dictionary of local parameters; of the form `param_str:param_dict`
        """

        self.global_params = global_params
        """List of strings, where each string is a model parameter that is 
        specified outside this class"""

        self.pos = pos
        """The `(x,y)` position of the center of this object"""

        self.string = string if string else "{node_type}: {name}".format(
            node_type=self.node_type, name=self.name
        )
        self.control_variables = control_variables

        self._free_ports = list(self.ports) if self.ports else []

        for arg, val in kwargs.items():
            if self.local_params and arg in self.local_params:
                self.local_params[arg] = val

    def __repr__(self):
        return "\'{}: {}\'".format(self.node_type, self.name)

    def __str__(self):
        return self.string

    def next_port(self):
        try:
            return self._free_ports[0]
        except IndexError:
            return None

    def reserve_port(self, port):
        self._free_ports.remove(port)

    def release_port(self, port):
        if port not in self.ports:
            raise ValueError("Invalid Port")
        self._free_ports.append(port)

    @staticmethod
    def factory(cls=None, *args, **kwargs):

        if not kwargs["ports"]:
            C = ManyPort
        if not cls:
            C = AtomicNode
        else:
            C = _find_subclass(cls, NodeBase)

        if not C:
            raise NotImplementedError(
                "Node class not found", cls
            )

        return C(*args, **kwargs)


def _find_subclass(name, base_class):

    for c in base_class.__subclasses__():
        if c.__name__ == name:
            return c
        else:
            sc = _find_subclass(name, c)
            if sc:
                return sc

#    return AtomicNode


class AtomicNode(NodeBase):
    def __init__(self, *args, constitutive_relations=None, **kwargs):
        super().__init__(*args, **kwargs)
        if constitutive_relations:
            self.constitutive_relations = constitutive_relations
        else:
            self.constitutive_relations = []

    def _generate_local_coordinates(self):
        xstr = [s + str(port) for port in self.ports
                for s in ["p_", "q_", "e_", "f_"]]
        return sp.sympify(xstr)


class CompositeNode(NodeBase):
    def __init__(self, *args, bond_graph, port_map, **kwargs):
        super().__init__(*args, **kwargs)

        self.bond_graph = bond_graph
        """Internal bond graph representation of the component"""
        self.port_map = port_map
        """Mapping between the ports exposed to the outer model, and the ports
        contained inside the internal bond graph representation"""


class ManyPort(AtomicNode):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ports = {}

    def release_port(self, port):
        del self.ports[port]

    def reserve_port(self, port):
        if (self.ports and port in self.ports) or not isinstance(port, int):
            raise ValueError("Invalid Port {}".format(port))

        self.ports[port] = None

    def next_port(self):
        n = 0

        if not self.ports:
            return n

        while n in self.ports:
            n += 1
        return n

#
# class OnePort(AtomicNode):
#     pass
#
#
# class TwoPort(AtomicNode):
#     pass


class IOPort(AtomicNode):
    def __init__(self, *args, sockets=None, **kwargs):
        super().__init__(*args,  **kwargs)



