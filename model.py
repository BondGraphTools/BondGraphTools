import math
import random

import component_manager as cm


def from_file(filename):
    pass


class BondGraph(object):
    _n = 0

    def __init__(self, name=None):

        self.nodes = {}
        self.ports = {}
        self.name = name if name else "Untitled Bond Graph"
        self.local_params = {}
        self.global_params = {}
        self.bonds = dict()

    def add_component(self, component,
                      pos=None, name=None, library=None, node_id=None):

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
            lib, comp = cm.find(component)
            build_args = cm.get_component_data(lib, comp)
        else:
            build_args = cm.get_component_data(library, component)

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

        if not from_port:
            from_port = self.nodes[fr].next_port()
        if not to_port:
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

        return node_id


class NodeBase(object):
    factories = {}

    def __init__(self,
                 name,
                 node_type,
                 ports=None,
                 pos=None,
                 local_params=None,
                 global_params=None,
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

        self._free_ports = list(self.ports) if self.ports else []

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
    def factory(cls, *args, **kwargs):

        C = find_subclass(cls, NodeBase)
        if not C:
            raise NotImplementedError(
                "Node class not found", cls
            )

        return C(*args, **kwargs)


def find_subclass(name, base_class):

    for c in base_class.__subclasses__():
        if c.__name__ == name:
            return c
        else:
            sc = find_subclass(name, c)
            if sc:
                return sc
    return None


class AtomicNode(NodeBase):
    def __init__(self, *args, constitutive_relations=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.constitutive_relations = constitutive_relations


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


class OnePort(AtomicNode):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        for p in self.ports:
            assert isinstance(p, int)


def _optimise(nodes, edges, N=1000):
    """
    e is a list of edge tuples [(i,j), (i,m), ... ]


    Returns:

    """
    dirs = [(a,b) for a in range(-1,2) for b in range(-1,2)]

    n = 0
    obj = 1000
    targets = set()
    w1 = 0.5
    x, y = zip(*nodes)
    T = 1
    xt = [p for p in x]
    yt = [p for p in y]
    lenx = len(x)
    step = math.ceil(math.sqrt(lenx))
    while n < N:
        # generate new candidate
        dist = 1 + step*int(math.exp(-n))

        for i in range(lenx):
            if i in targets:
                xt[i] = x[i] + random.randint(-step, step)
                yt[i] = y[i] + random.randint(-step, step)
            else:
                d = random.randint(0, 8)
                dx, dy = dirs[d]
                xt[i] = x[i] + dx*dist
                yt[i] = y[i] + dy*dist

        # make sure we're not sitting on top of eachother
        if any(xt[i] == xt[j] and yt[i] == yt[j]
                for i in range(lenx) for j in range(lenx) if i != j):
            continue

        targets = set()
        crossings = 0
        dist = 0
        node_dist = sum(
            sum(((x1-x2)**2 + (y1 - y2)**2)**0.5 for x1, y1 in zip(xt,yt))
                for x2,y2 in zip(xt,yt))

        bd_x = [xt[j] - xt[i] for i, j in edges]
        bd_y = [yt[j] - yt[i] for i, j in edges]

        for k, (i, j) in enumerate(edges):
            bx = bd_x[k]
            by = bd_y[k]

            dist += bx**2 + by**2

            for l in range(k, len(edges)):
                ip, jp = edges[l]

                rs = bd_x[k]*bd_y[l] - bd_x[l]*bd_y[k]
                if rs != 0:
                    t = ((x[ip] - x[i])*bd_y[l] - (y[ip] - y[i])*bd_x[l])/rs
                    u = ((x[ip] - x[i])*bx - (y[ip] - y[i])*by)/rs
                    if 0 <= u <= 1 and 0 <= t <= 1:
                        crossings += 1
                        targets &= {i, j, ip, jp}

        new_obj = dist**2*(1 + crossings)**2 + w1 * node_dist
        delta = new_obj - obj
        if delta <= 0 or math.exp(-delta/T) > random.uniform(0, 1):
            x = [p for p in xt]
            y = [p for p in yt]
            obj = new_obj
        n += 1

    xm = int(sum(x)/lenx)
    ym = int(sum(y)/lenx)
    out = [(xp - xm, yp - ym) for xp, yp in zip(x,y)]

    return out


def arrange(bond_graph):
    nodes = list()
    for node in bond_graph.nodes.values():
        try:
            x, y = node.pos
            nodes.append((x, y))
        except TypeError:
            nodes.append((0, 0))
    edges = set()
    for i, j, _, _ in bond_graph.bonds:
        if i < j:
            edges.add((i, j))
        else:
            edges.add((j, i))
    pos = _optimise(nodes, list(edges))
    for n, node in enumerate(bond_graph.nodes.values()):
        node.pos = pos[n]