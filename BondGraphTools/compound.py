
from .base import BondGraphBase, InvalidPortException, \
    InvalidComponentException
from .view import GraphLayout


class BondGraph(BondGraphBase):
    def __init__(self, name, components=None, **kwargs):

        super().__init__(name, **kwargs)
        self.components = dict()

        if components:
            for component in components:
                self.__add__(component)

        self.bonds = list()
        self.type = "Composite"
        self.view = GraphLayout()
        """Graphical Layout of internal components"""

        self.cv_map = dict()
        self._port_map = dict()

    def save(self, filename):
        raise NotImplementedError

    def __getitem__(self, item):
        return self.components[item]

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self.components
        else:
            return item in self.components.values()

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __add__(self, other):

        if isinstance(other, BondGraphBase):
            n = 0
            for k in self.components.keys():
                t_idx, n_idx = k.split('_')

                n_idx = int(n_idx)
                if t_idx == other.type and n_idx >= n:
                    n = n_idx + 1

            self.components[f"{other.type}_{n}"] = other

        else:
            raise TypeError("unsupported operand type(s) for +: '%s' and '%s'",
                            type(self),
                            type(other))
        return self

    @property
    def params(self):
        return {
            (c,p): v.params[p] for c,v in self.components.items() if v.params
            for p in v.params if not isinstance(p, (float,complex, int))
        }

    @property
    def ports(self):
        bonds = {v for bond in self.bonds for v in bond}

        return {(c, p): v.ports[p]
                for c, v in self.components.items() for p in v.ports
                if (v, p) not in bonds}

    @property
    def state_vars(self):
        return {(v, i): v.state_vars[i] for
                c, v in self.components.items() if v.state_vars
                for i in v.state_vars}

    @property
    def control_vars(self):
        return [(v, i) for
                v in self.components.values() if v.control_vars
                for i in v.control_vars]

    def connect(self, source, destination):
        """

        """
        bonds = {v for bond in self.bonds for v in bond}

        try:
            src, src_port = source
            if isinstance(src, str):
                src = self.components[src]
        except (ValueError, TypeError):
            src, src_port = self._find_port(source)

        try:
            dest, dest_port = destination
            if isinstance(dest, str):
                dest = self.components[dest]

        except (ValueError, TypeError):
            dest, dest_port = self._find_port(destination)

        if ((src, src_port) in bonds and src_port.isnumneric()) or (
                (dest, dest_port) in bonds and dest_port.isnumeric()):
            raise InvalidPortException("Could not join %s to %s: "
                                       "Port already in use",
                                       source, destination)
        bond = (src, src_port), (dest, dest_port)
        src.connect_port(src_port)
        dest.connect_port(dest_port)
        self.bonds.append(bond)

    def _find_port(self, component):

        if isinstance(component, str):
            try:
                comp = self.components[component]
            except AttributeError:
                raise InvalidComponentException(
                    "Could not find %s: not contained in %s", component, self)
        elif component not in self:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component, self)
        else:
            comp = component

        used_ports = {p for bond in self.bonds for (c, p) in bond
                      if c is comp}

        free_ports = set(comp.ports) - used_ports

        if not free_ports:
            try:
                port = comp.make_port()
            except AttributeError:
                raise InvalidComponentException(
                    "Could not find a free port on %s",component)
        elif len(free_ports) > 1:
            raise InvalidComponentException(
                "Could not find a unique free port on %s: "
                "specify a port ", component)
        else:
            port = free_ports.pop()

        return comp, port

    def disconnect(self, component_1, component_2=None,
                   port_1=None, port_2=None):
        """

        Args:
            component_1:
            component_2:
        """
        if component_1 in self.components.keys():
            c1 = self.components[component_1]
        elif component_1 in self.components.values():
            c1 = component_1
        else:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component_1, self
            )
        if port_1 and port_1 not in c1.ports:
            raise InvalidPortException("Could not find port %s on %s",
                                       port_1, component_1)

        if component_2 and component_2 in self.components.keys():
            c2 = self.components[component_2]
        elif component_2 and component_2 in self.components.values():
            c2 = component_2
        elif component_2:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component_2, self
            )
        if component_2 and port_2 and port_2 not in c2.ports:
            raise InvalidPortException("Could not find port %s on %s",
                                       port_2, component_2)

        def cmp(target, test):
            p, q = target
            r, s = test
            if p is r and (not q or (q is s)):
                return True
            else:
                return False

        def is_target(src, dest):
            if not c2:
                return cmp((c1, port_1), src) or cmp((c1, port_1), dest)
            else:
                return (
                   cmp((c1, port_1), src) and cmp((c2, port_2), dest)) \
                       or (cmp((c1, port_1), dest) and cmp((c2, port_2), src))

        target_bonds = [bond for bond in self.bonds if is_target(*bond)]

        for bond in target_bonds:
            self.bonds.remove(bond)
            for c, p in bond:
                c.release_port(p)

