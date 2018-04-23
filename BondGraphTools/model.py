"""
Bond Graph Model base files.
"""

import logging
import sympy
import copy

from .view import GraphLayout, Glyph
from .component_manager import get_component, base_id

logger = logging.getLogger(__name__)


def new(component, name=None, library=base_id, value=None):
    """
    Creates a new Bond Graph from a library component.

    Args:
        component(str or obj): The type of component to create.
         If a string is specified, the the component will be created from the
         appropriate libaray. If an existing bond graph is given, the bond
         graph will be cloned.
        name (str): The name for the new component
        library (str): The library from which to find this component (if
        component is specified by string).

    Returns: instance of :obj:`BondGraph`

    """
    if isinstance(component, str):
        build_args = get_component(component, library)

        if name:
            build_args.update({"name": name})
        if value:
            _update_build_params(build_args, value)

        obj = AtomicComponent(type=component, **build_args)

    elif isinstance(component, BondGraphBase):
        obj = copy.copy(component)
        if name:
            obj.name = name
        if value:
            _update_build_params(obj.__dict__, value)

    return obj


def _update_build_params(build_args, value):
    if isinstance(value, (list, tuple)):
        assignments = zip(value, build_args["params"].keys())
        for param, v in assignments:
            build_args["params"][param]["value"] = v
    elif isinstance(value, dict):
        for param, v in value:
            build_args["params"][param]["value"] = v
    else:
        p = next(iter(build_args["params"]))
        build_args["params"][p]["value"] = value


class BondGraphBase:
    def __init__(self, name, parent=None,
                 ports=None, params=None):
        """
        Base class definition for all bond graphs.

        Args:
            name: Assumed to be unique
            metadata (dict):
        """
        #super().__init__(name, parent)
        self.name = name
        self._ports = ports
        """ List of exposed Power ports"""
        self.control_vars = list()
        """ List of exposed control variables """
        self.params = params
        """ Dictionary of internal parameter and their values. The key is 
        the internal parameter, the value may be an exposed control value,
        a function of time, or a constant."""
        self.view = None

        if params:
            for param, value in self.params.items():
                if not value:
                    self.control_vars.append(param)

    @property
    def ports(self):
        return self._ports

    @property
    def state_vars(self):
        return NotImplementedError

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class AtomicComponent(BondGraphBase):
    """
    Atomic bond graph components are those defined by constitutive relations.
    """
    def __init__(self, type, constitutive_relations,
                 state_vars=None, **kwargs):
        super().__init__(**kwargs)

        self._state_vars = state_vars
        """ List of state variables"""

        self.view = Glyph()
        self.type = type
        self.constitutive_relations = constitutive_relations

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    @property
    def state_vars(self):
        return self._state_vars

    @property
    def model(self):
        model = []
        for rel_str in self.constitutive_relations:
            model.append(sympy.sympify(rel_str))

    def __add__(self, other):
        return BondGraph(
            name="{}+{}".format(self.name, other.name),
            components=[self, other]
        )

    def __hash__(self):
        return id(self)


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

        self.control_vars = []

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
    def ports(self):
        bonds = {v for bond in self.bonds for v in bond}

        return {(c, p): v.ports[p]
                for c, v in self.components.items() for p in v.ports
                if (v, p) not in bonds}

    @property
    def state_vars(self):
        return {(c, i): v.state_vars[i] for
                c, v in self.components.items() if v.state_vars
                for i in v.state_vars}

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
            raise InvalidPortException("Could not join {} to {}: "
                                       "Port already in use",
                                       source, destination)
        bond = (src, src_port), (dest, dest_port)
        self.bonds.append(bond)

    def _find_port(self, component):

        if isinstance(component, str):
            try:
                comp = self.components[component]
            except AttributeError:
                raise InvalidComponentException(
                    "Could not find {}: not contained in {}", component, self)
        elif component not in self:
            raise InvalidComponentException(
                "Could not find {}: not contained in {}", component, self)
        else:
            comp = component

        used_ports = {p for bond in self.bonds for (c, p) in bond
                      if c is comp and p.isnumeric()}

        free_ports = set(comp.ports) - used_ports
        if not free_ports:
            raise InvalidComponentException(
                "Could not find a free port on {}",component)
        elif len(free_ports) > 1:
            raise InvalidComponentException(
                "Could not find a unique free port on {}: "
                "specify a port ", component)

        return comp, free_ports.pop()

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
                "Could not find {}: not contained in {}", component_1, self
            )
        if port_1 and port_1 not in c1.ports:
            raise InvalidPortException("Could not find port {} on {}",
                                       port_1, component_1)

        if component_2 and component_2 in self.components.keys():
            c2 = self.components[component_2]
        elif component_2 and component_2 in self.components.values():
            c2 = component_2
        elif component_2:
            raise InvalidComponentException(
                "Could not find {}: not contained in {}", component_2, self
            )
        if component_2 and port_2 and port_2 not in c2.ports:
            raise InvalidPortException("Could not find port {} on {}",
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

        self.bonds = [bond for bond in self.bonds if not is_target(*bond)]


class InvalidPortException(Exception):
    pass


class InvalidComponentException(Exception):
    pass



