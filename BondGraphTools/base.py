"""
Bond Graph Model base files.
"""

import logging

from collections import namedtuple

import sympy as sp
from .exceptions import *

logger = logging.getLogger(__name__)


class BondGraphBase:
    def __init__(self, name=None, parent=None,
                 ports=None, description=None, params=None, metaclass=None):
        """
        Base class definition for all bond graphs.

        Args:
            name: Assumed to be unique
            metadata (dict):
        """

        # TODO: This is a dirty hack
        # Job for meta classes maybe?
        if not metaclass:
            self.__metaclass = "BondGraph"
        else:
            self.__metaclass = metaclass
        if not name:
            self.name = f"{self.metaclass}" \
                        f"{self.__class__.instances}"
        else:
            self.name = name
        self.parent = parent

        self.description = description
        # if ports:
        #     self._ports = {
        #         (int(p) if p.isnumeric() else p):v for p,v in ports.items()
        #     }
        # else:
        #     self._ports = {}
        """ List of exposed Power ports"""

        """ Dictionary of internal parameter and their values. The key is 
        the internal parameter, the value may be an exposed control value,
        a function of time, or a constant."""
        self.view = None

    def __repr__(self):
        return f"{self.metaclass}: {self.name}"

    def __new__(cls, *args, **kwargs):
        if "instances" not in cls.__dict__:
            cls.instances = 1
        else:
            cls.instances += 1

        return object.__new__(cls)

    def __del__(self):
        self.instances -= 1

    @property
    def metaclass(self):
        return self.__metaclass

    @property
    def template(self):
        raise NotImplementedError

    @property
    def constitutive_relations(self):
        raise NotImplementedError

    @property
    def uri(self):
        if not self.parent:
            return "/"
        else:
            parent_uri = self.parent.uri
            if parent_uri == "/":
                return f"{self.parent.uri}{self.name}"
            else:
                return f"{self.parent.uri}/{self.name}"

    @property
    def root(self):
        if not self.parent:
            return self
        else:
            return self.parent.root

    @property
    def basis_vectors(self):
        raise NotImplementedError

    def __hash__(self):
        return id(self)

    # def __eq__(self, other):
    #     return self.__dict__ == other.__dict__


Bond = namedtuple("Bond", ["tail", "head"])

class Port(object):
    """
    Basic object for representing ports;
    Looks and behaves like a namedtuple:
    component, index =  Port
    """

    def __init__(self, component, index):
        self.component = component
        """(`FixedPort`) The component that this port is attached to"""
        self.index = index
        """(int) The numberical index of this port"""
        self.is_connected = False
        """(bool) True if this port is plugged in."""

    def __iter__(self):
        return iter((self.component, self.index))

    def __len__(self):
        return 2

    def __getitem__(self, item):
        if item == 0:
            return self.component
        elif item == 1:
            return self.index
        else:
            raise KeyError

    def __hash__(self):
        return id(self)

    def __str__(self):
        if len(self.component.ports) > 1:
            return f"{self.component.name}.{self.index}"
        else:
            return f"{self.component.name}"

    def __repr__(self):
        return f"Port({self.component}, {self.index})"

    def __eq__(self, other):
        try:
            return ((self.component is other.component) and
                    (self.index == other.index))
        except AttributeError:
            try:
                c, p = other
                return c is self.component and p == self.index
            except AttributeError:
                pass
        return False

class FixedPort:
    """
    This class provides methods for interfacing with static ports on
    components.

    Args:
        ports (dict): The enumerated list of ports associated with
         this component.
    """

    def __init__(self, ports):
        self._ports = {}

        for port, data in ports.items():
            port_data = data if data else {}
            self._ports.update({Port(self, int(port)): port_data})

    def get_port(self, port=None):
        """
        Makes available a (or the) port for use.

        Args:
            port: (optional) The index or referece to the requested port.

        Returns: An instance of `Port`

        Raises: InvalidPortException
        """

        # If no port is specified, and there is only one port, grab it.
        if not port and not isinstance(port, int) and len(self._ports) == 1:
            for p in self._ports:
                if not p.is_connected:
                    return p
        # If it's a port object, then grab it
        elif port in self._ports and not port.is_connected:
            return port
        elif isinstance(port, int):
            p, = (pp for pp in self._ports if pp.index == port and
                  not pp.is_connected)
            if p:
                return p

        if port:
            raise InvalidPortException(f"Could not find port: {self}.{port}")
        else:
            raise InvalidPortException(f"Could not find a free port: {self}")

    def _port_vectors(self):
        return {
            sp.symbols((f"e_{port.index}", f"f_{port.index}")): port
            for port in self._ports
        }

    @property
    def ports(self):
        """A dictionary of the active ports"""
        return self._ports

class PortExpander(FixedPort):
    """
    Class for handling templated virtual ports; for example summing junctions.
    Additionally, attributes can be attached to the generated ports by via a
    dictionary of attributes.

    Args:
        ports (dict): The different port classes (keys) and corresponding
         dictionary of attributes and values (maybe be None)

    Keyword Args:
        static_ports: Ports to be created as per `FixedPort`

    See Also: `PortTemplate`
    """
    def __init__(self, ports, static_ports=None):
        if static_ports:
            super().__init__(static_ports)
        else:
            super().__init__({})

        self._templates = {PortTemplate(self, p, v) for p,v in ports.items()}

        if len(self._templates) == 1:
            self._default_template, = self._templates
        else:
            self._default_template = False
        self.max_index = len(static_ports) if static_ports else 0

    def get_port(self, port=None):
        """
        Tries to fetch the specified port, or tries to generate the port from
        a template.

        Args:
            port: The Port, port index, template, or template string for
             requested port

        Returns: An instance of `Port`.
        """
        if not port and not isinstance(port, int):
            # port is None, so lets try and make a new one
            try:
                return self._default_template.spawn()
            except AttributeError:
                raise InvalidPortException("You must specify a port")
        elif port in self._templates:
            # this is a template, so lets
            template, = {t for t in self._templates if t == port}
            return template.spawn()

        elif isinstance(port, str):
            template, = {t for t in self._templates if t.index == port}
            return template.spawn()

        elif isinstance(port, int) and port\
                not in [p.index for p in self._ports]:
            try:
                return self._default_template.spawn(port)
            except AttributeError:
                raise InvalidPortException("You must specify a port")
        try:    # suppose we've got port or a port tuple that exists
            return super().get_port(port)
        except InvalidPortException as ex:
            pass

        raise InvalidPortException(f"Could not create new port:{port}")

    def _port_vectors(self):
        return {
            sp.symbols((f"e_{port.index}", f"f_{port.index}")): port
            for port in self._ports if port.is_connected
        }

class ExpandedPort(Port):
    def __init__(self, *args, port_class):
        super().__init__(*args)
        self.port_class = port_class

    def __str__(self):
        if self.port_class:
            return f"{self.component.name}.{self.port_class}"
        else:
            return f"{self.component.name}"


class PortTemplate(object):
    """
    Template class for generating new ports of a specific type.

    Args:
          parent (`PortExpander`): The associated object which is to hold
           the port instances.
          index (srt): The label or identifier

    """
    def __new__(cls, parent, index, data=None):
        self = object.__new__(cls)
        self.parent = parent
        self.index = index
        self.ports = []
        self.data = data if data else {}
        return self

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        try:
            p, s = other
            if p == self.parent and s == self.index:
                return True
        except TypeError:
            if other is self:
                return True
        return False

    def spawn(self, index=None):
        """
        Generates a new port from this template.

        Args:
            index: The desired index of this port.

        Returns:
            `Port`

        Raises: InvalidPortException
        """
        if not index:
            index = self.parent.max_index
        elif index in [p.index for p in self.parent.ports]:
            raise InvalidPortException("Could not create port: index "
                                       "already exists")
        port = ExpandedPort(self.parent, index, port_class=self.index)
        port.__dict__.update({k:v for k,v in self.data.items()})
        self.parent._ports[port] = self.index
        self.ports.append(port)
        self.parent.max_index = max(index, self.parent.max_index) + 1
        return port

class LabeledPort(Port):
    def __init__(self,*args,name=None, **kwargs):
        self.name = name
        Port.__init__(self, *args, **kwargs)

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        if isinstance(other, str) and other == self.name:
            return True
        else:
            return super().__eq__(self, other)

class LabeledPortManager(FixedPort):

    def __init__(self, ports=None):
        if ports:
            super().__init__(ports)
        else:
            super().__init__({})
        self.max_index = len(self._ports)

    def get_port(self, port=None):
        if isinstance(port, str):
            try:
                p, = {pp for pp in self.ports if pp.name == port}
                return p
            except ValueError:
                idx = self.max_index
                self.max_index = + 1
                new_port = LabeledPort(self, idx, name=port)
                self._ports[new_port] = None
                return new_port
        elif not port and len(self._ports) == 0:
            idx = self.max_index
            self.max_index += 1
            new_port = LabeledPort(self, idx, name=str(idx))
            self._ports[new_port] = None
            return new_port
        else:
            return super().get_port(port)
