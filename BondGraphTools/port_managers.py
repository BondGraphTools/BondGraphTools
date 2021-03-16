"""port_managers.py

This module contains the interface for creating and connecting up
models.
"""


import logging

import sympy as sp

from BondGraphTools.exceptions import InvalidPortException
from BondGraphTools.base import Port

logger = logging.getLogger(__name__)


class PortManager(object):
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
            port: (optional) The index or reference to the requested port.

        Returns: An instance of `Port`

        Raises: InvalidPortException
        """
        # If no port is specified, and there is only one port, grab it.
        target_port = None

        if port is None and len(self._ports) == 1:
            p = list(self._ports.keys())[0]
            target_port = p
        # If it's a port object, then grab it
        elif port in self._ports:
            target_port = port
        elif isinstance(port, int):
            target_port = next((pp for pp in self._ports if pp.index == port),
                               None)

        if target_port and not target_port.is_connected:
            return target_port

        elif target_port and target_port.is_connected:
            raise InvalidPortException(f"Port is already connected: "
                                       f"{self}.{port}")
        else:
            raise InvalidPortException(f"Could not find port: {self}.{port}")

    def _port_vectors(self):
        return {
            sp.symbols((f"e_{port.index}", f"f_{port.index}")): port
            for port in self._ports
        }

    @property
    def ports(self):
        """A dictionary of the active ports"""
        return self._ports


class PortExpander(PortManager):
    """
    Class for handling templated virtual ports; for example summing junctions.
    Additionally, attributes can be attached to the generated ports by via a
    dictionary of attributes.

    Args:
        ports (dict): The different port classes (keys) and corresponding
         dictionary of attributes and values (maybe be None)

    Keyword Args:
        static_ports: Ports to be created as per `PortManager`

    See Also: `PortTemplate`
    """

    def __init__(self, ports, static_ports=None):
        if static_ports:
            super().__init__(static_ports)
        else:
            super().__init__({})

        self._templates = {PortTemplate(self, p, v) for p, v in ports.items()}

        if len(self._templates) == 1:
            self._default_template, = self._templates
        else:
            self._default_template = False
        self.max_index = len(static_ports) if static_ports else 0

    def get_port(self, port=None):
        logger.info("Looking for port %s", port)
        if port is not None:
            this_port = None

            if port in self._templates:
                logger.info("Port template found")
                this_port = next((p for p in port.ports if not p.is_connected),
                                 None)
                logger.info("Corresponding port found: %s", this_port)

            if isinstance(port, (str, int)):
                this_port = next(
                    (p for p in self._templates if p.index == port), None
                )

            if this_port is None:
                logger.info("Raising exception..")
                raise InvalidPortException("No free ports of type %s", port)
            else:
                return this_port

        return super().get_port(port)

    def new_port(self, port=None):
        logger.info("trying to create a new port")
        if port is None:
            return self._spawn(self._default_template)
        elif isinstance(port, int):
            return self._spawn(self._default_template, port)
        elif port in self._templates:
            return self._spawn(port)

        raise InvalidPortException(f"Could not create new port:{port}")

    def _spawn(self, port_template, index=None):

        if index is None:
            index = self.max_index
        port = port_template.spawn(index)
        self._ports[port] = port_template.index
        self.max_index = max(index, self.max_index) + 1
        return port

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

    def __init__(self, parent, index, data=None):

        self.parent = parent
        self.index = index
        self.ports = []
        self.data = data if data else {}

    def __hash__(self):
        return id(self)

    def spawn(self, index=None):
        """
        Generates a new port from this template.

        Args:
            index: The desired index of this port.

        Returns:
            `Port`

        Raises: InvalidPortException
        """
        port = ExpandedPort(self.parent, index, port_class=self.index)
        port.__dict__.update({k: v for k, v in self.data.items()})
        self.ports.append(port)

        return port


class LabeledPort(Port):
    """See Also: `Port`, `LabeledPortManager`"""

    def __init__(self, *args, name=None, **kwargs):
        self.name = name
        """The name of this port"""
        Port.__init__(self, *args, **kwargs)

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        if isinstance(other, str) and other == self.name:
            return True
        else:
            return super().__eq__(other)


class LabeledPortManager(PortManager):
    """Interface for labelled ports

    See Also: `PortManager`

    """

    def __init__(self, ports=None):
        if ports:
            super().__init__(ports)
        else:
            super().__init__({})
        self.max_index = len(self._ports)

    def new_port(self, port=None):

        if (isinstance(port, str)
                and port not in {port.name for port in self._ports}):
            port_idx = self.max_index
            port_name = port
            self.max_index += 1
        else:
            port_idx = self.max_index
            port_name = str(port_idx)
            current_names = {port.name for port in self._ports}
            while port_name in current_names:
                port_idx += 1
                port_name = str(port_idx)

        new_port = LabeledPort(self, port_idx, name=port_name)
        self._ports[new_port] = None
        return new_port

    def get_port(self, port=None):
        if isinstance(port, str):
            try:
                p, = {pp for pp in self.ports if pp.name == port}
                return p
            except ValueError:
                raise InvalidPortException("Could not find port %s", port)

        return super().get_port(port)
