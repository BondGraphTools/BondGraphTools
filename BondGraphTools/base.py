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
        if ports:
            self._ports = {
                (int(p) if p.isnumeric() else p):v for p,v in ports.items()
            }
        else:
            self._ports = {}
        """ List of exposed Power ports"""

        """ Dictionary of internal parameter and their values. The key is 
        the internal parameter, the value may be an exposed control value,
        a function of time, or a constant."""
        self.view = None

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

    # @property
    # def max_ports(self):
    #     raise NotImplementedError

    @property
    def constitutive_relations(self):
        raise NotImplementedError

    @property
    def uri(self):
        if not self.parent:
            return ""
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


Port = namedtuple("Port", ["component", "index"])
Bond = namedtuple("Bond", ["tail", "head"])


class FixedPort:
    """
    Args:
        ports
    """

    def __init__(self, ports):
        self._ports = {}

        for port, data in ports.items():
            port_data = data if data else {}
            port_data.update({"is_connected":False})
            self._ports.update({Port(self, int(port)): port_data})

    def get_port(self, port=None):
        if not port and not isinstance(port, int) and len(self._ports) == 1:
            for p in self._ports:
                if not self._ports[p]["is_connected"]:
                    return p

        elif port in self._ports and not self._ports[port]["is_connected"]:
            return port

        elif isinstance(port, int):
            p, = (pp for pp, data in self._ports.items() if pp.index == port and
                    not data["is_connected"])
            if p:
               return p

        raise InvalidPortException

    def _pre_connect_hook(self, port):
        pass

    def _post_connect_hook(self, port):
        self._ports[port]["is_connected"] = True
        pass

    def _pre_disconnect_hook(self, port):
        pass

    def _post_disconnect_hook(self, port):
        self._ports[port]["is_connected"] = False


    def _port_vectors(self):
        return {
            sp.symbols((f"e_{port.index}", f"f_{port.index}")): port
            for port in self._ports
        }

    @property
    def ports(self):
        return self._ports

class PortExpander(FixedPort):

    def __init__(self,ports, static_ports=None):
        if static_ports:
            super().__init__(static_ports)
        else:
            super().__init__({})

        self._port_sets = {p: [] for p in ports}
        self.__max_idx = len(static_ports) if static_ports else 0

    def get_port(self, port=None):

        try:
            return super().get_port(port)
        except InvalidPortException:
            pass

        if not port and len(self._port_sets) == 1:
            return self.get_port(port=next(iter(self._port_sets)))
        elif port in self._port_sets:
            return self._new_port(port)
        else:
            raise InvalidPortException("Could not get a port")


    def _new_port(self, port):
        idx = self.__max_idx
        self.__max_idx += 1
        new_port = Port(self, idx)
        port_data = {"is_connected": False}
        self._port_sets[port].append(idx)
        self._ports.update({new_port:port_data})
        return new_port


    def _port_vectors(self):
        return {
            sp.symbols((f"e_{port.index}", f"f_{port.index}")): port
            for port in self._ports if self._ports[port]["is_connected"]
        }
