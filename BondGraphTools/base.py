"""Base classes for bond graph models and connections."""

import logging
from collections import namedtuple

logger = logging.getLogger(__name__)


class BondGraphBase(object):
    """
    Base class definition for all bond graphs.

    Args:
        name:               Name of this model, assumed to be unique.
        parent:             Parent model, or none.
        metamodel:          Metamodel class.

    """

    def __init__(self,
                 name=None,
                 parent=None,
                 metamodel=None,
                 **kwargs):
        # TODO: This is a dirty hack
        # Job for meta classes maybe?
        if not metamodel:
            self.__metamodel = "BG"
        else:
            self.__metamodel = metamodel

        if not name:
            self.name = f"{self.metamodel}" \
                        f"{self.__class__.instances}"
        else:
            self.name = name

        self.parent = parent
        """Model that contains this object."""
        self.view = None

    def __repr__(self):
        return f"{self.metamodel}: {self.name}"

    def __new__(cls, *args, **kwargs):  # noqa: D102
        if "instances" not in cls.__dict__:
            cls.instances = 1
        else:
            cls.instances += 1

        return object.__new__(cls)

    def __del__(self):
        self.instances -= 1

    @property
    def metamodel(self):
        """The meta-model type of this object."""  # noqa: D401
        return self.__metamodel

    @property
    def template(self):
        """The model template from which this was created."""  # noqa: D401
        raise NotImplementedError

    @property
    def constitutive_relations(self):
        """The equations governing this objects dynamics."""   # noqa: D401
        raise NotImplementedError

    @property
    def uri(self):
        """Model reference locator."""  # noqa: D401
        if not self.parent:
            return f"{self.name}:"
        else:
            return f"{self.parent.uri}/{self.name}"

    @property
    def root(self):
        """The root of the tree to which this object belongs."""  # noqa: D401
        if not self.parent:
            return self
        else:
            return self.parent.root

    @property
    def basis_vectors(self):
        """The input/ouput, dynamic and interconnect vectors."""  # noqa: D401
        raise NotImplementedError

    def __hash__(self):
        return id(self)


class Bond(namedtuple("Bond", ["tail", "head"])):
    """Stores the connection between two ports.

    Head and tail are specified to determine orientation

    Attributes:
        head: The 'harpoon' end of the power bond and direction of positive $f$
        tail: The non-harpoon end, and direction of negative $f$

    """

    def __contains__(self, item):
        if isinstance(item, BondGraphBase):
            return self.head[0] is item or self.tail[0] is item

        try:
            c, i = item
            return any(c == comp and i == idx for comp, idx in self)
        except TypeError:
            return False


class Port(object):
    """Basic object for representing ports.

    Looks and behaves like a namedtuple:
        component, index =  Port

    Args:
        component:      The owner of this port.
        index:          The index of this port in the component.

    """

    def __init__(self, component, index):
        self.component = component
        """(`PortManager`) The component that this port is attached to"""
        self.index = index
        """(int) The numerical index of this port"""
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

    def __contains__(self, item):
        return item is self.component

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
            return ((self.component is other.component) and # noqa
                    (self.index == other.index))
        except AttributeError:
            try:
                c, p = other
                return c is self.component and p == self.index
            except (AttributeError, TypeError):
                pass
        return False
