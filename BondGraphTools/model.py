"""
Bond Graph Model base files.
"""

import logging


from .view import GraphLayout, Glyph
from .component_manager import get_component, base_id

logger = logging.getLogger(__name__)


def new(component, name=None, library=base_id, **kwargs):
    """
    Creates a new Bond Graph from a library component.

    Args:
        component:
        name:
        library:
        **kwargs:

    Returns:

    """

    return None


class BondGraph:
    def __init__(self, name, parent=None, metadata=None, **kwargs):
        """
        Base class definition for all bond graphs.

        Args:
            name: Assumed to be unique
            metadata (dict):
        """
        self.name = name

        self.metadata = metadata

        self.ports = None
        """ List of exposed Power ports"""
        self.state_vars = None
        """ List of state variables"""
        self.control_vars = None
        """ List of exposed control variables """
        self.params = None
        """ Dictionary of internal parameter and their values. The key is 
        the internal parameter, the value may be an exposed control value,
        a function of time, or a constant."""

        self.parent = parent
        self.view = None

    @property
    def state_equations(self):
        raise NotImplementedError


class AtomicComponent(BondGraph):
    """
    Atomic bond graph components are those defined by constitutive relations.
    """
    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        self.view = Glyph()
        self.model = kwargs["constitutive_relation"]


class CompositeBondGraph(BondGraph):

    def __init__(self, name, **kwargs):

        super().__init__(name, **kwargs)

        self.components = list()
        """Add and remove components via list"""
        self.bonds = list()
        """Add and remove bonds via list"""
        self.view = GraphLayout()
        """Graphical Layout of internal components"""

        self.control_vars = []

        self.cv_map = dict()

    def save(self, filename):
        raise NotImplementedError

    def __contains__(self, item):
        pass