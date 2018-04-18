"""
Bond Graph Model base files.
"""

import logging


from .view import GraphLayout, Glyph
from .component_manager import get_component, base_id

logger = logging.getLogger(__name__)


def new(component_type, name=None, library=base_id, **kwargs):
    """
    Creates a new Bond Graph from a library component.

    Args:
        component_type:
        name:
        library:
        **kwargs:

    Returns:

    """

    component_data = get_component(component_type, library)

    return AtomicComponent(

    )


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
        self.state_vars = None
        self.params = None
        self.observables = None
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
        self.mapping = list()


class CompositeBondGraph(BondGraph):

    def __init__(self, name, **kwargs):

        super().__init__(name, **kwargs)

        self.components = list()
        """Add and remove components via list"""
        self.bonds = list()
        """Add and remove bonds via list"""
        self.view = GraphLayout()
        """Graphical Layout of internal components"""

        self.parameter_map = dict()
        """Mapping between class level parameter values and component values"""

    def save(self, filename):
        raise NotImplementedError
