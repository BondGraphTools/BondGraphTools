from collections import OrderedDict
from .base import BondGraphBase, InvalidPortException, \
    ModelParsingError
from .compound import BondGraph
from .view import Glyph


class BaseComponent(BondGraphBase):
    """
    Atomic bond graph components are those defined by constitutive relations.
    """

    def __init__(self, type, constitutive_relations,
                 state_vars=None, params=None, **kwargs):
        super().__init__(**kwargs)

        self._state_vars = state_vars

        self._params = params

        self.view = Glyph()
        self.type = type
        self._constitutive_relations = constitutive_relations

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    @property
    def control_vars(self):
        if self.params:
            return [param for param, value in self.params.items() if not value]
        else:
            return []

    @property
    def params(self):
        return self._params

    @property
    def state_vars(self):
        return self._state_vars

    @property
    def model(self):

        vects = self._generate_basis_vectors()

        models = []

        for var in self.state_vars:
            var_type, _ = var.split("_")

            if var_type == "q":
                ef_var = f'f_{port}'
            elif var_type == "p":
                ef_var = f'e_{port}'
            else:
                raise ModelParsingError(
                    "Error parsing model {}: "
                    "state variable {} must be either p or q",
                    self.type, var
                )
            models.append({vects[ef_var]: -1, vects[f"d{var}"]: 1})

        for relation in self._constitutive_relations:
            pass

    def _generate_basis_vectors(self):

        vects = OrderedDict()

        for var in self.state_vars:
            vects[var] = (self, var)

        for port in self.ports:
            if not port.isnumeric():
                continue
            vects[f"e_{port}"] = ((self, port), 'e')
            vects[f"f_{port}"] = ((self, port), 'f')

        for var in self.state_vars:
            vects[f"d{var}"] = (self, f"d{var}")

        for control in self.control_vars:
            vects[control] = (self, control)

        return vects

    def __add__(self, other):
        return BondGraph(
            name="{}+{}".format(self.name, other.name),
            components=[self, other]
        )

    def __hash__(self):
        return id(self)

