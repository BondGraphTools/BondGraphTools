from collections import OrderedDict
import sympy as sp

from .base import BondGraphBase, ModelParsingError, InvalidPortException
from .compound import BondGraph
from .view import Glyph
from .algebra import extract_coefficients


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
        return self._params if self._params else {}

    @property
    def state_vars(self):
        return self._state_vars if self._state_vars else []

    @property
    def constitutive_relations(self):
        """

        Returns:

        """
        models = self._build_relations()

        for var in self.state_vars:
            var_type, port = var.split("_")

            if var_type == "q":
                ef_var = f'f_{port}'
            elif var_type == "p":
                ef_var = f'e_{port}'
            else:
                raise ModelParsingError(
                    "Error parsing model %s: "
                    "state variable %s must be either p or q",
                    self.type, var
                )

            models.append(sp.sympify(f"d{var} - {ef_var}"))

        subs = []
        for param, value in self.params.items():
            if isinstance(value, (int, float, complex)):
                subs.append((sp.symbols(param), value))
            elif isinstance(value, str) and value.isnumeric():
                subs.append((sp.symbols(param), float(value)))
            elif isinstance(value, dict) and "value" in value:
                if isinstance(value["value"], (int, float, complex)):
                    subs.append((sp.symbols(param), value["value"]))
                elif isinstance(value["value"], str) \
                        and value["value"].isnumeric():
                    subs.append((sp.symbols(param), float(value["value"])))

        return [model.subs(subs) for model in models]

        # for each relation, pull out the linear part

    @property
    def basis_vectors(self):

        tangent_space = dict()
        port_space = dict()
        control_space = dict()

        for var in self.state_vars:
            tangent_space[sp.symbols((var, f"d{var}"))] = (self, var)

        for port in self.ports:
            if not isinstance(port, int):
                continue

            port_space[sp.symbols((f"e_{port}", f"f_{port}"))] = (self, port)

        for control in self.control_vars:
            control_space[sp.symbols(control)] = (self, control)

        return tangent_space, port_space, control_space

    def _build_relations(self):
        rels = []
        for string in self._constitutive_relations:
            iloc = 0
            iloc = string.find("_i", iloc)

            if iloc < 0:
                # just a plain old string; sympy can take care of it
                rels.append(sp.sympify(string))
                continue

            sloc = string.rfind("sum(", 0, iloc)

            if sloc < 0:
                # we have a vector equation here.
                for port_id in self.ports:
                    if isinstance(port_id, int):
                        rels.append(
                        sp.sympify(
                            string.replace("_i", "_{}".format(port_id))))

            else:
                tiers = 0

                next_open = string.find("(", sloc + 4)
                eloc = string.find(")", sloc + 4)
                while next_open > 0 or tiers > 0:
                    if next_open < eloc:
                        tiers += 1
                        next_open = string.find("(", next_open)
                    else:
                        tiers -= 1
                        eloc = string.find(")", eloc)

                if eloc < 0:
                    raise ValueError("Unbalanced brackets", string)

                substring = string[sloc + 4: eloc]
                terms = [substring.replace("_i", "_{}".format(p))
                         for p in self.ports if isinstance(p, int)]
                symstr = string[0:sloc] + "(" + " + ".join(terms) + string[
                                                                    eloc:]
                rels.append(
                    sp.sympify(symstr)
                )

        return [r for r in rels if r != 0]

    def __add__(self, other):
        return BondGraph(
            name="{}+{}".format(self.name, other.name),
            components=[self, other]
        )

    def __hash__(self):
        return id(self)

    def make_port(self):
        raise AttributeError(
            "Cannot add a port to %s", self)


class NPort(BaseComponent):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.__fixed_ports = (int(p) for p in self.ports if isinstance(p, int))

    @property
    def ports(self):
        return {p:v for p,v in self._ports.items() if isinstance(p, int)}

    def make_port(self):
        n = 0
        while n in self._ports:
            n += 1
        self._ports[n] = None

        return n

    def delete_port(self, port):
        del self._ports[port]

    def release_port(self, port):
        if port not in self.__fixed_ports:
            self.delete_port(port)
