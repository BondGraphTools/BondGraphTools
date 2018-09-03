import logging
import sympy as sp

from .base import BondGraphBase
from .exceptions import *
from .compound import BondGraph
from .view import Glyph

logger = logging.getLogger(__name__)


class BaseComponent(BondGraphBase):
    """
    Atomic bond graph components are those defined by constitutive relations.
    """

    def __init__(self, type, constitutive_relations,
                 state_vars=None, params=None, **kwargs):
        super().__init__(**kwargs)

        self._state_vars = state_vars

        self._params = params

        self.view = Glyph(self)
        self.type = type
        self._constitutive_relations = constitutive_relations

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    @property
    def control_vars(self):
        if self.params:
            return [param for param, value in self.params.items()
                    if ((not value) and not isinstance(value,(int,float, complex)))
                    or (isinstance(value, dict) and "value"
                    not in value)]
        else:
            return []

    @property
    def params(self):
        return self._params if self._params else {}

    def set_param(self, param, value):
        if isinstance(self._params[param], dict):
            self._params[param]["value"] = value
        else:
            self._params[param] = value

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
            components=[self, other])

    def __hash__(self):
        return super().__hash__()

    def make_port(self, **kwargs):
        raise AttributeError(
            "Cannot add a port to %s", self)


class NPort(BaseComponent):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self._fixed_ports = set(int(p) for p in self.ports if isinstance(p, int))

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    def __add__(self, other):
        return BondGraph(
            name="{}+{}".format(self.name, other.name),
            components=[self, other]
        )

    @property
    def ports(self):
        return {p: v for p, v in self._ports.items() if isinstance(p, int)}

    def make_port(self, port=None, value=None):

        if port and not isinstance(port, int):
            raise InvalidPortException("Could not make port %s", port)

        if port and port not in self.ports:
            n = port
        else:
            n = 0
            while n in self.ports:
                n += 1

        self._ports[n] = value

        return n

    def delete_port(self, port):
        del self._ports[port]

    def release_port(self, port):
        if port not in self._fixed_ports:
            self.delete_port(port)


class NPortWeighted(NPort):
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    def __add__(self, other):
        return BondGraph(
            name="{}+{}".format(self.name, other.name),
            components=[self, other]
        )

    def make_port(self, port=None, value=1):
        return super().make_port(port=port, value=value)

    def release_port(self, port):
        if port not in self._fixed_ports:
            self.delete_port(port)

            #del self._params[f"c_{port}"]

    def _build_relations(self):

        rels = super()._build_relations()
        subs = []

        for port in self.ports:
            if port in self._fixed_ports:
                pair = (f"c_{port}", self.ports[port]["value"])
            else:
                pair = (f"c_{port}", self.ports[port])

            subs.append(pair)
        return [rel.subs(subs) for rel in rels]