import logging
import sympy as sp

from .base import *
from .exceptions import *
from .compound import BondGraph
from .view import Glyph

logger = logging.getLogger(__name__)

class BaseComponent(BondGraphBase, FixedPort):
    """
    Atomic bond graph components are those defined by constitutive relations.
    """

    def __init__(self, metaclass, constitutive_relations,
                 state_vars=None, params=None, **kwargs):

        self._metaclass = metaclass
        ports = kwargs.pop("ports")
        super().__init__(**kwargs)
        FixedPort.__init__(self, ports)
        self._state_vars = state_vars

        self._params = params

        self.view = Glyph(self)
        self._constitutive_relations = constitutive_relations

    def __eq__(self, other):
        return self.__dict__ == self.__dict__

    @property
    def template(self):
        return f"{self.__library__}/{self.__component__}"

    @property
    def metaclass(self):
        return self._metaclass

    @property
    def control_vars(self):
        """See `BondGraphBase`"""

        def is_const(value):
            if isinstance(value, (int, float, complex)):
                return True
            elif isinstance(value, sp.Symbol):
                return True
            else:
                return False

        out = []

        for p, v in self.params.items():
            try:
                if is_const(v) or is_const(v["value"]):
                    continue
            except (KeyError, TypeError):
                pass

            out.append(p)
        return out

    # @property
    # def max_ports(self):
    #     return len(self._ports)

    @property
    def params(self):
        """See `BondGraphBase`"""
        return self._params if self._params else {}

    def set_param(self, param, value):
        """
        Warning: Scheduled to be deprecated
        """
        if isinstance(self._params[param], dict):
            self._params[param]["value"] = value
        else:
            self._params[param] = value

    @property
    def state_vars(self):
        """See `BondGraphBase`"""
        return self._state_vars if self._state_vars else []

    @property
    def constitutive_relations(self):
        """See `BondGraphBase`"""
        models = self._build_relations()
        # for var in self.state_vars:
        #     var_type, port = var.split("_")
        #
        #     if var_type == "q":
        #         ef_var = f'f_{port}'
        #     elif var_type == "p":
        #         ef_var = f'e_{port}'
        #     else:
        #         raise ModelParsingError(
        #             "Error parsing model %s: "
        #             "state variable %s must be either p or q",
        #             self.metaclass, var
        #         )
        #
        #     models.append(sp.sympify(f"d{var} - {ef_var}"))

        subs = []

        def _value_of(v):
            if isinstance(v, (int, float, complex, sp.Symbol)):
                return v
            elif not v:
                raise KeyError
            elif isinstance(v, str):
                v_out, = v.split(" ")
                return sp.S(v_out)
            elif isinstance(v, dict):
                return _value_of(v["value"])
            else:
                raise ValueError(f"Invalid Parameter")

        for param, value in self.params.items():
            try:
                v = _value_of(value)
                subs.append((sp.symbols(param), v))
            except KeyError:
                pass
            except ValueError as ex:
                raise ValueError(f"{self}.{oaram}: {ex.args}")

        return [model.subs(subs) for model in models]

        # for each relation, pull out the linear part

    @property
    def basis_vectors(self):
        """See `BondGraphBase.basis_vectors`"""
        port_space = self._port_vectors()

        tangent_space = dict()

        control_space = dict()

        for var in self.state_vars:
            tangent_space[sp.symbols((var, f"d{var}"))] = (self, var)

        # for port in self.ports:
        #     if not isinstance(port, int):
        #         continue
        #
        #     port_space[sp.symbols((f"e_{port}", f"f_{port}"))] = (self, port)

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

    def __hash__(self):
        return super().__hash__()


class BaseComponentSymmetric(BaseComponent):
    """
    Refer to `BaseComponent`.

    Instances of this class are multi-port components which are able to have
    connections made without specifying ports.
    """
    def get_port(self, port=None):
        """See `BaseComponent`"""
        if not port and not isinstance(port, int):
            p = [port for port in self.ports if not port.is_connected]

            if not p:
                raise InvalidPortException("No free ports")
            return super().get_port(p[0])

        else:
            return super().get_port(port)


class EqualEffort(BondGraphBase, PortExpander):
    """Implements 0-junction.

    Attributes:
         template:
         view:
         basis_vectors:
         constitutive_relations:

    """
    def __init__(self, **kwargs):

        PortExpander .__init__(self, {None: None})
        BondGraphBase.__init__(self, **kwargs)
        self.view = Glyph(self)

    @property
    def template(self):
        return "base/0"

    @property
    def basis_vectors(self):
        return {}, self._port_vectors(), {}

    @property
    def constitutive_relations(self):

        vars = list(self._port_vectors())
        e_0, f_0 = vars.pop()
        partial_sum = f_0

        relations = []

        while vars:
            e_i, f_i = vars.pop()
            relations.append(
                e_i - e_0
            )
            partial_sum += f_i

        relations.append(partial_sum)

        return relations

class EqualFlow(BondGraphBase, PortExpander):

    def __init__(self, **kwargs):
        PortExpander.__init__(self, {"non_inverting": {"weight": 1},
                                     "inverting": {"weight": -1}})
        BondGraphBase.__init__(self, **kwargs)
        self.view = Glyph(self)

    @property
    def non_inverting(self):
        t, = (tp for tp in self._templates if tp.index == "non_inverting")
        return t

    @property
    def inverting(self):
        t, = (tp for tp in self._templates if tp.index == "inverting")
        return t

    @property
    def template(self):
        return "base/1"

    @property
    def basis_vectors(self):
        return {},  self._port_vectors(), {}

    def get_port(self, port=None):
        try:
            return super().get_port(port)
        except InvalidPortException as ex:
            if not port:
                raise InvalidPortException("You must specify a port")

    @property
    def constitutive_relations(self):

        relations = []

        var = list(self._port_vectors().items())
        (e_0, f_0), port = var.pop()

        sigma_0 = port.weight
        partial_sum = sigma_0*e_0

        while var:
            (e_i, f_i), port = var.pop()
            sigma_i = port.weight
            partial_sum += sigma_i*e_i
            relations.append(sigma_i*f_i - sigma_0*f_0)

        relations.append(partial_sum)
        return relations
