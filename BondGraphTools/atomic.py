"""This module contains class definitions for atomic components; those which
cannot be decomposed into other components.
"""

import logging
from BondGraphTools.base import BondGraphBase
from BondGraphTools.exceptions import InvalidPortException
from BondGraphTools.view import Glyph
from BondGraphTools.port_managers import PortManager, PortExpander
import sympy as sp

logger = logging.getLogger(__name__)


class Component(BondGraphBase, PortManager):
    """Components defined by constitutive relations."""

    def __init__(self, metamodel, constitutive_relations,
                 state_vars=None, params=None, **kwargs):

        self._metamodel = metamodel
        ports = kwargs.pop("ports")
        super().__init__(**kwargs)
        PortManager.__init__(self, ports)
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
    def metamodel(self):
        return self._metamodel

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
                raise ValueError("Invalid Parameter")

        for param, value in self.params.items():
            try:
                v = _value_of(value)
                subs.append((sp.symbols(param), v))
            except KeyError:
                pass
            except ValueError as ex:
                raise ValueError(f"({self}, {param}): {ex.args}")

        return [model.subs(subs) for model in models]

    @property
    def basis_vectors(self):
        """See `BondGraphBase.basis_vectors`"""
        port_space = self._port_vectors()

        tangent_space = dict()

        control_space = dict()

        for var in self.state_vars:
            tangent_space[sp.symbols((var, f"d{var}"))] = (self, var)

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


class SymmetricComponent(Component):
    """
    Refer to `Component`.

    Instances of this class are multi-port components which are able to have
    connections made without specifying ports.
    """

    def get_port(self, port=None):
        """See `Component`"""
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

        vects = list(self._port_vectors())
        e_0, f_0 = vects.pop()
        partial_sum = f_0

        relations = []

        while vects:
            e_i, f_i = vects.pop()
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
        return {}, self._port_vectors(), {}

    def get_port(self, port=None):
        try:
            return super().get_port(port)
        except InvalidPortException:
            if not port:
                raise InvalidPortException("You must specify a port")

    @property
    def constitutive_relations(self):

        relations = []

        var = list(self._port_vectors().items())
        (e_0, f_0), port = var.pop()

        sigma_0 = port.weight
        partial_sum = sigma_0 * e_0

        while var:
            (e_i, f_i), port = var.pop()
            sigma_i = port.weight
            partial_sum += sigma_i * e_i
            relations.append(sigma_i * f_i - sigma_0 * f_0)

        relations.append(partial_sum)
        return relations
