import logging
import sympy as sp

from .base import *
from .exceptions import *
from .view import GraphLayout
from .algebra import smith_normal_form, adjacency_to_dict, \
    inverse_coord_maps, reduce_model, get_relations_iterator

logger = logging.getLogger(__name__)

class BondGraph(BondGraphBase, LabeledPortManager):
    def __init__(self, name, components=None, **kwargs):

        BondGraphBase.__init__(self, name, **kwargs)
        LabeledPortManager.__init__(self)
        self.components = set()

        if components:
            for component in components:
                self.add(component)

        self._bonds = BondSet()

        self.view = GraphLayout(self)
        """Graphical Layout of internal components"""

        self._port_map = dict()
        self._model_changed = True

    @property
    def template(self):
        return None

    @property
    def bonds(self):
        return list(self._bonds)

    def __truediv__(self, other):
        try:
            if not self.parent:
                test_uri = f"/{other}"
            else:
                test_uri = f"{self.uri}/{other}"
            c, = (c for c in self.components if c.uri == test_uri)
            return c
        except (ValueError, TypeError):

            raise ValueError(f"Cannot find {other}")

    @property
    def metaclass(self):
        return "BG"

    @bonds.setter
    def bonds(self, arg):
        raise AttributeError("Use add/remove functions.")

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        if self.__dict__ != other.__dict__:
            return False

        for c1, c2 in zip(self.components,
                          other.componets):
            if c1 != c2:
                return False

        return True

    @property
    def internal_ports(self):
        return [p for c in self.components for p in c.ports]

    def map_port(self, label, e, f):
        """

        Args:
            label:
            e:
            f:

        Returns:

        """

        port = self.get_port(label)
        self._port_map[port] = (e,f)

    def add(self, *args):
        """
        Warning: Scheduled to be deprecated
        """

        def validate(component):
            if not isinstance(component, BondGraphBase):
                raise InvalidComponentException("Invalid component class")
            if component is self:
                raise InvalidComponentException("Cannot add a model to itself")
            elif component.root is self.root:
                raise InvalidComponentException(
                    "Component already exists in model")

        work_list = []
        for arg in args:
            if isinstance(arg, BondGraphBase):
                validate(arg)
                work_list.append(arg)
            elif isinstance(arg, list):
                for item in arg:
                    validate(item)
                    work_list.append(item)
            else:
                raise InvalidComponentException(f"Invalid Component: {arg}")

        for item in work_list:
            item.parent = self
            self.components.add(item)


    def remove(self, component):
        """
        Warning: Scheduled to be deprecated
        """
        if [b for b in self._bonds if b.head.component is component or
                b.tail.component is component]:
            raise InvalidComponentException("Component is still connected")
        if component not in self.components:
            raise InvalidComponentException("Component not found")

        component.parent = None
        self.components.remove(component)

    def set_param(self, param, value):
        """
        Warning: Scheduled to be deprecated
        """
        c, p = self.params[param]
        c.set_param(p, value)

    @property
    def params(self):
        j = 0
        out = dict()

        excluded = {
            v for pair in self._port_map.values() for v in pair
        }

        for v in self.components:
            try:
                params = v.params
            except AttributeError:
                continue
            for p in params:
                param = (v, p)
                if param not in excluded:
                    out.update({j: param})
                    j+=1
        return out

    @property
    def state_vars(self):
        j = 0
        out = dict()
        for v in self.components:
            try:
                x_local = v.state_vars
            except AttributeError:
                continue

            for i in x_local:
                out.update({f"x_{j}": (v, i)})
                j += 1

        return out

    @property
    def control_vars(self):
        j = 0
        out = dict()
        excluded = {
            v for pair in self._port_map.values() for v in pair
         }

        for v in self.components:
            try:
                for i in v.control_vars:
                    cv = (v,i)
                    if cv not in excluded:
                        out.update({f"u_{j}": cv})
                        j += 1
            except AttributeError:
                pass
        return out

    @property
    def basis_vectors(self):
        """
        Basis vectors for the state space (X), port space (J),
        and control space (U) from an external point of view.

        For the state space dictionaries are of the form
        ```
            X = {
                sympy.Symbol('x_i'): (object, var)
            }
        ```
        We assume the object is a subclass of BondGraphBase
        and the var refers to the variable name in the objects local
        co-ordinate system and may be a string or a sympy.Symbol

        For the port space, dictionaries are of the form
        ```
            J = {
                (sympy.Symbol(e_i), sympy.Symbol(f_i)): Port(obj, idx)
            }
        ```
        where Port is an instance of `Port`.

        Finally for the cotrol variables we have
        ```
            U = {
            sympy.Symbol(u_i):(object, var)
            }
        ```
        Where object and var are specified as per the state space.
        """

        tangent_space = dict()
        control_space = dict()

        for var, var_id in self.state_vars.items():
            tangent_space[sp.symbols((f"{var}", f"d{var}"))] = var_id

        port_space = self._port_vectors()

        for var, var_id in self.control_vars.items():
            control_space[sp.symbols(f"{var}")] = var_id

        return tangent_space, port_space, control_space

    @property
    def constitutive_relations(self):
        coordinates, mappings, lin_op, nlin_op, constraints = self.system_model()
        inv_tm, inv_js, _ = mappings
        out_ports = [idx for p, idx in inv_js.items() if p in self.ports]
        logger.debug("Getting IO ports: %s",out_ports)
        js_size = len(inv_js)  # number of ports
        ss_size = len(inv_tm)  # number of state space coords

        coord_vect = sp.Matrix(coordinates)
        relations = [
            sp.Add(l,r) for i, (l,r) in enumerate(zip(
                lin_op*coord_vect,nlin_op))
            if not ss_size <= i < ss_size + 2*js_size - 2*len(out_ports)
        ]
        if isinstance(constraints, list):
            for constraint in constraints:
                logger.debug("Adding constraint %s", repr(constraint))
                if constraint:
                    relations.append(constraint)
        else:
            logger.warning("Constraints %s is not a list. Discarding",
                           repr(constraints))
        subs = []

        for local_idx, c_idx in enumerate(out_ports):
            p, = {pp for pp in self.ports if pp.index == local_idx}
            label = p.index
            subs.append(sp.symbols((f"e_{c_idx}", f"e_{label}")))
            subs.append(sp.symbols((f"f_{c_idx}", f"f_{label}")))

        return [r.subs(subs).simplify().nsimplify() for r in relations if r]

    def system_model(self, control_vars=None):
        mappings, coordinates = inverse_coord_maps(
            *self._build_internal_basis_vectors()
        )
        inv_tm, inv_js, inv_cv = mappings

        js_size = len(inv_js) # number of ports
        ss_size = len(inv_tm) # number of state space coords
        cv_size = len(inv_cv)
        n = len(coordinates)

        size_tuple = (ss_size, js_size, cv_size, n)

        lin_dict = adjacency_to_dict(inv_js, self.bonds, offset=ss_size)

        nlin_dict = {}

        try:
            row = max(row + 1 for row, _ in lin_dict.keys())
        except ValueError:
            row = 0

        inverse_port_map = {}

        for port, (cv_e, cv_f) in self._port_map.items():
            inverse_port_map[cv_e] = ss_size + 2*inv_js[port]
            inverse_port_map[cv_f] = ss_size + 2*inv_js[port] + 1

        for component in self.components:
            relations = get_relations_iterator(
                component, mappings, coordinates, inverse_port_map
            )

            for linear, nonlinear in relations:
                lin_dict.update({(row, k): v
                                 for k, v in linear.items()})
                nlin_dict.update({(row, 0): nonlinear})
                row += 1

        linear_op = sp.SparseMatrix(row, n, lin_dict)
        nonlinear_op = sp.SparseMatrix(row, 1, nlin_dict)
        coordinates, linear_op, nonlinear_op, constraints = reduce_model(
                linear_op, nonlinear_op, coordinates, size_tuple,
            control_vars=control_vars
        )

        return coordinates, mappings, linear_op, nonlinear_op, constraints

    def _build_internal_basis_vectors(self):
        tangent_space = dict()
        control_space = dict()
        port_space = {}
        # bond_space = dict()
        #
        # for i, bond in enumerate(self.bonds):
        #     bond_space[sp.symbols((f"e_{i}", f"f_{i}"))] = bond

        mapped_cvs = {
            var for pair in self._port_map.values() for var in pair
        }

        for component in self.components:
            c_ts, c_ps, c_cs = component.basis_vectors

            for var_id in c_ts.values():
                i = len(tangent_space)
                tangent_space[sp.symbols((f"x_{i}", f"dx_{i}"))] = var_id

            for cv in c_cs.values():
                if cv not in mapped_cvs:
                    i = len(control_space)
                    control_space[sp.symbols(f"u_{i}")] = cv

            for port in c_ps.values():
                i = len(port_space)
                port_space[sp.symbols((f"e_{i}", f"f_{i}"))] = port

        n = len(port_space)
        external_ports = {
            sp.symbols((f"e_{n + i}", f"f_{n + i}")): port
            for i, port in enumerate(self._port_map)
        }
        port_space.update(external_ports)

        return tangent_space, port_space, control_space



def _is_label_invalid(label):
    if not isinstance(label, str):
        return True

    for token in [" ", ".", "/"]:
        if len(label.split(token)) >1:
            return True

    return False

class BondSet(set):
    """
    Container class for internal bonds.
    """
    def add(self, bond):
        tail = bond.tail
        head = bond.head
        super().add(bond)
        head.is_connected = True
        tail.is_connected = True

    def remove(self, bond):
        tail = bond.tail
        head = bond.head
        if bond in self:
            super().remove(bond)
        else:
            super().remove(Bond(head, tail))
        head.is_connected = False
        tail.is_connected = False

    def __contains__(self, item):
        p1, p2 = item
        return super().__contains__(item) or super().__contains__((p2,p1))
