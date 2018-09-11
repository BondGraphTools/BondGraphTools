import logging
import sympy as sp

from .base import BondGraphBase, Bond, Port
from .exceptions import *
from .view import GraphLayout
from .algebra import smith_normal_form, adjacency_to_dict, \
    inverse_coord_maps, reduce_model, get_relations_iterator

logger = logging.getLogger(__name__)

class BondGraph(BondGraphBase):
    def __init__(self, name, components=None, **kwargs):

        super().__init__(name, **kwargs)
        self.components = set()

        if components:
            for component in components:
                self.add(component)

        self._bonds = BondSet()

        self._internal_ports = []
        self.view = GraphLayout(self)
        """Graphical Layout of internal components"""
        self.cv_map = dict()
        self._port_map = dict()
        self._model_changed = True

    @property
    def bonds(self):
        return list(self._bonds)

    @property
    def metaclass(self):
        return "BondGraph"

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

    def add(self, *args):

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
        if [b for b in self._bonds if b.head.component is component or
                b.tail.component is component]:
            raise InvalidComponentException("Component is still connected")
        if component not in self.components:
            raise InvalidComponentException("Component not found")

        component.parent = None
        self.components.remove(component)

    @property
    def params(self):
        j = 0
        out = dict()
        for v in self.components:
            try:
                params = v.params
            except AttributeError:
                continue
            for p in params:
                out.update({j: (v, p)})
                j+=1
        return out


    # def set_param(self, param, value):
    #     c, v = param
    #     self.components[c].set_param(v, value)

    @property
    def ports(self):


        # bonds = {v for bond in self.bonds for v in bond}
        # j = 0
        # out = {p:v for p, v in self._ports.items()}
        #
        # for v in self.components:
        #     for p in v.ports:
        #         while j in out:
        #             j += 1
        #         if (v, p) not in bonds:
        #             out.update({j: (v, p)})
        #             j += 1
        # print(len(self._bonds))

        return {}

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
        for v in self.components:
            try:
                for i in v.control_vars:
                    out.update({f"u_{j}":(v,i)})
                    j += 1
            except AttributeError:
                pass
        return out

    @property
    def basis_vectors(self):
        tangent_space = dict()
        port_space = dict()
        control_space = dict()

        for var, var_id in self.state_vars.items():
            tangent_space[sp.symbols((f"{var}", f"d{var}"))] = var_id

        for port, port_id in self.ports.items():
            port_space[sp.symbols((f"e_{port}", f"f_{port}"))] = port_id

        for var, var_id in self.control_vars.items():
            control_space[sp.symbols(f"{var}")] = var_id

        return tangent_space, port_space, control_space

    @property
    def constitutive_relations(self):

        coordinates, mappings, lin_op, nlin_op, constraints = self.system_model()

        inv_tm, inv_js, _ = mappings
        js_size = len(inv_js)  # number of ports
        ss_size = len(inv_tm)  # number of state space coords

        # Rx = (sp.eye(lin_op.rows) - lin_op)
        # Rxs = [r - nlin_op[i] for i, r in enumerate(Rx.dot(coordinates))]
        # subs = [pair for i, pair in enumerate(zip(coordinates, Rxs)) if
        #         Rx[i,i] == 0]
        # logger.info("Substituting with subs: %s", subs)
        # nlin_op = nlin_op.subs(subs)
        #
        # for k,c in enumerate(coordinates):
        #     if Rx[k,k] == 0:
        #         subs.append(
        #             (c, c - nlin_op[k])
        #
        #         )
        coord_vect = sp.Matrix(coordinates)
        relations = [
            sp.Add(l,r) for i, (l,r) in enumerate(zip(
                lin_op*coord_vect,nlin_op))
            if not ss_size <= i < ss_size + 2*js_size
        ]
        if isinstance(constraints, list):
            for constraint in constraints:
                logger.info("Adding constrain %s", repr(constraint))
                if constraint:
                    relations.append(constraint)
        else:
            logger.warning("Constraints %s is not a list. Discarding",
                           repr(constraints))

        return [r.nsimplify() for r in relations if r]

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

        for component in self.components:
            relations = get_relations_iterator(component, mappings, coordinates)
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
        port_space = dict()
        input_port_space = {
            sp.symbols((f"e_{i}", f"f_{i}")): (self, i)
            for i in self._internal_ports
        }
        for component in self.components:
            c_ts, c_ps, c_cs = component.basis_vectors

            for var_id in c_ts.values():
                i = len(tangent_space)
                tangent_space[sp.symbols((f"x_{i}", f"dx_{i}"))] = var_id

            for port in c_ps.values():
                i = len(port_space) + len(input_port_space)
                port_space[sp.symbols((f"e_{i}", f"f_{i}"))] = port

            for cv in c_cs.values():
                i = len(control_space)
                control_space[sp.symbols(f"u_{i}")] = cv

        port_space.update(input_port_space)
        return tangent_space, port_space, control_space


    # def make_port(self, port=None):
    #     if port and not isinstance(port, int):
    #         raise InvalidPortException("Could not make port %s", port)
    #
    #     if port and port not in self.ports:
    #         n = port
    #     else:
    #         n = 0
    #         while n in self.ports:
    #             n += 1
    #
    #     self._ports[n] = (self, n)
    #     self._internal_ports.append(n)
    #     return n

    # def delete_port(self, port):
    #     if port in self._ports:
    #         del self._ports[port]
    #         del self._internal_ports[port]
    #
    # def find(self, name, component=None):
    #     """
    #     Searches through this model for a component of with the specified name.
    #     If the type of component is specified by setting c_type, then only
    #     components of that type are considered.
    #
    #     Args:
    #         name (str): The name of the object to search for.
    #         component (str): (Optional) the class of components in which to
    #          search.
    #
    #     Returns:
    #         None if no such component exists, or the component object if it
    #      does.
    #
    #     Raises:
    #     """
    #     out = [obj for obj in self.components if
    #            (not component or obj.type == component) and
    #            obj.name == name]
    #     if len(out) > 1:
    #         raise NotImplementedError("Object is not unique")
    #     elif len(out) == 1:
    #         return out.pop()
    #     else:
    #         return None


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

