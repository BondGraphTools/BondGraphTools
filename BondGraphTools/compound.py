
import sympy as sp

from .base import BondGraphBase
from .exceptions import *
from .view import GraphLayout
from .algebra import smith_normal_form, adjacency_to_dict, \
    inverse_coord_maps, _handle_constraints


class BondGraph(BondGraphBase):
    def __init__(self, name, components=None, **kwargs):

        super().__init__(name, **kwargs)
        self.components = dict()

        if components:
            for component in components:
                self.__add__(component)

        self.bonds = list()
        self._internal_ports = []
        self.type = "Composite"
        self.view = GraphLayout(self)
        """Graphical Layout of internal components"""

        self.cv_map = dict()
        self._port_map = dict()

    def save(self, filename):
        raise NotImplementedError

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        if self.__dict__ != other.__dict__:
            return False

        for c1, c2 in zip(self.components.values(),
                          other.componets.values()):
            if c1 != c2:
                return False

        return True

    def __getitem__(self, item):
        return self.components[item]

    def __contains__(self, item):
        if isinstance(item, str):
            return item in self.components
        elif isinstance(item, BondGraphBase):
            for comp in self.components.values():
                if item is comp:
                    return True

        return False

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __add__(self, other):

        if isinstance(other, BondGraphBase):
            n = 0
            for k in self.components.keys():
                t_idx, n_idx = k.split('_')

                n_idx = int(n_idx)
                if t_idx == other.type and n_idx >= n:
                    n = n_idx + 1

            self.components[f"{other.type}_{n}"] = other

        else:
            raise TypeError("unsupported operand type(s) for +: '%s' and '%s'",
                            type(self),
                            type(other))

        return self
    # def __iadd__(self, other):
    #     return self.__add__(other)

    @property
    def params(self):
        return {
            (c,p): v.params[p] for c,v in self.components.items() if v.params
            for p in v.params if not isinstance(p, (float,complex, int))
        }

    @property
    def ports(self):
        bonds = {v for bond in self.bonds for v in bond}
        j = 0
        out = {p:v for p, v in self._ports.items()}

        for v in self.components.values():
            for p in v.ports:
                while j in out:
                    j += 1
                if (v, p) not in bonds:
                    out.update({j: (v, p)})
                    j += 1
        return out

    @property
    def state_vars(self):
        j = 0
        out = dict()
        for v in self.components.values():
            for i in v.state_vars:
                out.update({f"x_{j}": (v, i)})
                j += 1
        return out

    @property
    def control_vars(self):
        j = 0
        out = dict()
        for v in self.components.values():
            for i in v.control_vars:
                out.update({f"u_{j}":(v,i)})
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

        coordinates, mappings, lin_op, nlin_op = self._system_rep()

        inv_tm, inv_js, _ = mappings
        js_size = len(inv_js)  # number of ports
        ss_size = len(inv_tm)  # number of state space coords

        if nlin_op:
            Rx = (sp.eye(lin_op.rows) - lin_op).dot(coordinates)
            Rx = [r - nlin_op[i] for i,r in enumerate(Rx)]
            subs = [pair for pair in zip(coordinates, Rx)]
            nlin_op = nlin_op.subs(subs)

        relations = []
        for row in range(ss_size):
            rel = lin_op[row, :].dot(coordinates)

            if nlin_op:
                rel += nlin_op[row]
            if rel:
                relations.append(rel)

        for row in range(ss_size, ss_size + 2*js_size, 2):
            if lin_op[row, row] == 0 or lin_op[row+1, row+1] == 0:
                rels = lin_op[row:row+2, :].dot(coordinates)
                relations += [rel for rel in rels if rel]

        rels = lin_op[(ss_size + 2 * js_size):, :].dot(coordinates)
        relations += [rel for rel in rels if rel]

        return relations

    def _system_rep(self):
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
        nlin_vect = []
        try:
            lin_row = max(row + 1 for row, _ in lin_dict.keys())
        except ValueError:
            lin_row = 0

        nlin_row = 0

        for component in self.components.values():
            relations = component.get_relations_iterator(mappings, coordinates)
            for linear, nonlinear in relations:
                if not nonlinear:
                    lin_dict.update({(lin_row, k): v
                                     for k, v in linear.items()})
                    lin_row += 1
                else:

                    nlin_dict.update({(nlin_row, c):v for c,v in linear.items()})
                    nlin_vect.append(nonlinear)
                    nlin_row += 1

        if nlin_dict:
            lin_dict.update({
                (lin_row + row, col): value
                for (row, col), value in nlin_dict.items()})

            nonlinear_op = sp.zeros(lin_row, 1).col_join(
                sp.Matrix(nlin_row, 1, nlin_vect)
            )

            linear_op, nonlinear_op = smith_normal_form(
                sp.SparseMatrix(lin_row + nlin_row, n, lin_dict), nonlinear_op)

        else:

            linear_op = smith_normal_form(
                sp.SparseMatrix(lin_row + nlin_row, n, lin_dict)
            )

            nonlinear_op = None

        # if nlin_dict:
        #     nonlinear_op = _build_nonlinear_operator((nlin_dict, nlin_vect),
        #                                              coordinates)
        # else:
        #     nonlinear_op = None

        # if we have a constraint like x_0 + x_1 - u_0 = 0
        # turn it into dx_0 + dx_1 - du_0 = 0

        coordinates, linear_op, nonlinear_op = _handle_constraints(
            linear_op, nonlinear_op, coordinates, size_tuple)

        return coordinates, mappings, linear_op, nonlinear_op

    def _build_internal_basis_vectors(self):
        tangent_space = dict()
        control_space = dict()
        port_space = dict()
        input_port_space = {
            sp.symbols((f"e_{i}", f"f_{i}")): (self, i)
            for i in self._internal_ports
        }
        for component in self.components.values():
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

    def _validate_port(self, target, port=None, **kwargs):

        if isinstance(target, tuple):
            return self._validate_port(*target)

        # target is either str or component
        if isinstance(target, str):
            comp = self.components[target]
            port = port
        elif isinstance(target, int) and target in self._internal_ports:
            return self, target
        elif target is self and port in self.ports:
            return target, port
        elif isinstance(target, BondGraphBase) and \
            target in self.components.values():
            comp = target
        else:
            raise InvalidComponentException(
                "Could not find component %s", target
            )

        if port in comp.ports:
            return comp, port
        elif port and port not in comp.ports:
            return comp, comp.make_port(port=port)
        else:
            return self._find_port(comp)

    def connect(self, source, destination):
        """

        """
        bonds = {v for bond in self.bonds for v in bond}

        src, src_port = self._validate_port(source)

        dest, dest_port = self._validate_port(destination)

        if ((src, src_port) in bonds) or ((dest, dest_port) in bonds):
            raise InvalidPortException("Could not join %s to %s: "
                                       "Port already in use",
                                       source, destination)

        bond = (src, src_port), (dest, dest_port)

        src.connect_port(src_port)
        dest.connect_port(dest_port)
        self.bonds.append(bond)

    def _find_port(self, component):

        if isinstance(component, str):
            try:
                comp = self.components[component]
            except AttributeError:
                raise InvalidComponentException(
                    "Could not find %s: not contained in %s", component, self)
        elif component not in self:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component, self)
        else:
            comp = component

        used_ports = {p for bond in self.bonds for (c, p) in bond
                      if c is comp}

        free_ports = set(comp.ports) - used_ports

        if not free_ports:
            try:
                port = comp.make_port()
            except AttributeError:
                raise InvalidPortException(
                    "Could not find a free port on %s",component)
        elif len(free_ports) > 1:
            raise InvalidPortException(
                "Could not find a unique free port on %s: "
                "specify a port ", component)
        else:
            port = free_ports.pop()

        return comp, port

    def disconnect(self, component_1, component_2=None,
                   port_1=None, port_2=None):
        """

        Args:
            component_1:
            component_2:
        """
        if component_1 in self.components.keys():
            c1 = self.components[component_1]
        elif component_1 in self.components.values():
            c1 = component_1
        else:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component_1, self
            )
        if port_1 and port_1 not in c1.ports:
            raise InvalidPortException("Could not find port %s on %s",
                                       port_1, component_1)

        if component_2 and component_2 in self.components.keys():
            c2 = self.components[component_2]
        elif component_2 and component_2 in self.components.values():
            c2 = component_2
        elif component_2:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component_2, self
            )
        if component_2 and port_2 and port_2 not in c2.ports:
            raise InvalidPortException("Could not find port %s on %s",
                                       port_2, component_2)

        def cmp(target, test):
            p, q = target
            r, s = test
            if p is r and (not q or (q is s)):
                return True
            else:
                return False

        def is_target(src, dest):
            if not c2:
                return cmp((c1, port_1), src) or cmp((c1, port_1), dest)
            else:
                return (
                   cmp((c1, port_1), src) and cmp((c2, port_2), dest)) \
                       or (cmp((c1, port_1), dest) and cmp((c2, port_2), src))

        target_bonds = [bond for bond in self.bonds if is_target(*bond)]

        for bond in target_bonds:
            self.bonds.remove(bond)
            for c, p in bond:
                c.release_port(p)

    def make_port(self, port=None):
        if port and not isinstance(port, int):
            raise InvalidPortException("Could not make port %s", port)

        if port and port not in self.ports:
            n = port
        else:
            n = 0
            while n in self.ports:
                n += 1

        self._ports[n] = (self, n)
        self._internal_ports.append(n)
        return n

    def delete_port(self, port):
        if port in self._ports:
            del self._ports[port]
            del self._internal_ports[port]




