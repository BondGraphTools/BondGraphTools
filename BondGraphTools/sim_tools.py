import numpy as np
from scipy.integrate import ode
import sympy as sp
from .base import ModelException
from .algebra import inverse_coord_maps,smith_normal_form

def simulate(bond_graph,
             timespan, initial_state,input=None, delta_t=0.001,
             dtype=np.float32):

    if bond_graph.ports:
        raise ModelException(
            "Cannot Simulate %s: unconnected ports %s",
            bond_graph, bond_graph.ports)

    if bond_graph.control_vars and not input:
        raise ModelException("Control variable not specified")

    t = np.linspace(*timespan, 1/delta_t)

    x = np.empty(shape=(len(initial_state), len(t)), dtype=dtype)
    for i, xi_0 in enumerate(initial_state):
        x[i,0] = xi_0

    func, jacobian = _build_ode(bond_graph)

    return t, x


def _build_ode(bond_graph, input=None):

    mappings, coords = inverse_coord_maps(*bond_graph.basis_vectors)
    ss_map, js_map, cv_map = mappings

    m = len(ss_map)
    if len(js_map) > 0:
        raise NotImplementedError

    lin_dict = {}
    row = 0
    for lin, nlin in bond_graph.get_relations_iterator(mappings, coords):

        if nlin:
            raise NotImplementedError

        for col, value in lin.items():
            lin_dict.update({(row, col): value})

        row += 1
    if row != m:
        raise NotImplementedError



    return None, None

