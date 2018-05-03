import numpy as np
import functools

from  diffeqpy import de

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

    subs = []
    for c in coords:
        try:
            name, idx = str(c).split('_')
            subs += [(str(c), f'{name}[{idx}]')]
        except ValueError:
            pass

    differential = []
    julia_string  = "def f(dx,x,p,t):\n"
    return_string = "    return ["
    for i, rel in enumerate(bond_graph.constitutive_relation):
        if 'dx' in rel or 'dq' in rel or 'dp' in rel:
            differential.append(True)
        else:
            differential.append(False)

        eqn = functools.reduce(lambda x, y: x.replace(*y), subs, rel)
        julia_string += f'    res{i} = {eqn}\n'
        return_string += 'res{i},'

    julia_string += return_string + ']'
    func = de.eval(julia_string)



    return func, None


