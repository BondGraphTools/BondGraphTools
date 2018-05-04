import numpy as np


from diffeqpy import de
import julia
from .base import ModelException
from .algebra import inverse_coord_maps,smith_normal_form

j = julia.Julia()

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

    func = _build_dae(bond_graph)

    return t, x


def _build_dae(bond_graph, control_vars=None):

    mappings, coords = inverse_coord_maps(*bond_graph.basis_vectors)
    ss_map, js_map, cv_map = mappings

    m = len(ss_map)
    k = len(cv_map)
    if len(js_map) > 0:
        raise NotImplementedError

    # construct julia coords

    julia_string = "function f(dX, X, p, t)\n"
    derivatives = set(coords[0:m])
    differential_vars = []
    end_string = "    return ["
    for i, relation in enumerate(bond_graph.constitutive_relations):
        differential_vars.append(
            derivatives ^ relation.atoms() != set()
        )
        temp_string = str(relation)
        for k, v in cv_map.items():
            if isinstance(control_vars, dict):
                temp_string = temp_string.replace(f"u_{v}", control_vars[k])
            elif isinstance(control_vars, list):
                temp_string = temp_string.replace(f"u_{v}", control_vars[v])
            elif len(cv_map) == 1 and isinstance(control_vars, str):
                temp_string = temp_string.replace(f"u_{v}", control_vars)
            else:
                raise NotImplementedError

        for k, idx in ss_map.items():
            temp_string = temp_string.replace(f"x_{idx}", f"X[{idx+1}]")

        julia_string += f"    res{i+1} = {temp_string}\n"
        if i > 0:
            end_string += ', '

        end_string += f"res{i+1}"

    end_string += "]\nend"
    julia_string += end_string
    func = j.eval(julia_string)

    return func, differential_vars


