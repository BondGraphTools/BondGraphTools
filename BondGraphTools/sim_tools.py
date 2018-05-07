import numpy as np

from .base import ModelException
from .algebra import inverse_coord_maps
import logging
logger = logging.getLogger(__name__)
j = None
de = None


def start_julia():
    global j
    global de

    logger.info("Starting Julia Interpreter.")
    from diffeqpy import de as de

    import julia

    j = julia.Julia()
    logger.info("Julia Interpreter Loaded.")


def simulate(system,
             timespan, initial, control_vars=None):
    """
    Simulate the system dynamics.

    Args:
        system obj:`BondGraph`:
        timespan tuple(float):
        initial list(float):
        control_vars (str,list(str), dict(str)):

    Returns:

    """
    if system.ports:
        raise ModelException(
            "Cannot Simulate %s: unconnected ports %s",
            system, system.ports)

    if system.control_vars and not input:
        raise ModelException("Control variable not specified")

    func, diffs = _build_dae(system, control_vars)
    dx0, x0 = initial
    DX0 = np.array(dx0, dtype=np.float64)
    X0 = np.array(x0, dtype=np.float64)
    tspan = tuple(float(t) for t in timespan)
    if not de:
        start_julia()
    problem = de.DAEProblem(func, DX0, X0, tspan, differential_vars=diffs)
    sol = de.solve(problem)

    return np.transpose(sol.t), np.transpose(sol.u)


def _build_dae(system, control_vars=None):

    mappings, coords = inverse_coord_maps(*system.basis_vectors)
    ss_map, js_map, cv_map = mappings

    m = len(ss_map)
    k = len(cv_map)
    if len(js_map) > 0:
        raise NotImplementedError("Bond Graph has unconnected Ports")

    # construct julia coords

    julia_string = "function f(dX, X, p, t)\n"
    derivatives = set(coords[0:m])
    differential_vars = []
    end_string = "    return ["
    for i, relation in enumerate(system.constitutive_relations):
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
                raise NotImplementedError("No control var for %s",
                                          control_vars)

        for k, idx in ss_map.items():
            temp_string = temp_string.replace(f"x_{idx}", f"X[{idx+1}]")

        julia_string += f"    res{i+1} = {temp_string}\n"
        if i > 0:
            end_string += ', '

        end_string += f"res{i+1}"

    end_string += "]\nend"
    julia_string += end_string

    if not j: start_julia()

    func = j.eval(julia_string)

    return func, differential_vars


