import numpy as np
import sympy as sp
from .config import config
import os, sys

from .exceptions import ModelException
from .algebra import inverse_coord_maps, create_ds
import logging
logger = logging.getLogger(__name__)


def simulate(system,
             timespan,
             x0,
             dx0=None,
             dt=0.1,
             control_vars=None):
    """
    Simulate the system dynamics.

    Args:
        system obj:`BondGraph`:
        timespan tuple(float):
        initial list(float):
        control_vars (str,list(str), dict(str)):

    Returns: t, X

    """

    de = config.de
    j = config.julia
    if system.ports:
        raise ModelException(
            "Cannot Simulate %s: unconnected ports %s",
            system, system.ports)

    if system.control_vars and not control_vars:
        raise ModelException("Control variable not specified")

    tspan = tuple(float(t) for t in timespan)
    X0 = np.array(x0, dtype=np.float64)
    assert len(X0) == len(system.state_vars)

    func_str, diffs = to_julia_function_string(system, control_vars)

    if dx0:
        DX0 = np.array(dx0, dtype=np.float64)
    else:
        DX0 = np.zeros(X0.shape, dtype=np.float64)

    func = j.eval(func_str)
    problem = de.DAEProblem(func, DX0, X0, tspan, differential_vars=diffs)

    sol = de.solve(problem, dense=True, saveat=float(dt))

    if sol.retcode not in ("Default", "Success"):
        raise SolverException("Integration error: Solver returned %s "
                              % sol.retcode, sol)

    t = np.transpose(sol.t)

    return np.resize(t, (len(t), 1)), np.transpose(sol.u).T

class Simulation(object):
    def __init__(self, model,
                 timespan=None,
                 x0=None,
                 dx_0=None,
                 control_vars=None):

        self.sol = None




    def run(self, x0, timespan):
        pass


def to_julia_function_string(model, control_vars=None, in_place=False):
    """
    Produces a Julia function string from the given model.

    We expect that that control_vars is a dict with the same keys,
    or list of the same size, as the model.control_vars

    Args:
        model:
        control_vars:

    Returns:
        (string, list)
        A string containing the function definition, and a list of bools
        identifing which variables contain derivatives.
    """

    dX = sp.IndexedBase('dX')
    X = sp.IndexedBase('X')
    subs = []

    iterable = {}
    if isinstance(control_vars, list) \
            and len(control_vars) == len(model.control_vars):

        iterable = {k:v for k,v in zip(list(model.control_vars.keys()),
                                        control_vars)}
    elif isinstance(control_vars, dict):
        iterable = {}
        for u in model.control_vars:
            if u in control_vars:
                iterable[u] = control_vars[u]
            else:
                iterable[u] = control_vars[str(u)]

    for i, x in enumerate(model.state_vars):
        subs.append((x, X[i+1]))
        subs.append((sp.S(f'dx_{i}'), dX[i+1]))

    for i, (cv, val) in enumerate(iterable.items()):

        u_t = sp.sympify(val).subs(subs)
        try:
            du_t = u_t.diff('t') + sum([
                dX[i+1]*(u_t.diff(X[i+1]).subs(subs))
                for i in range(len(model.state_vars))
            ])
        except:
            raise NotImplementedError
        subs.append((cv, u_t))

        # check if val is a number, if so, directly substitute.

        # check if val is a differential function of space and time

        # check if val is a graph


    differential_vars = []
    if in_place:
        out_str = "function f(res, dX, X, p, t)\n"
    else:
        k = len(model.constitutive_relations)
        out_str = "function f(dX, X, p, t)\n"

        out_str += f"    res = zeros({k})\n"
    for relation in model.constitutive_relations:

        eqn_str = str(relation.subs(subs))
        eqn_str = eqn_str.replace('**', '^')
        if 'dX' in eqn_str:
            differential_vars.append(True)
        else:
            differential_vars.append(False)

        out_str += f"    res[{len(differential_vars)}] = {eqn_str}\n"

    if in_place:
        out_str += "end\n"
    else:
        out_str+="    return res\n end\n"

    return (out_str,
            differential_vars)


class SolverException(Exception):
    pass