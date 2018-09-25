import logging

import numpy as np
import sympy as sp
from sympy.core import SympifyError

from .config import config
from .exceptions import ModelException

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
    x_subs = []
    dx_subs = []

    for i, x in enumerate(model.state_vars):
        x_subs.append((x, X[i+1]))
        dx_subs.append((sp.S(f'dx_{i}'), dX[i+1]))


    cv_strings, dcv_strings = generate_control_strings(
        list(model.control_vars.keys()),
        control_vars,
        x_subs,
        dx_subs
    )

    differential_vars = []
    subs = x_subs + dx_subs

    function_header = ""
    function_body = ""
    function_footer = ""

    if in_place:
        function_header = "function f(res, dX, X, p, t)\n"
    else:
        k = len(model.constitutive_relations)
        function_header = "function f(dX, X, p, t)\n"
        function_header += f"    res = zeros({k})\n"

    for cv, dcv in zip(cv_strings, dcv_strings):
        function_header += cv
        if dcv:
            function_header += dcv

    assert model.constitutive_relations

    for relation in model.constitutive_relations:

        eqn_str = str(relation.subs(subs))
        eqn_str = eqn_str.replace('**', '^')
        if 'dX' in eqn_str:
            differential_vars.append(True)
        else:
            differential_vars.append(False)

        function_body += f"    res[{len(differential_vars)}] = {eqn_str}\n"

    if in_place:
        function_footer += "end\n"
    else:
        function_footer += "    return res\n end\n"

    out_str = function_header + function_body + function_footer

    return (out_str,
            differential_vars)


def generate_control_strings(cv, cv_substitutions, x_subs, dx_subs):

    if isinstance(cv_substitutions, dict):
        pairs = [(cv_i, cv_substitutions[cv_i]) for cv_i in cv]
    elif isinstance(cv_substitutions, list):
        pairs = list(zip(cv, cv_substitutions))
    elif len(cv) == 1 and (isinstance(cv_substitutions, (str, sp.Expr))):
        pairs = [cv[0], cv_substitutions]
    elif not cv and not cv_substitutions:
        return [], []
    else:
        raise NotImplementedError

    cv_strings = []
    dcv_strings = []
    subs = x_subs + dx_subs
    partial_pairs = [(X_i,DX_i) for (_,X_i), (_,DX_i) in zip(x_subs, dx_subs)]
    for var, val in pairs:
        try:

            u_i = sp.sympify(val).subs(subs)
            du_i = u_i.diff(sp.S('t')) + sum([
                u_i.diff(X)*DX for X,DX in partial_pairs
            ])
            cv_strings.append(f"    {str(var)} = {u_i}\n")
            dcv_strings.append(f"    d{str(var)} = {du_i}\n")
        except SympifyError as ex:
            s = val

            for x, X in reversed(x_subs):
                s = s.replace(str(x), str(X))

            cv_strings.append(f"    {str(var)} = {s}\n")
            dcv_strings.append(None)

    return cv_strings, dcv_strings


class SolverException(Exception):
    pass