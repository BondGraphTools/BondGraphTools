"""Tools for running mdoel simulations"""


import logging

import numpy as np
import sympy as sp
from sympy.core import SympifyError
from scipy.optimize import broyden1

from .exceptions import ModelException, SolverException

logger = logging.getLogger(__name__)


def _fetch_ic(x0, dx0, system, func, eps=0.001):
    if isinstance(x0, list):
        assert len(x0) == len(system.state_vars)
        X0 = np.array(x0, dtype=np.float64)
    elif isinstance(x0, dict):
        X0 = np.array(
            [np.NaN for _ in system.state_vars], dtype=np.float64
        )
        for k, v in x0.items():
            *_ , idx = str(k).split('_')
            idx = int(idx)
            X0[idx] = v
    elif isinstance(x0, (int, float, complex)) and len(system.state_vars) == 1:
        X0 = np.array([x0], dtype=np.float64)
    elif isinstance(x0, np.ndarray) and x0.shape == (len(system.state_vars), ):
        X0 = x0
    else:
        raise ModelException(f"Invalid Initial Conditions: {x0}")

    if dx0:
        DX0 = np.array(dx0, dtype=np.float64)
    else:
        DX0 = np.zeros(X0.shape, dtype=np.float64)

    # if we don't have consistent initial conditions; find them if we can
    # fail if we can't
    f = lambda y: func(y, X0, 0, 0)
    if np.linalg.norm(f(DX0)) > eps:

        DX0 = broyden1(f, DX0)
        if np.linalg.norm(f(DX0)) > 0.001:
            raise ModelException(
                f"Inconsistent initial conditions: "
                f"Could not find dx0 for the given x0 {x0}")

    return X0, DX0

def simulate(system,
             timespan,
             x0,
             dx0=None,
             dt=0.1,
             control_vars=None):
    """Simulate the system dynamics.

    This method integrates the dynamics of the system over the specified
    interval of time, starting at the specified initial state.

    The solver used is a differential-algebraic integrator which respects
    conservation laws and algebraic constraints. It is expected that the
    initial state satisfies the systems inherent algebraic constrains;
    inconsistent initial conditions will raise exceptions.

    The initial values of derivatives can be specified and the solver will
    ensure they are consistent with the initial state, or change them if they
    are not.

    Currently, control variables can take the form of numbers or a strings
    and are assigned via a dictionary or list.

    Permissible strings:

        * numerical constants such as `1.0`, `pi`
        * time `t`
        * state variables; for example `x_0`
        * arithmetic operators such as `+`,`-`, `*`, `/`, as well as `^`
          (power operator), `%` (remainder)
        * elementary math functions such as `sin`, `exp`, `log`
        * ternary if; for example `t < 0 ? 0 : 1` which implements the Heaviside

    Args:
        system :obj:`BondGraph`: The system to simulate
        timespan: A pair (`list` or `tuple`) containing the start and end points
                  of the simulation.
        x0: The initial conditions of the system.
        dx0 (Optional): The initial rates of change of the system. The default
                        value (`None`) indicates that the system should be
                        initialised from the state variable initial conditions.
        dt: The time step between reported (not integrated) values.
        control_vars: A `dict`, `list` or `tuple` specifing the values of the
                      control variables.
    Returns:
        t: numpy array of timesteps
        x: numpy array of state values

    Raises:
        ModelException, SolverException
    """
    from .config import config
    de = config.de
    j = config.julia
    if system.ports:
        raise ModelException(
            "Cannot Simulate %s: unconnected ports %s",
            system, system.ports)

    if system.control_vars and not control_vars:
        raise ModelException("Control variable not specified")

    tspan = tuple(float(t) for t in timespan)


    func_str, diffs = to_julia_function_string(system, control_vars)
    func = j.eval(func_str)
    X0, DX0 = _fetch_ic(x0, dx0, system, func)


    problem = de.DAEProblem(func, DX0, X0, tspan, differential_vars=diffs)

    sol = de.solve(problem, dense=True, saveat=float(dt))

    if sol.retcode not in ("Default", "Success"):
        raise SolverException("Integration error: Solver returned %s "
                              % sol.retcode, sol)

    t = np.transpose(sol.t)

    return np.resize(t, (len(t), 1)), np.transpose(sol.u).T


def to_julia_function_string(model, control_vars=None, in_place=False):
    """
    Produces a Julia function string from the given model.

    We expect that that control_vars is a dict with the same keys,
    or list of the same size, as the model.control_vars

    Args:
        model:
        control_vars:
        in_place:

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


    cv_strings, dcv_strings = _generate_control_strings(
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


def _generate_control_strings(cv, cv_substitutions, x_subs, dx_subs):

    if isinstance(cv_substitutions, dict):
        pairs = [(cv_i, cv_substitutions[cv_i]) for cv_i in cv]
    elif isinstance(cv_substitutions, list):
        pairs = list(zip(cv, cv_substitutions))
    elif len(cv) == 1 and (isinstance(cv_substitutions, (str, sp.Expr))):
        pairs = [cv[0], cv_substitutions]
    elif not cv and not cv_substitutions:
        return [], []
    else:
        raise NotImplementedError(
            f"Could not substitute {cv_substitutions} into {cv}"
        )

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
