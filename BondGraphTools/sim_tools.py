"""Tools for running mdoel simulations"""


import logging

import numpy as np
import sympy as sp
from sympy.core import SympifyError
from scipy.optimize import broyden1
from scikits.odes.dae import dae
from .exceptions import ModelException, SolverException

logger = logging.getLogger(__name__)


def _fetch_ic(x0, dx0, system, func, t0, eps=0.001):
    if isinstance(x0, list):
        assert len(x0) == len(system.state_vars)
        X0 = np.array(x0, dtype=np.float64)
    elif isinstance(x0, dict):
        X0 = np.array(
            [np.NaN for _ in system.state_vars], dtype=np.float64
        )
        for k, v in x0.items():
            _, idx = str(k).split('_')
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

    def f(y):
        res = np.empty_like(X0)
        func(t0, X0, y, res)
        return res

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

    if system.ports:
        raise ModelException(
            "Cannot Simulate %s: unconnected ports %s",
            system, system.ports)

    if system.control_vars and not control_vars:
        raise ModelException("Control variable not specified")

    samples = int(1 / dt) + 1
    t = np.linspace(*timespan, samples)

    res, X =_bondgraph_to_residuals(system, control_vars)
    X0, DX0 = _fetch_ic(x0, dx0, system, res, t[0])

    solver_name = 'ida'
    dae_solver = dae(solver_name, res)
    sol = dae_solver.solve(t, X0, DX0)

    return t.reshape((samples, 1)), np.transpose(sol.values.y).T


def _to_function(string, X, DX, substitutions):
    f = sp.sympify(string).subs(substitutions)

    f_n = sp.lambdify((sp.S('t'), X, DX), f, "numpy")
    return f_n


def _bondgraph_to_residuals(model, control_vars=None):
    dX = sp.IndexedBase('dX')
    X = sp.IndexedBase('X')
    U = sp.IndexedBase('U')
    x_subs = []
    dx_subs = []
    u_subs = []
    u_func = []
    n = len(model.state_vars)
    m = 0

    for i, x in enumerate(model.state_vars):
        x_subs.append((x, X[i]))
        dx_subs.append((sp.S(f'dx_{i}'), dX[i]))

    if len(model.control_vars) > 0:
        u_func_dict = {}
        u_constants = {}
        if isinstance(control_vars, list):
            u_func_dict.update({
                i: f for i, f in enumerate(control_vars)}
            )
        elif isinstance(control_vars, dict):
            u_func_dict.update({
                int(v[2:]): f for v, f in control_vars.items()
            })
        elif len(model.control_vars) == 1:
            u_func_dict[0] = control_vars
        else:
            raise TypeError(f"Control argument {control_vars} is invalid")

        test_x = np.zeros(shape=(n,), dtype=np.float32)
        for idx, f in u_func_dict.items():
            try:
                if isinstance(f, str):
                    f = _to_function(f, X, dX, dx_subs + x_subs)
                    u_func_dict[idx] = f
                if isinstance(f, (float, int, sp.Number)):
                    u_constants[idx] = f
                if n == 1:
                    r = f(0, 0, 0)
                else:
                    r = f(0, test_x, test_x)
                assert isinstance(r, (float, int, sp.Number)), "Invalid output from control"
            except Exception as ex:
                message = f"Invalid control function for var: {idx}.\n " \
                           "Control functions should be of the form:\n" \
                           f"{idx} = f(t, x, dx/dt)"
                raise ModelException(message)

        for i, u in enumerate(model.control_vars):
            if i in u_constants:
                u_subs.append((u, u_constants[i]))
                continue
            u_subs.append((u, U[m]))
            try:
                u_func.append(u_func_dict[i])
            except KeyError:
                raise ModelException(f"Control variable {u} must be specified")
            m += 1
    rels = [r.subs(dx_subs).subs(x_subs).subs(u_subs)
            for r in model.constitutive_relations]

    if len(rels) != n:
        raise ModelException("Model simplification error: system is under-determined")

    Fsym = sp.symarray('F', shape=n)
    for i, r in enumerate(rels):
        Fsym[i] = r

    t = sp.S('t')
    _r = np.empty(shape=(n,), dtype=np.float64)
    if not u_func:
        F = sp.lambdify((t, X, dX), Fsym)

        def residual(_t,_x,_dx,_res):
            _r = F(_t,_x, _dx)
            for i in range(n):
                _res[i] = _r[i]
    else:
        _u = np.empty(shape=(m,), dtype=np.float64)
        Fsym_u = sp.lambdify((t, X, dX, U), Fsym)

        def residual(_t, _x, _dx, _res):
            _u = [u_f(_t, _x, _dx) for u_f in u_func]
            _r = Fsym_u(_t, _x, _dx, _u)
            for i in range(n):
                _res[i] = _r[i]
    return residual, X

