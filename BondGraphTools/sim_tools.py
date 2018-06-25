import numpy as np
import sympy as sp

from .exceptions import ModelException
from .algebra import inverse_coord_maps, create_ds
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
             timespan,
             x0,
             dx0=None,
             control_vars=None):
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

    if system.control_vars and not control_vars:
        raise ModelException("Control variable not specified")

    if not de:
        start_julia()


    tspan = tuple(float(t) for t in timespan)
    X0 = np.array(x0, dtype=np.float64)
    assert len(X0) == len(system.state_vars)
    try:
        func = _build_ode(system, control_vars)
        problem = de.ODEProblem(func, X0, tspan)
    except NotImplementedError:
        func, diffs = _build_dae(system, control_vars)
        if dx0:
            DX0 = np.array(dx0, dtype=np.float64)
        else:
            DX0 = np.zeros(X0.shape, dtype=np.float64)
        problem = de.DAEProblem(func, DX0, X0, tspan, differential_vars=diffs)

    sol = de.solve(problem)
    t = np.transpose(sol.t)

    return np.resize(t, (len(t), 1)), np.transpose(sol.u).T


def _build_ode(system, control_vars=None):
    coords, mappings, linear, nonlinear, constraints = system.system_model()

    ss_map, js_map, cv_map = mappings
    m = len(ss_map)
    offset = m + 2*len(js_map)

    A = linear[0:m, 0:m]
    B = linear[0:m, m:offset]

    if not B.is_zero or not (A - sp.eye(m)).is_zero:
        raise NotImplementedError("DAE's not yet implemented")
    x, subs, string_subs = _generate_cv_subs(mappings, control_vars)

    L = -linear[0:m, offset:offset+m]

    Lu = linear[0:m, offset+m:].dot(coords[offset + m:])

    if isinstance(Lu, sp.Symbol):
        Lu = [Lu]

    Nu = nonlinear[0:m, :]
    N = [-sp.Add(left, right).subs(subs) for left, right in zip(Lu, Nu)]

    # DX = LX + N(X, t)
    julia_string = """function dxdt(dX, X, p, t)\n"""

    for var, var_string in string_subs.items():
        julia_string += f"    {var} = {var_string}\n"

    for i in range(m):
        julia_string += f"    dX[{i+1}] ="
        lx = sp.simplify(L[i,:].dot(x))
        nl = sp.sympify(N[i])

        if lx:
            julia_string += f"{repr(lx)}"

        julia_string += repr(nl) + "\n"

    julia_string += "end"
    julia_string = julia_string.replace("**","^")
    func = j.eval(julia_string)

    return func


def _generate_cv_subs(mappings, control_vars=None):

    ss_map, js_map, cv_map = mappings
    m = len(ss_map)
    k = len(cv_map)

    x = [sp.symbols(f"x_{i}") for i in range(m)]
    X = [sp.symbols(f"X[{i+1}]") for i in range(m)]
    subs = list(zip([sp.symbols(f"dx_{i}") for i in range(m)],
                    [sp.symbols(f"dX[{i+1}]") for i in range(m)]))

    subs += list(zip(x, X))
    t = sp.S('t')
    string_subs = {}

    if isinstance(control_vars, (float, int, complex)) and\
            len(k) == 1:
        subs += [
            (sp.Symbol('u_0'), control_vars),
            (sp.Symbol('du_0'), 0)
        ]

    elif isinstance(control_vars, list) and len(control_vars) == k:
        for i, cv_string in enumerate(control_vars):
            u = sp.Symbol(f'u_{i}')
            du = sp.Symbol(f'du_{i}')
            try:
                fx = sp.sympify(cv_string)
                dfx = sum(fx.diff(x_i) for x_i in x) + fx.diff(t)
                subs.append((u, fx))
                subs.append((du, dfx))
            except sp.SympifyError:
                u_str = f"u{i+1}"
                for i in reversed(range(m)):
                    cv_string = cv_string.replace(
                        f"x_{i}", f"X[{i+1}]"
                    )
                string_subs[u_str] = cv_string
                subs.append(sp.symbols(f"u_{i}, {u_str}"))

    elif isinstance(control_vars, dict):
        for key, cv_string in control_vars:
            u = sp.Symbol(f'{key}')
            du = sp.Symbol(f'd{key}')
            try:
                fx = sp.sympify(cv_string)
                dfx = sum(fx.diff(x_i) for x_i in x) + fx.diff(t)
                pair = [(du, dfx), (u, fx)]
                for pp in pair:
                    subs = [
                        s.subs(pp) for s in subs
                    ]
                    subs.append(pp)

            except sp.SympifyError:
                u_str = f"u{len(string_subs)}"
                for i in reversed(range(m)):
                    cv_string = cv_string.replace(
                        f"x_{i}", f"X[{i+1}]"
                    )

                string_subs[u_str] = cv_string
                subs.append(sp.symbols(f"{key}, {u_str}"))
    else:
        raise ValueError("Invalid control variables: %s", repr(control_vars))

    return X, subs, string_subs


def _build_dae(system, control_vars=None):

    mappings, coords = inverse_coord_maps(*system.basis_vectors)
    ss_map, js_map, cv_map = mappings

    m = len(ss_map)

    if len(js_map) > 0:
        raise NotImplementedError("Bond Graph has unconnected Ports")

    derivatives = set(coords[0:m])
    differential_vars = []

    # construct julia coords

    # x = [sp.symbols(f"x_{i}") for i in range(m)]
    # subs = list(zip([sp.symbols(f"dx_{i}") for i in range(m)],
    #             [sp.symbols(f"dX[{i+1}]") for i in range(m)]))
    #
    # subs += list(zip(x, [sp.symbols(f"X[{i+1}]") for i in range(m)]))

    # subs, cv_text = _generate_cv_subs(control_vars, subs)
    x, subs, string_subs = _generate_cv_subs(mappings, control_vars)

    julia_string = "function f(dX, X, p, t)\n"
    end_string = "    return ["
    i = 0

    for var, var_string in string_subs.items():
        julia_string += f"    {var} = {var_string}\n"

    for relation in system.constitutive_relations:
        r = relation.subs(subs)
        if not r:
            continue

        differential_vars.append(
            derivatives & relation.atoms() != set()
        )

        temp_string = str(r)

        julia_string += f"    res{i+1} = {temp_string}\n"

        if i > 0:
            end_string += ', '
        end_string += f"res{i+1}"
        i += 1
    assert len(differential_vars) == i
    end_string += "]\nend"
    julia_string += end_string
    logger.warning("Julia Function")
    logger.warning(julia_string)
    if not j: start_julia()

    func = j.eval(julia_string)

    return func, differential_vars



def julia():
    global j


class Simulation(object):
    def __init__(self, model,
                 timespan=None,
                 x0=None,
                 dx_0=None,
                 control_vars=None):

        coords, mapping, linear, nonlinear, constraints = model.system_model(
            control_vars=control_vars
        )

        self._state_map, self._port_map, self._cv_map = mapping

        d_system, port_func, constraints = create_ds(
            coords, mapping, linear, nonlinear, constraints
        )

        self._solver = None
        self._julia_func = None

    def run(self, x0, timespan):
        pass




