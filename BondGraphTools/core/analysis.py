from functools import reduce

import sympy as sp
import numpy as np
from scipy.optimize import fsolve

from .algebra import bondgraph_to_equations, null_space_basis

DTYPE = np.float64

class Simulation(object):

    def __init__(self, bond_graph, **kwargs):

        system, coords, params, labels = bondgraph_to_equations(bond_graph)

        self._x_s, self._dx_s, u = coords
        self._t_s = sp.symbols("t")
        self._u_s = {us: None for us in u}
        self._params = {p: None for p in params}

        self._A, self._B, J, self._G, self._F = system

        n, m = J.shape

        self._A = self._A.col_join(J)
        self._B = self._B.col_join(sp.SparseMatrix.zeros(n, m))
        self._G = self._G.col_join(sp.SparseMatrix.zeros(n, 1))
        self._F = self._F.col_join(sp.SparseMatrix.zeros(n, 1))

        for k, v in kwargs.items():
            if k in ("parameters", "params"):
                self.parameters = v
            elif k in ("control_variables", "cv"):
                self.control_variables = v

    @property
    def parameters(self):
        return {k:v for k,v in self._params.items()}

    @parameters.setter
    def parameters(self, params):
        if isinstance(params, list):
            for i, k in enumerate(self._params.keys()):
                self._params[k] = params[i]
        elif isinstance(params, dict):
            self._params.update(params)

    @property
    def control_variables(self):
        return {k:v for k,v in self._u_s.items()}

    @control_variables.setter
    def control_variables(self, values):
        if isinstance(values, list):
            for i, k in enumerate(self._u_s):
                fn = values[i]
                if isinstance(fn, str):
                    fn = sp.sympify(fn)

                for atom in fn.atoms():
                    if not atom.is_number and atom != self._t_s:
                        raise Exception

                self._u_s[k] = fn

    def run(self, tspan, h):


        A = self._A.subs(self._params.items())
        B = self._B.subs(self._params.items())
        F = self._F.subs(self._params.items())
        G = self._G



        return

    #     self._P = null_space_basis(J)
    #
    #     self._A = A*self._P
    #     self._B = B*self._P
    #
    #     self.__ys = sp.Matrix(
    #         [sp.symbols("y_{}".format(i)) for i in range(self._A.cols)]
    #     )
    #     self.__xs = sp.Matrix(
    #         [sp.symbols("x_{}".format(i)) for i in range(self._A.cols)]
    #     )
    #
    #     Pl = self._P*self.__xs
    #     Pdl = self._P*self.__ys
    #
    #     for i in range(len(self._x_s)):
    #         self._F = self._F.subs(
    #             self._x_s[i], Pl[i]).subs(
    #             self._dx_s[i], Pdl[i]
    #         )
    #         self._G = self._G.subs(self._x_s[i], Pl[i]).subs(
    #             self._dx_s[i], Pdl[i])
    #
    #     self.labels = labels
    #
    #     self._x0 = None
    #     self._dx0 = None
    #     self._u = None
    #
    #     self._t = None
    #     self._yn = None
    #     self._xn = None
    #
    #
    #     print()
    #
    #
    # def run(self, tspan, output_vars=None, dt=0.0001):
    #
    #     """
    #     Args:
    #         tspan(tuple): Tuple containing start and end times.
    #         output_vars: List containing the co-ordinates to produce.
    #
    #     Returns:
    #         None if no co-ordinates are requested, otherwise the timeseries
    #         tuple containing (t, x[])
    #     """
    #     return None
    #
    #
    # def _onestep_BDF(self, t_span, dt=0.001):
    #
    #     t_start, t_stop = t_span
    #
    #     A = self._A.subs(self._params.items())
    #     B = self._B.subs(self._params.items())
    #     F = self._F.subs(self._params.items())
    #     G = self._G
    #     DFx, DFy = self._generate_partials(F)
    #
    #     T = np.linspace(0.0, t_stop, int(1/dt))
    #
    #     h = 0.5*(T[1] - T[0])
    #     assert h > 0
    #     m = len(T)
    #     _, n = A.shape
    #     u = np.zeros((len(self._u_s), m), dtype=DTYPE)
    #
    #     X, Y = self._initialise(m)
    #
    #     for i, k in enumerate(self._u_s):
    #         for j in range(m):
    #             u[i, j] = DTYPE(
    #                 self._u_s[k].subs(self._t_s, T[j])
    #             )
    #
    #     args = list(self.__xs) + list(self.__ys) +list(self._u_s) + [self._t_s]
    #
    #     func = sp.lambdify(
    #         args,((A + DFy + h * (B + DFx)).inv()) * (
    #                 (h * (B + DFx) - DFy)*self.__ys + B * self.__xs + self._G)
    #     )
    #     print(func(*args))
    #     for i, t in enumerate(T[1:]):
    #         arg = X[:, i-1].tolist() + Y[:, i-1].tolist() + u[:, i].tolist() + [t]
    #         v = func(*arg)
    #         Y[:, i] = v.T
    #         X[:, i] = X[:, i-1] + h*(Y[:,i] + Y[:,i-1])
    #
    #     return X, Y, T
    #
    #
    #     # algo:
    #     # For a given x_n,t
    #     # solve for y_{n+1}; (A + dt/2 B) y +Bx Gu + F(dty/2 ) = 0
    #     # update y_{n+1} = 2*dt*x_n + y_n
    #
    # def _initialise(self, m):
    #     n = self._A.cols
    #     return np.zeros((n,m), dtype=DTYPE), np.zeros((n,m), dtype=DTYPE)
    #
    #
    #
    # def _generate_partials(self, F=None):
    #
    #     rows = self._A.rows
    #     cols = self._A.cols
    #     DFy = sp.MutableSparseMatrix.zeros(rows, cols)
    #     DFx = sp.MutableSparseMatrix.zeros(rows, cols)
    #     if not F:
    #         F = self._F
    #     for row in range(rows):
    #         for col in range(cols):
    #             DFy[row, col] = F[row].diff(self.__ys[col])
    #             DFx[row, col] = F[row].diff(self.__xs[col])
    #
    #     return DFx, DFy
    #

