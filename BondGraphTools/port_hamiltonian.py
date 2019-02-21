
import sympy as sp

from .atomic import Component


class PortHamiltonian(Component):
    """
    Port Hamiltonians are specified by a energy storage function $H(x)$ which
    generate the dynamics of the state variables $x$ and the associated
    Dirac structure $(e, f)$.

    Reserved variable names:
        $x, y, z$ or indexed notations thereof.

    Variables in the reserved list are taken to be stateful.
    $x$,$y$,$z$ are taken to be in some sense canonical co-ordinates
    of some differentiable manifold. Where $p$ is taken to be in the
    co-tangent space in the usual manner.


    Args:
        hamiltonian (str): The Hamiltonian storage function from which to
                           generate the component.


    See Also: Component
    """
    __vars = {"x", "y", "z"}

    def __init__(self, hamiltonian,
                 *args, **kwargs):

        cr, x, params, ports = self._generate_relations(hamiltonian)
        kwargs["constitutive_relations"] = cr
        kwargs["state_vars"] = x

        if "params" not in kwargs:
            kwargs["params"] = {}
        else:
            for p in params:
                if p not in kwargs["params"] and str(p) not in kwargs["params"]:
                    kwargs["params"][str(p)] = None
        kwargs["ports"] = ports
        # for each state variable, assign
        # e_i = D_{x_i} H(x)
        # f_i = dx_i

        super().__init__(*args, **kwargs)
        self.hamiltonian = hamiltonian

    @staticmethod
    def _generate_relations(hamiltonian):
        Hx = sp.sympify(hamiltonian)
        params = {}
        state_vars = {}
        relations = []
        subs = []
        variables = sorted([str(a) for a in Hx.atoms() if not a.is_number])

        for atom in variables:
            if atom[0] not in PortHamiltonian.__vars:
                params[atom] = None
            elif atom[0] in ("q", "p"):
                raise ValueError("Co-ordinates q and p are protected."
                                 "Use indexed x, y instead")
            else:
                # Need to make sure we don't mess with the co_ordinates
                i = len(state_vars)
                state_vars[f"q_{i}"] = str(atom)
                subs.append((atom, sp.S(f"q_{i}")))

        Hx = Hx.subs(subs)

        for i in range(len(state_vars)):
            q = sp.S(f"q_{i}")
            # todo: this is dirty, fix me
            relations.append(str(Hx.diff(q).simplify() - sp.S(f"e_{i}")))
            relations.append(f"d{q} - f_{i}")

        ports = {i: None for i in range(len(state_vars))}

        return relations, state_vars, params, ports



