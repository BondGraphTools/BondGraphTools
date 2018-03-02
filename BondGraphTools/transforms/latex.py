import functools
import sympy as sp
from sympy.physics.mechanics import dynamicsymbols


def bondgraph_to_sympy(graph):
    """

    Args:
        graph:

    Returns:
        (dx, x, J, A, B , NL)

        So that the bond graph's motion is given by
        A*dx + B*x = 0
        Jx = 0
        NL(x,dx) = 0
    """
    xdot, x, p, relations, labels = _construct_coordinates(graph)

    #dy, y, J, D,  NL = _matricize(xdot, x, p)
    return _matricize(xdot, x, relations)


def reduce(dx, x, A, B, J, NL, targets=None):
    """
    Projects the system onto the dynamics of the linear subspace spanned by
    the solutions to Jx = 0

    Assumptions:
        A, B are in rref
    Args:
        dx:
        x:
        A:
        B:
        J:
        NL:
        targets:

    Returns:

    """
    n = J.cols
    permute = sp.eye(n)

    if not targets:
        targets = []

    # Permute the co-ordinates so that the matrix perm mapping
    # y = perm*x
    # now has the target co-ordinates in the last positions
    # and derivatives at the top.

    top_swaps = 0
    bot_swaps = n

    for i, (xi, dxi) in enumerate(zip(x, dx)):
        if not A[:, i].is_zero:
            permute.row_swap(top_swaps, i)
            top_swaps += 1
        # we then find the projection onto the nullspace in these new
        # co-ordinates and solve for the orthognal projector.

    P = _smith_normal_form(J * permute.T)
    R = (sp.eye(n) - P)

    nonlinearity = NL
    dy = permute*dx
    y = permute*x
    system = A * permute * R * dy + B * permute * R * y

    Ry = R*y
    for i in range(y.rows):
        nonlinearity.subs(y[i], Ry[i])

    system.simplify()
    nonlinearity.simplify()
    return system, nonlinearity


def _matricize(xdot, x, relations):
    t = sp.symbols("t")

    Y = sp.Matrix(len(x) + len(xdot), 1, x + xdot)
    ydot = sp.Matrix(len(x) + len(xdot), 1, xdot + [sp.diff(dx, t) for dx in xdot])
    N = len(Y)
    NL = sp.Matrix(0, 1, [])
    D = sp.Matrix(0, N, [])
    J = sp.Matrix(0, N/2, [])
    for rel in relations:
        H = sp.hessian(rel, Y)
        if H.is_zero:
            #linear
            row_x = sp.Matrix(1, N/2, [sp.diff(rel, x1) for x1 in x])
            row_dotx = sp.Matrix(1, N/2, [sp.diff(rel, x1) for x1 in xdot])
            if row_x.is_zero:
                J = J.col_join(row_dotx)
            elif row_dotx.is_zero:
                J = J.col_join(row_x)
            else:
                D = D.col_join(
                    row_dotx.row_join(row_x)
                )
        else:
            # if H[0:int(N/2), 0:int(N/2)].is_zero:
            #     row = sp.Matrix(1, 1, [rel])
            # else:
            #     df = sp.Matrix([sp.diff(rel, y) for y in Y])
            #     row = sp.Matrix(1, 1, [df.dot(ydot)])
            # print(row)
            # print(NL)
            NL = NL.col_join(sp.Matrix(1,1,[rel]))

    # J, _ = J.rref()
    # D, _ = D.rref()
    A = D[:, 0:int(N/2)]
    B = D[:, int(N/2):]
    return ydot[int(N/2):, :], ydot[0:int(N/2), :], A, B, J, NL


