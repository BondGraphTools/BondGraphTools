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


def reduce(dx, x, A, B, J, NL):
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

    Returns:

    """

    nonlinearity = NL
    P = _smith_normal_form(J) # projection onto the nullspace of J
    n, _ = P.shape
    R = (sp.eye(n) - P)
    system = A * R * dx + B * R * x

    Rx = R*x
    for i in range(x.rows):
        nonlinearity.subs(x[i], Rx[i])

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

    J, _ = J.rref()
    D, _ = D.rref()
    A = D[:, 0:int(N/2)]
    B = D[:, int(N/2):]
    return ydot[int(N/2):,:], ydot[0:int(N/2),:], A, B, J, NL


def _smith_normal_form(matrix):
    """
    Assume n >= m
    Args:
        matrix:

    Returns:
    n x n smith normal form of the matrix.
    Particularly for projection onto the nullspace of M and the orthogonal
    complement
    that is, for a matrix M,
    P = _smith_normal_form(M) is a projection operator onto the nullspace of M
    """
    M, _ = matrix.rref()
    m, n = M.shape
    Mp = sp.Matrix(0, n, [])
    row = 0
    ins = 0

    while row < m:
        col = row
        while col + ins < n and M[row , col+ ins] == 0:
            col += 1
        if col > row + ins:
            Mp = Mp.col_join(sp.zeros(col - row, n))
            ins += col - row
        Mp = Mp.col_join(M.row(row))
        row += 1

    m, n = Mp.shape

    if m<n:
        Mp = Mp.col_join(sp.zeros(n-m, n))
    return Mp


def _construct_coordinates(graph):
    # coordinates
    E = []
    F = []
    P = []
    Q = []
    A = []

    relations = []

    n = 0
    globals_map = {}
    an = 0
    bn = 0
    coord_labels = {}
    params_labels = {}
    globals_labels = {}

    for node_id, node in graph.nodes.items():

        subs = []

        one_port = len(node.ports) == 1

        for port_id in node.ports:
            if one_port:
                label = node.name
            else:
                label = "{};{}".format(node.name, port_id)
            en = dynamicsymbols("e_{{{}}}".format(label))
            fn = dynamicsymbols("f_{{{}}}".format(label))
            pn = dynamicsymbols("p_{{{}}}".format(label))
            qn = dynamicsymbols("q_{{{}}}".format(label))
            coord_labels[(node_id, port_id)] = (node.name, [en, fn, pn, qn])
            E.append(en)
            F.append(fn)
            P.append(pn)
            Q.append(qn)
            subs += [("e_{}".format(port_id), en),
                    ("f_{}".format(port_id), fn),
                    ("p_{}".format(port_id), pn),
                    ("q_{}".format(port_id), qn)]

            n += 1

        if node.local_params:
            for param, v in node.local_params.items():
                if isinstance(v, int) or isinstance(v, float):
                    subs.append((param, v))
                else:
                    p_an = sp.symbols("{}_{{{}}}".format(param, node.name))
                    params_labels[(node.name, param)] = p_an
                    A.append(p_an)
                    subs.append((param, p_an))
                    an += 1


        if node.global_params:
            for param in node.global_params:
                try:
                    p_bn = globals_map[param]
                except KeyError:
                    p_bn = sp.symbols("{}".format(param))
                    globals_map[param] = p_bn

                subs.append((param, p_bn))

        node_cr = _build_relations(node)
        relations += [r.subs(subs) for r in node_cr]

    subs, (E, F, P, Q) = _construct_edge_subs(graph, coord_labels, (E, F, P, Q))

    relations = [r.subs(subs).simplify() for r in relations]
    B = list(globals_map.values())
    DX = F + E
    X = Q + P

    return DX, X, A + B, relations, (coord_labels, params_labels, globals_map)


def _build_relations(node):
    rels = []
    for string in node.constitutive_relations:
        iloc = 0
        iloc = string.find("_i", iloc)

        if iloc < 0:
            # just a plain old string; sympy can take care of it
            rels.append(sp.sympify(string))
            continue

        sloc = string.rfind("sum(", 0, iloc)

        if sloc < 0:
            # we have a vector equation here.
            for port_id in node.ports:
                rels.append(
                    sp.sympify(
                        string.replace("_i", "_{}".format(port_id))))
            continue

        else:
            tiers = 0

            next_open = string.find("(", sloc + 4)
            eloc = string.find(")", sloc + 4)
            while next_open > 0 or tiers > 0:
                if next_open < eloc:
                    tiers += 1
                    next_open = string.find("(", next_open)
                else:
                    tiers -= 1
                    eloc = string.find(")", eloc)

            if eloc < 0:
                raise ValueError("Unbalanced brackets", string)

            substring = string[sloc+4: eloc]
            terms = [substring.replace("_i", "_{}".format(p))
                     for p in node.ports]
            symstr = string[0:sloc] + "(" + " + ".join(terms) + string[eloc:]
            rels.append(
                sp.sympify(symstr)
            )

    return [r for r in rels if r != 0]


def _construct_edge_subs(graph, coord_labels, coords):
    E, F, P, Q = coords
    edge_list = [((n1, p1), (n2, p2)) for n1, n2, p1, p2 in graph.bonds]

    node_list = sorted([(k,v) for k,v in graph.nodes.items()],
                       key=lambda x: len(x[1].ports))
    subs = []

    while node_list:
        from_node, node = node_list.pop()
        edges = [e for e in edge_list if e[0][0] == from_node or e[1][0] == from_node]

        while edges:
            edge = edges.pop()
            (n1, po1), (n2, po2) = edge
            edge_list.remove(edge)

            if n1 == from_node:
                _,  (e1, f1, p1, q1) = coord_labels[(n1, po1)]
                _,  (e2, f2, p2, q2) = coord_labels[(n2, po2)]
                dirn = -1  # power leaving
            else:
                _, (e1, f1, p1, q1) = coord_labels[(n2, po2)]
                _, (e2, f2, p2, q2) = coord_labels[(n1, po1)]
                dirn = 1  # power coming in

            subs += [
                (e1, e2), (f1, -f2), (p1, p2), (q1, -q2)
            ]
            E.remove(e1)
            F.remove(f1)
            P.remove(p1)
            Q.remove(q1)

    return subs, (E, F, P, Q)