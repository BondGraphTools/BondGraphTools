
import functools
import sympy as sp

def bondgraph_to_equations(graph):

    coords, params, relations, labels = _build_coordinates(graph)

    x, dx, u = coords
    A, B, G, F = _matricize(coords, params, relations)

    # Adx + Bx = Gu + F(x)

    return None


def _build_coordinates(graph):
    # coordinates
    E = []
    F = []
    P = []
    Q = []
    U = []
    params = []
    ctrl_params = {}
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
            en = sp.symbols("e_{{{}}}".format(label))
            fn = sp.symbols("f_{{{}}}".format(label))
            pn = sp.symbols("p_{{{}}}".format(label))
            qn = sp.symbols("q_{{{}}}".format(label))
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

        if node.control_variables:
            for v in node.control_variables:
                un = sp.symbols("u_{{{}}}".format(str(len(ctrl_params))))
                ctrl_params[node_id, v] = un
                subs += [(v, un)]
                U.append(un)

        if node.local_params:
            for param, v in node.local_params.items():
                if isinstance(v, int) or isinstance(v, float):
                    subs.append((param, v))
                else:
                    p_an = sp.symbols("{}_{{{}}}".format(param, node.name))
                    params_labels[(node.name, param)] = p_an
                    params.append(p_an)
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
    params += list(globals_map.values())
    DX = F + E
    X = Q + P

    return (X, DX, U), params, relations, (coord_labels, params_labels, globals_map)


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

def _matricize(coords, params, relations):

    x, dx, u = coords

    A_dict = {}
    B_dict = {}
    J_dict = {}
    G_dict = {}
    F_dict = {}

    de_row = 0
    j_row = 0
    for rel in relations:

        rel_is_linear = True
        a_terms = {}
        b_terms = {}
        coeffs = extract_coefficients(rel, constants=params)
        for c in coeffs:
            
            if c in dx:
                a_terms[dx.index(c)] = coeffs[c]
            elif c in x:
                b_terms[x.index(c)] = coeffs[c]
            elif c.atoms() & set(u):
                rel_is_linear = False
                try:
                    G_dict[(de_row, 0)] -= c*coeffs[c]
                except KeyError:
                    G_dict[(de_row, 0)] = -c * coeffs[c]
            else:
                rel_is_linear = False
                try:
                    F_dict[(de_row, 0)] -= c*coeffs[c]
                except KeyError:
                    F_dict[(de_row, 0)] = -c * coeffs[c]

        if rel_is_linear and (a_terms and not b_terms):
            for col, v in a_terms.items():
                J_dict[(j_row, col)] = v
            j_row += 1
        elif rel_is_linear and (b_terms and not a_terms):
            for col, v in b_terms.items():
                J_dict[(j_row, col)] = v
            j_row += 1

        else:
            for col, coeff in a_terms.items():
                A_dict[(de_row,col)] = coeff

            for col, coeff in b_terms.items():
                B_dict[(de_row,col)] = coeff
            de_row += 1

    A = sp.SparseMatrix(de_row, len(dx), A_dict)
    B = sp.SparseMatrix(de_row, len(x), B_dict)
    J = sp.SparseMatrix(j_row, len(x), J_dict)
    G = sp.SparseMatrix(de_row, 1, G_dict)
    F = sp.SparseMatrix(de_row, 1, F_dict)

    return A, B, J, G, F

def _construct_edge_subs(graph, coord_labels, coords, keep_dirs=False):
    """
    Generates the substitution list for the bond graph.

    When the bond graph comes in, it is specified in terms of the power across
    each port of each node. This method uses the bond wiring to determine a
    substitution list to chance the co-ordinate system so that it is in terms
    of the power balance for each bond; halving the co-ordinate space.

    Args:
        graph: The bond graph
        coord_labels: dictionary of co-ordinate information
        coords: The nodal co-ordinate space.
        keep_dirs (False): Set to true if to preserve arrowhead direciton

    Returns:
        a list of substition rules, and a new set of co-ordinates

    """

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

            # ed is  e-discard, ek is e-keep
            if n1 == from_node or keep_dirs:
                _,  (ed, fd, pd, qd) = coord_labels[(n1, po1)]
                _,  (ek, fk, pk, qk) = coord_labels[(n2, po2)]
            else:
                _, (ed, fd, pd, qd) = coord_labels[(n2, po2)]
                _, (ek, fk, pk, qk) = coord_labels[(n1, po1)]


            subs += [
                (ed, ek), (fd, -fk), (pd, pk), (qd, -qk)
            ]
            E.remove(ed)
            F.remove(fd)
            P.remove(pd)
            Q.remove(qd)

    return subs, (E, F, P, Q)

def smith_normal_form(matrix):
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
    Mp = sp.MutableSparseMatrix(0, n, [])
    row = 0
    ins = 0

    while row < m:
        col = row
        while col + ins < n and M[row, col + ins] == 0:
            col += 1
        if col >= row + ins:
            Mp = Mp.col_join(sp.zeros(col - row, n))
            ins += col - row
        Mp = Mp.col_join(M.row(row))
        row += 1

    m, n = Mp.shape

    if m < n:
        Mp = Mp.col_join(sp.zeros(n-m, n))
    return Mp





def is_linear(equation, basis):
    """
    Check to see if an equation is linear.

    Args:
        equation: Sympy Equation
        basis: A vector of basis co-ordinates

    Returns:
        True if the equation is linear, false if not

    """
    if isinstance(equation, sp.Add):
        for arg in equation.args:
            if not is_linear(arg, basis):
                return False
    elif isinstance(equation, sp.Mul):
        if sum([a.atoms() & set(basis) != set() for a in equation.args]) > 1:
            return False
        for arg in equation.args:
            if not is_linear(arg, basis):
                return False
    elif isinstance(equation, sp.Number) or isinstance(equation, sp.Symbol):
        return True
    else:
        return False
    return True




def linearise(equation, X, X0=None):
    """

    Args:
        equation:
        X:
        X0:

    Returns:

    """
    if not X0:
        X0 = [0 for _ in X]

    grad = [sp.diff(equation, x) for x in X]
    return [functools.reduce(lambda eq, pair: eq.subs(*pair),
            zip(X, X0), partial) for partial in grad]


def extract_coefficients(equaiton, constants=None):

    coeff_dict = {}

    for term in equaiton.expand().args:

        products = flatten(term.as_coeff_mul())

        coeff = 1
        base = 1

        for factor in products:
            if factor.is_number or factor in constants:
                coeff = sp.Mul(coeff, factor)
            else:
                base = sp.Mul(base, factor)

        coeff_dict[base] = coeff
    return coeff_dict


def flatten(sequence):
    for item in sequence:
        if isinstance(item, (list, tuple)):
            for subitem in flatten(item):
                yield subitem
        else:
            yield item