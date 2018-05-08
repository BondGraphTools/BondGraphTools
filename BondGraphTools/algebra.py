import sympy as sp


def extract_coefficients(equation, local_map, global_coords):

    coeff_dict = {}
    nonlinear_terms = sp.S(0)
    subs = [(k, global_coords[v]) for k, v in local_map.items()]

    subs.sort(key=lambda x: str(x[1])[-1], reverse=True)

    for term in equation.expand().args:
        prod_iter = flatten(term.as_coeff_mul())
        coeff = next(prod_iter)
        base = list(prod_iter)
        if not base:
            coeff_dict[len(global_coords)-1] = coeff
        elif len(base) == 1 and base[0] in local_map:
            coeff_dict[local_map[base[0]]] = coeff
        else:
            new_term = term
            new_term = new_term.subs(subs)
            nonlinear_terms = sp.Add(new_term, nonlinear_terms)

    return coeff_dict, nonlinear_terms


def flatten(sequence):
    for item in sequence:
        if isinstance(item, (list, tuple)):
            for subitem in flatten(item):
                yield subitem
        else:
            yield item


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
    # M, _ = matrix.rref()
    # m, n = M.shape
    # M = sp.SparseMatrix(m, n, M)
    # m_dict = {}
    # current_row = 0
    #
    # row_map = {}
    #
    # current_row = 0
    #
    # for row, c_idx, entry in M.RL:
    #     if row not in row_map:
    #         row_map[row] = c_idx
    #         r_idx = c_idx
    #
    #     else:
    #         r_idx = row_map[row]
    #
    #     m_dict[(r_idx, c_idx)] = entry
    #
    # return sp.SparseMatrix(n, n, m_dict)

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
        Mp = Mp.col_join(sp.zeros(n - m, n))
    elif m > n and Mp[:, n:m].is_zero:
        Mp = Mp[:n, :n]

    return Mp


def adjacency_to_dict(nodes, edges, offset=0):
    """
    matrix has 2*#bonds rows
    and 2*#ports columes
    so that MX = 0 and X^T = (e_1,f_1,e_2,f_2)

    Args:
        index_map: the mapping between (component, port) pair and index

    Returns: Matrix M

    """
    M = dict()

    for i, (node_1, node_2) in enumerate(edges):
        j_1 = offset + 2 * nodes[node_1]
        j_2 = offset + 2 * nodes[node_2]
        # effort variables
        M[(2 * i, j_1)] = - 1
        M[(2 * i, j_2)] = 1
        # flow variables
        M[(2 * i + 1, j_1 + 1)] = 1
        M[(2 * i + 1, j_2 + 1)] = 1

    return M


def inverse_coord_maps(tangent_space, port_space, control_space):
    inverse_tm = {
        coord_id: index for index, coord_id
        in enumerate(tangent_space.values())
    }
    inverse_js = {
        coord_id: index for index, coord_id
        in enumerate(port_space.values())
    }
    inverse_cm = {
        coord_id: index for index, coord_id
        in enumerate(control_space.values())
    }

    coordinates = [dx for _, dx in tangent_space]

    for e, f in port_space:
        coordinates += [e, f]
    for x, _ in tangent_space:
        coordinates.append(x)
    for u in control_space:
        coordinates.append(u)
    coordinates.append(sp.S("c"))

    return (inverse_tm, inverse_js, inverse_cm), coordinates


def generate_relations(coords, matrix, nonlinear_op):

    """
    System in the form

    Ax = 0
    F(x) = 0

    Args:
        coords:
        matrix:
        nonlinear_op:

    Returns:

    """

