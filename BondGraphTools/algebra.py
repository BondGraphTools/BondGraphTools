import sympy as sp


def extract_coefficients(equation, local_map, global_coords):

    coeff_dict = {}
    nonlinear_terms = sp.S(0)
    subs = [(k, global_coords[v]) for k,v in local_map.items()]
    print(local_map)
    for term in equation.expand().args:
        prod_iter = flatten(term.as_coeff_mul())
        coeff = next(prod_iter)
        base = list(prod_iter)
        if not base:
            coeff_dict[len(global_coords)-1] = coeff
        elif len(base) == 1 and base[0] in local_map:
            coeff_dict[local_map[base[0]]] = coeff
        else:
            nonlinear_terms += term.subs(subs)

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
    M, _ = matrix.rref()
    n = max(M.shape)

    m_dict = {}
    current_row = 0

    for row, col, entry in M.RL:
        if row > current_row:
            current_row = col

        m_dict[(current_row, col)] = entry

    return sp.SparseMatrix(n, n, m_dict)
