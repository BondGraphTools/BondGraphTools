
import logging
import sympy as sp

from .exceptions import SymbolicException

logger = logging.getLogger(__name__)

class NotInvertible(Exception):
    pass


def extract_coefficients(equation, local_map, global_coords):

    coeff_dict = {}
    nonlinear_terms = sp.S(0)
    subs = [(k, global_coords[v]) for k, v in local_map.items()]

    subs.sort(key=lambda x: str(x[1])[-1], reverse=True)
    logger.info("Extracting coefficients from %s", repr(equation))
    logger.info("Using local-to-global substitutions %s", repr(subs))

    terms = equation.expand().args
    if not terms:
        if equation in local_map:
            coeff_dict[local_map[equation]] = sp.S(1)
        elif equation.is_number:
            coeff_dict[len(global_coords)-1] = equation
        else:
            nonlinear_terms = equation
    else:
        for term in terms:
            factors = list(flatten(term.as_coeff_mul()))
            logger.info("Factors: %s", repr(factors))
            coeff = sp.S(1)
            base = []
            while factors:
                factor = factors.pop()
                if factor.is_number:
                    coeff *= factor
                else:
                    base.append(factor)
            logger.info("Base: %s", repr(base))
            if not base:
                coeff_dict[len(global_coords)-1] = coeff
            elif len(base) == 1 and base[0] in local_map:
                coeff_dict[local_map[base[0]]] = coeff
            else:
                new_term = term
                new_term = new_term.subs(subs)
                nonlinear_terms = sp.Add(new_term, nonlinear_terms)

    logger.info("Linear terms: %s", repr(coeff_dict))
    logger.info("Nonlinear terms: %s", repr(nonlinear_terms))

    return coeff_dict, nonlinear_terms


def _handle_constraints(linear_op, nonlinear_op, coordinates,
                        size_tuple):
    """
    Args:
        linear_op: Linear part of the constitutive relations. Assumed to be
        a matrix in smith normal form.
        nonlinear_op: The corresponding nonlinear part; a symbolic vector with
        the same number of rows.
        coordinates:
        size_tuple:

    Returns:

    """

    rows_added = 0
    added_cvs = []

    logger.info("Handling algebraic constraints")

    ss_size, js_size, cv_size, n = size_tuple
    offset = 2*js_size + ss_size
    for row in range(offset, linear_op.rows):
        logger.info("Testing row %s: %s + %s", repr(row),
                    repr(linear_op[row, :].dot(coordinates)),
                    repr(nonlinear_op[row]) if nonlinear_op else '')

        if nonlinear_op:
            nonlinear_constraint = nonlinear_op[row]
        else:
             nonlinear_constraint = None
        if linear_op[row, offset:-1].is_zero and not nonlinear_constraint:
            continue

        state_constraint = linear_op[row, offset: offset + ss_size]
        control_constraint = linear_op[row, offset + ss_size:-1]

        if nonlinear_constraint:
            logger.warning("Nonlinear constraints not yet implemented")

            iX = set(coordinates[0:offset+ss_size]) & nonlinear_constraint.atoms()

            # case 1: constraint is f(x,u, c)= 0 for ONE x in iX
            #         and f(x) is invertible;

            # case 2: function




        row = state_constraint.row_join(sp.SparseMatrix(1, offset + cv_size, {}))

        cv_dict = {}
        if not control_constraint.is_zero:
            logging.info("Found higher order control constraint")
            for cv_col in range(control_constraint.cols):
                const = control_constraint[cv_col]
                if not const:
                    continue

                try:
                    idx = added_cvs.index(cv_col)
                except ValueError:
                    idx = len(added_cvs)
                    added_cvs.append(cv_col)
                    linear_op= linear_op.col_insert(
                        -1, sp.SparseMatrix(linear_op.rows, 1, {}))
                    coord = coordinates[offset + ss_size + cv_col]
                    d_coord = sp.Symbol(f"d{str(coord)}")
                    coordinates.insert(-1, d_coord)
                    cv_size += 1
                    n+=1
                cv_dict[(0,idx)] = const

        row = row.row_join(sp.SparseMatrix(1, len(added_cvs) + 1, cv_dict))

        linear_op = linear_op.col_join(row)
        rows_added += 1

    if rows_added:
        if nonlinear_op:
            nonlinear_op = nonlinear_op.col_join(
                sp.SparseMatrix(rows_added, 1, {})
        )
            linear_op, nonlinear_op = smith_normal_form(linear_op, nonlinear_op)
        else:
            linear_op = smith_normal_form(linear_op)

    return coordinates, linear_op, nonlinear_op


def flatten(sequence):
    for item in sequence:
        if isinstance(item, (list, tuple)):
            for subitem in flatten(item):
                yield subitem
        else:
            yield item


def augmented_rref(matrix, augment=None):

    if augment:
        assert matrix.rows == augment.rows

    pivot = 0

    for col in range(matrix.cols):

        if matrix[pivot, col] == 0:
            j = None
            v_max = 0
            for row in range(pivot, matrix.rows):
                val = matrix[row, col]
                v = abs(val)
                if v > v_max:
                    j = row
                    v_max = v
            if not j:
                continue  # all zeros below, skip on to next column
            else:
                matrix.row_swap(pivot, j)
                if augment:
                    temp = augment[pivot, :]
                    augment[pivot, :] = augment[j, :]
                    augment[j, :] = temp

        a = matrix[pivot, col]
        for i in range(matrix.rows):
            if i != pivot and matrix[i, col] != 0:
                b = matrix[i, col]/a
                matrix[i, :] = matrix[i, :] - b * matrix[pivot, :]
                if augment:
                    augment[i, :] = augment[i, :] - b * augment[pivot, :]

        matrix[pivot, :] = matrix[pivot, :] / a
        if augment:
            augment[pivot, :] = augment[pivot, :] / a

        pivot += 1
        if pivot >= matrix.rows:
            break

    return matrix, augment


def smith_normal_form(matrix, augment=None):
    """
    Assume n >= m
    Args:
        matrix:
        augment:

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

    if augment:
        M, augment = augmented_rref(matrix.copy(), augment.copy())
    else:
        M, _ = augmented_rref(matrix.copy())
    # else:
    #     M, _ = matrix.rref()

    m, n = M.shape

    Mp = sp.MutableSparseMatrix(n, n, {})

    if augment:
        k = augment.cols
        Ap = sp.MutableSparseMatrix(n, k, {})

    for row in range(m):
        for col in range(row, n):
            if M[row, col] != 0:
                Mp[col, :] = M[row, :]
                if augment:
                    Ap[col, :] = augment[row, :]
                break
    if augment:
        return Mp, Ap
    else:
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

    coordinates.append(sp.S(1))

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

