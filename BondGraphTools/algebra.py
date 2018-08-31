
import logging
import sympy as sp
import numpy as np
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
            if len(base) == 1 and base[0] in local_map:
                coeff_dict[local_map[base[0]]] = coeff
            else:
                new_term = term
                new_term = new_term.subs(subs)
                nonlinear_terms = sp.Add(new_term, nonlinear_terms)

    logger.info("Linear terms: %s", repr(coeff_dict))
    logger.info("Nonlinear terms: %s", repr(nonlinear_terms))

    return coeff_dict, nonlinear_terms


def _generate_substitutions(linear_op, nonlinear_op, constraints, coords, size_tup):

    # Lx + F(x) = 0 =>  Ix = (I - L)x - F(x) = Rx - F(x)
    # Since L is in smith normal form (rref and square)
    # If (Rx)_{ii} = 0, and F(x)_i doesn't depend upon x_i
    # then we have x_i = (Rx)_i - F_i(x)
    c_atoms = set(coords)
    atoms = nonlinear_op.atoms() & c_atoms
    for constraint in constraints:
        atoms |= (constraint.atoms() & c_atoms)

    if not atoms:
        logger.info("No substitutions required")
        return []

    Rx = (sp.eye(linear_op.rows) - linear_op)

    ss_size, js_size, cv_size, n = size_tup
    substitutions = []
    coords_vect = sp.Matrix(coords)
    for i in reversed(range(2*(ss_size  + js_size))):
        co = coords[i]
        if Rx[i,i] == 0 and co in atoms and not co in nonlinear_op[i].atoms():

            eqn = (Rx[i,:]*coords_vect)[0] - nonlinear_op[i]
            pair = (coords[i], eqn)
            logger.info("Generating substition %s = %s",
                        repr(coords[i]), repr(eqn))
            substitutions = [
                (c, s.subs(*pair)) for c, s in substitutions
            ]
            substitutions.append(pair)

    return substitutions


def _process_constraints(linear_op,
                         nonlinear_op,
                         constraints,
                         coordinates,
                         size_tup):

    initial_constraints = []
    ss_size, js_size, cv_size, n = size_tup
    offset = 2 * js_size + ss_size

    coord_atoms = set(coordinates[0:offset+ss_size])

    coords_vect = sp.Matrix(coordinates)
    cv_constraints = list(
        linear_op[offset+ss_size:n,:]*coords_vect +
        nonlinear_op[offset+ss_size:n,0]
    )
    constraints += [cons for cons in cv_constraints if cons]
    linear_op = linear_op[:offset+ss_size, :]
    nonlinear_op = nonlinear_op[:offset+ss_size, :]

    while constraints:
        constraint, _ = sp.fraction(constraints.pop())
        logger.info("Processing constraint: %s",repr(constraint))
        atoms = constraint.atoms() & set(coord_atoms)

        # todo: check to see if we can solve f(x) = u => g(u) = x
        if len(atoms) == 1:
            c = atoms.pop()
            logger.info("Attempting to find inverse")
            solns = list(sp.solveset(constraint, c))
            if len(solns) == 1:
                idx = coordinates.index(c)
                sol = solns.pop()

                linear_op = linear_op.col_join(
                    sp.SparseMatrix(1, linear_op.cols, {(0, idx): 1})
                )
                nonlinear_op = nonlinear_op.col_join(
                    sp.SparseMatrix(1, 1, {(0,0): -sol})
                )
                constraint = c - sol
        else:
            logger.warning("..skipping %s", repr(constraint))
            initial_constraints.append(constraint)
        try:
            partials = [constraint.diff(c) for c in coordinates]
        except Exception as ex:
            logger.exception("Could not differentiate %s with respect to %s",
                         repr(constraint),repr(coordinates)
             )
            raise ex

        if any(p != 0 for p in partials[0:offset]):
            logger.warning("Cannot yet reduce order of %s", repr(constraint))
            initial_constraints.append(constraint)
        else:
            ss_derivs = partials[offset: offset + ss_size]
            cv_derivs = partials[offset + ss_size:]
            factor = 0
            lin_dict = {}
            nlin = 0
            for idx, coeff in enumerate(ss_derivs):
                if factor == 0 and coeff != 0:
                    factor = 1 / coeff
                    lin_dict.update({(0, idx): 1})
                elif factor != 0 and coeff != 0:
                    new_coeff = sp.simplify(coeff / factor)
                    if new_coeff.is_number:
                        lin_dict.update({(0, idx): new_coeff})
                    else:
                        nlin += new_coeff * coordinates[idx]
            if factor == 0:
                raise SymbolicException("Invalid Constraint %s",
                                        repr(constraint))
            for idx, coeff in enumerate(cv_derivs):
                if coeff != 0:
                    cv = coordinates[offset+ss_size+idx]
                    dvc = sp.Symbol(f"d{str(cv)}")
                    try:
                        dc_idx = coordinates.index(dvc)
                    except ValueError:
                        dc_idx = len(coordinates)
                        coordinates.append(dvc)
                        cv_size += 1
                        n += 1
                        linear_op = linear_op.row_join(
                            sp.SparseMatrix(linear_op.rows, 1, {})
                        )
                    eqn = coeff/factor
                    if eqn.is_number:
                        lin_dict.update({(0, dc_idx): eqn})
                    else:
                        nlin += eqn*dvc
            linear_op = linear_op.col_join(
                sp.SparseMatrix(1,linear_op.cols, lin_dict)
            )
            nonlinear_op = nonlinear_op.col_join(
                    sp.SparseMatrix(1,1,{(0,0):nlin})
            )

    linear_op, nonlinear_op, new_constraints = smith_normal_form(
        matrix=linear_op,
        augment=nonlinear_op)

    return linear_op, nonlinear_op, new_constraints + initial_constraints, \
           coordinates, (ss_size, js_size, cv_size, n)


def create_ds(coords, mapping, linear, nonlinear, constraints):

    ss_size = len(mapping[0])
    js_size = len(mapping[1])

    ##
    # We'd hope that the Linear Operator is in block form
    #
    # L =  [[A_1 B_1 C_1],
    #       [0   B_2 C_2],
    #
    # X = [[dx, l, x]]
    #
    # So that L.dot(X) + F(X,t) = 0
    #

    A_1 = linear[0:ss_size, 0:ss_size]
    B_1 = linear[0:ss_size, ss_size: 2*js_size + ss_size]
    C_1 = linear[0:ss_size, 2*js_size + ss_size:2*(js_size + ss_size)]
    D_1 = linear[0:ss_size, 2*(js_size + ss_size):]
    F_1 = nonlinear[0:ss_size,:]

    B_2 = linear[ss_size: 2*js_size + ss_size, ss_size: 2*js_size + ss_size]
    C_2 = linear[ss_size: 2 * js_size + ss_size,
          2 * js_size + ss_size:2 * (js_size + ss_size)]
    D_2 = linear[ss_size: 2 * js_size + ss_size, 2 * (js_size + ss_size):]
    F_2 = nonlinear[ss_size: 2 * js_size + ss_size, :]

    assert (B_2 - sp.eye(B_2.rows)).is_zero
    assert (A_1 - sp.eye(A_1.rows)).is_zero

    ds = (
        np.array(A_1).astype(np.float64),
        np.array(C_1).astype(np.float64),
        np.array(D_1).astype(np.float64)
    )

    def port_func(x,t):

        pass


    # dX = -C_1*X + -D_1*U - F_1(X)
    # J = -C_2*X - -D_2*U - F_2(X)
    #
    # x=(y,z) s.t. [dy,dz] = [0, f(y,z,u)]
    # set x = [[ M_1^T ]  (y
    #          [ M_2^T ]]  z)
    #


def _generate_cv_substitutions(subs_pairs, mappins, coords):
    state_map, port_map, control_map = mappins
    ss_size = len(state_map)

    cv_offset = 2*(ss_size + len(port_map))

    control_vars = {str(c) for c in coords[cv_offset:]}
    print(control_vars)
    subs = []
    for var, fx_str in subs_pairs.items():

        if var in control_vars:
            u = sp.S(var)
        elif var in control_map:
            u = sp.S(f"u_{control_map[var]}")
        else:
            raise SymbolicException("Could not substitute control variable %s",
                                    str(var))
        fx = sp.sympify(fx_str)

        subs.append((u, fx))

    return subs


def reduce_model(linear_op, nonlinear_op, coordinates, size_tuple,
                 control_vars=None):
    """
    Args:
        linear_op: Linear part of the constitutive relations.
        nonlinear_op: The corresponding nonlinear part; a symbolic vector with
        the same number of rows.
        coordinates: a list of all the relevant co-ordinates
        size_tuple:

    Returns:


    Todo:
        refactor so as to remove size_tuple;
        the co-ordinates should arrive partitioned!
    """
    #
    linear_op, nonlinear_op, constraints = smith_normal_form(
        matrix=linear_op,
        augment=nonlinear_op)

    rows_added = 0
    added_cvs = []
    cv_diff_dict = {}
    lin_dict = {}
    nlin_dict = {}

    logger.info("Handling algebraic constraints")

    ###
    # First; take care of control variables
    #

    #
    # Then substitute as much of the junction space as possible.
    #

    subs_list = _generate_substitutions(
        linear_op, nonlinear_op, constraints, coordinates, size_tuple
    )
    logger.info("Applying substitutions")

    nonlinear_op = nonlinear_op.subs(subs_list)
    constraints = [c.subs(subs_list) for c in constraints]

    logger.info("Reducing purely algebraic constraints")
    # second, reduce the order of all nonlinear constraints
    linear_op, nonlinear_op, constraints, coordinates, size_tuple =\
        _process_constraints(linear_op, nonlinear_op,
                             constraints, coordinates, size_tuple)
    logger.info("Applying substitutions, round 2")
    subs_list = _generate_substitutions(
        linear_op, nonlinear_op, constraints, coordinates, size_tuple
    )
    nonlinear_op = nonlinear_op.subs(subs_list)
    constraints = [c.subs(subs_list) for c in constraints]
    ##
    # Split the constraints into:
    # - Linear constraints; ie Lx = 0
    # - Nonlinear Constraints Lx + F(x) = 0
    #
    # Linear constraints are rows with more than 1 non-zero
    # that are not in the derivative subspace, and have a zero nonlinear part
    #
    # ## New Code
    ss_size, js_size, cv_size, n = size_tuple
    offset = 2 * js_size + ss_size
    for row in reversed(range(linear_op.rows, offset)):
        atoms = nonlinear_op[row].atoms()
        if not atoms & set(coordinates) and linear_op[row].nnz() > 1:
            logger.info("Linear constraint in row %s", repr(row))
            for idx in range(ss_size):
                v = linear_op[row, idx + offset]
                if v:
                    lin_dict.update({(rows_added,idx): v})
            for idx in range(cv_size):
                v = linear_op[row, idx + offset+ss_size]
                if v:
                    cv_diff_dict.update({(rows_added, idx): v})

    for row in range(offset, linear_op.rows):
        logger.info("Testing row %s: %s + %s", repr(row),
                    repr(linear_op[row, :] * sp.Matrix(coordinates)),
                    repr(nonlinear_op[row]) if nonlinear_op else '')

        nonlinear_constraint = nonlinear_op[row]
        F_args = set(coordinates[0:offset + ss_size]) & \
                 nonlinear_constraint.atoms()
        if linear_op[row, offset:-1].is_zero and not nonlinear_constraint:
            continue

        state_constraint = linear_op[row, offset: offset + ss_size]
        control_constraint = linear_op[row, offset + ss_size:]

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
                    linear_op= linear_op.row_join(sp.SparseMatrix(linear_op.rows, 1, {}))
                    coord = coordinates[offset + ss_size + cv_col]
                    d_coord = sp.Symbol(f"d{str(coord)}")
                    coordinates.append(d_coord)
                    cv_size += 1
                    n += 1

                cv_dict[(0,idx)] = const

        row = row.row_join(sp.SparseMatrix(1, len(added_cvs), cv_dict))
        jac_dx = [nonlinear_constraint.diff(c) for c in coordinates[:ss_size]]
        jac_junciton = [
            nonlinear_constraint.diff(c)
            for c in coordinates[ss_size:offset]
        ]
        jac_x = [
            nonlinear_constraint.diff(c)
            for c in coordinates[offset:
            offset+ss_size]
        ]
        jac_cv = [
            nonlinear_constraint.diff(c)
            for c in coordinates[offset + ss_size:]
        ]

        nlin_row = sp.S(0)

        if any(x!=0 for x in jac_dx):
            logger.warning("Second order constriants not implemented: %s",
                           jac_dx)
        elif any(x!=0 for x in jac_junciton):
            logger.warning("First order junciton constriants not implemented: %s",
                           jac_cv)
        elif any(x!=0 for x in jac_cv):
            logger.warning("First order control constriants not implemented: %s",
                           jac_cv)
        elif any(x!=0 for x in jac_x):
            logger.info("First order constriants: %s", jac_x)
            fx = sum(x*y for x,y in zip(jac_x, coordinates[:ss_size]))
            logger.info(repr(fx))
            p, q = sp.fraction(sp.simplify(fx))
            if row.is_zero:
                lin_dict, nlin = extract_coefficients(
                    p, {c:i for i,c in enumerate(coordinates)},
                    coordinates)

                for k, v in lin_dict.items():
                    row[0, k] += v

                nlin_row += nlin

            else:
                nlin_row += fx

        nonlinear_op = nonlinear_op.col_join(sp.SparseMatrix(1,1,[nlin_row]))

        linear_op = linear_op.col_join(row)
        rows_added += 1

    if rows_added:
        linear_op, nonlinear_op, constraints = \
            smith_normal_form(linear_op, nonlinear_op)

    return coordinates, linear_op, nonlinear_op, constraints


def flatten(sequence):
    for item in sequence:
        if isinstance(item, (list, tuple)):
            for subitem in flatten(item):
                yield subitem
        else:
            yield item


def augmented_rref(matrix, augment=0):

    pivot = 0
    m = matrix.cols - augment
    for col in range(m):
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

        a = matrix[pivot, col]

        for i in range(matrix.rows):
            if i != pivot and matrix[i, col] != 0:
                b = matrix[i, col]/a
                matrix[i, :] += - b * matrix[pivot, :]

        matrix[pivot, :] *= 1 / a

        pivot += 1

        if pivot >= matrix.rows:
            break
    return matrix


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
        M = matrix.row_join(augment)
        k = augment.cols
    else:
        M = matrix
        k = 0
    m, n = M.shape
    M = augmented_rref(M, k)

    Mp = sp.MutableSparseMatrix(n-k, n, {})

    constraints = []
    for row in range(m):
        leading_coeff = -1
        for col in range(row, n-k):
            if M[row, col] != 0:
                leading_coeff = col
                break
        if leading_coeff < 0:
            if not M[row, n-k:].is_zero:
                constraints.append(sum(M[row,:]))
        else:
            Mp[leading_coeff, :] = M[row, :]

    if augment:
        return Mp[:,:-k], Mp[:, -k:], constraints
    else:
        return Mp, sp.SparseMatrix(m,k,{}), constraints


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

    return (inverse_tm, inverse_js, inverse_cm), coordinates
