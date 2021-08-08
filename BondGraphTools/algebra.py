"""
This module contains methods for performing symbolic model reduction.
"""


import logging
import sympy

from .exceptions import SymbolicException

logger = logging.getLogger(__name__)

__all__ = [
    "extract_coefficients",
    "reduce_model",
    "flatten",
    "smith_normal_form",
    "augmented_rref"
]


def extract_coefficients(equation: sympy.Expr,
                         local_map: dict,
                         global_coords: list) -> tuple:
    """

    Args:
        equation: The equation in local coordinates.
        local_map: The mapping from local coordinates to the index of a
            global coordinate.
        global_coords: The list of global co-ordinates.

    Returns:
        The linear and nonlinear parts of the equation in the global
        co-ordinate system.

    Extracts the coordinates from the given equation and maps them into
    the global coordinate space.
    Equations are assumed to come in as sympy expressions of the form
    :math:`\\Phi(x) = 0`.
    local_map is a dictionary mappings

    .. math::

       M: \\rightarrow i

    where :math:`x` are the local co-ordinates and the keys of local_map, and
    the values are the indices :math:`i` such that `global_coord[i]` is the
    corresponding global coordinate. The result is :math:`L,N` such that:

    .. math::

       Ly + N(y) = 0

    """

    coefficients = {}
    nonlinear_terms = sympy.S(0)
    subs = [(k, global_coords[v]) for k, v in local_map.items()]
    local_variables = set(local_map.keys())

    subs.sort(key=lambda x: str(x[1])[-1], reverse=True)
    logger.debug("Extracting coefficients from %s", repr(equation))
    logger.debug("Using local-to-global substitutions %s", repr(subs))

    remainder = equation
    mappings = list(local_map.items())

    while mappings and remainder:
        variable, index = mappings.pop()
        coefficient = equation.coeff(variable)
        if not coefficient:
            continue
        remainder -= coefficient * variable
        if coefficient.atoms() & local_variables:
            nonlinear_terms += (coefficient * variable).subs(subs)
        else:
            coefficients[index] = coefficient
    nonlinear_terms = sympy.expand(nonlinear_terms + remainder.subs(subs))

    logger.debug("Linear terms: %s", repr(coefficients))
    logger.debug("Nonlinear terms: %s", repr(nonlinear_terms))
    return coefficients, nonlinear_terms


def _generate_substitutions(linear_op, nonlinear_op,
                            constraints, coords, size_tup):

    # Lx + F(x) = 0 =>  Ix = (I - L)x - F(x) = Rx - F(x)
    # Since L is in smith normal form (rref and square)
    # If (Rx)_{ii} = 0, and F(x)_i doesn't depend upon x_i
    # then we have x_i = (Rx)_i - F_i(x)
    c_atoms = set(coords)
    atoms = nonlinear_op.atoms() & c_atoms
    for constraint in constraints:
        atoms |= (constraint.atoms() & c_atoms)

    if not atoms:
        logger.debug("No substitutions required")
        return []

    Rx = (sympy.eye(linear_op.rows) - linear_op)

    state_size, network_size, inout_size, n = size_tup
    substitutions = []
    coords_vect = sympy.Matrix(coords)
    for i in reversed(range(2 * (state_size + network_size))):
        co = coords[i]
        if Rx[i, i] == 0 and co in atoms and co not in nonlinear_op[i].atoms():

            eqn = (Rx[i, :] * coords_vect)[0] - nonlinear_op[i]
            pair = (coords[i], eqn)
            logger.debug("Generating substitution %s = %s",
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
    state_size, network_size, inout_size, n = size_tup
    offset = 2 * network_size + state_size

    coord_atoms = set(coordinates[0:offset + state_size])
    linear_op = linear_op[:offset + state_size, :]
    nonlinear_op = nonlinear_op[:offset + state_size, :]

    while constraints:
        constraint, _ = sympy.fraction(constraints.pop())
        logger.debug("Processing constraint: %s", repr(constraint))
        atoms = constraint.atoms() & set(coord_atoms)

        if len(atoms) == 1:
            c = atoms.pop()
            logger.debug("Attempting to find inverse")
            solns = list(sympy.solveset(constraint, c))

            if len(solns) == 1:
                idx = coordinates.index(c)
                sol = solns.pop()

                linear_op = linear_op.col_join(
                    sympy.SparseMatrix(1, linear_op.cols, {(0, idx): 1})
                )
                nonlinear_op = nonlinear_op.col_join(
                    sympy.SparseMatrix(1, 1, {(0, 0): -sol})
                )
                constraint = c - sol
        else:
            logger.warning("..skipping %s", repr(constraint))
            initial_constraints.append(constraint)
        try:
            partials = [constraint.diff(c) for c in coordinates]
        except Exception as ex:
            logger.exception("Could not differentiate %s with respect to %s",
                             repr(constraint), repr(coordinates)
                             )
            raise ex

        if any(p != 0 for p in partials[0:offset]):
            logger.warning("Cannot yet reduce order of %s", repr(constraint))
            initial_constraints.append(constraint)
        else:
            ss_derivs = partials[offset: offset + state_size]
            cv_derivs = partials[offset + state_size:]
            factor = 0
            lin_dict = {}
            nlin = 0
            for idx, coeff in enumerate(ss_derivs):
                if factor == 0 and coeff != 0:
                    factor = 1 / coeff
                    lin_dict.update({(0, idx): 1})
                elif factor != 0 and coeff != 0:
                    new_coeff = sympy.simplify(coeff / factor)
                    if new_coeff.is_number:
                        lin_dict.update({(0, idx): new_coeff})
                    else:
                        nlin += new_coeff * coordinates[idx]
            for idx, coeff in enumerate(cv_derivs):
                if coeff != 0:
                    cv = coordinates[offset + state_size + idx]
                    dvc = sympy.Symbol(f"d{str(cv)}")
                    try:
                        dc_idx = coordinates.index(dvc)
                    except ValueError:
                        dc_idx = len(coordinates)
                        coordinates.append(dvc)
                        inout_size += 1
                        n += 1
                        linear_op = linear_op.row_join(
                            sympy.SparseMatrix(linear_op.rows, 1, {})
                        )
                    eqn = coeff / factor
                    if eqn.is_number:
                        lin_dict.update({(0, dc_idx): eqn})
                    else:
                        nlin += eqn * dvc
            linear_op = linear_op.col_join(
                sympy.SparseMatrix(1, linear_op.cols, lin_dict)
            )
            nonlinear_op = nonlinear_op.col_join(
                sympy.SparseMatrix(1, 1, {(0, 0): nlin})
            )

    linear_op, nonlinear_op, new_constraints = smith_normal_form(
        matrix=linear_op,
        augment=nonlinear_op)

    return linear_op, nonlinear_op, new_constraints + initial_constraints, \
        coordinates, (state_size, network_size, inout_size, n)


def _generate_cv_substitutions(subs_pairs, mappins, coords):
    state_map, port_map, control_map = mappins
    state_size = len(state_map)

    cv_offset = 2 * (state_size + len(port_map))

    control_vars = {str(c) for c in coords[cv_offset:]}
    subs = []
    for var, fx_str in subs_pairs.items():

        if var in control_vars:
            u = sympy.S(var)
        elif var in control_map:
            u = sympy.S(f"u_{control_map[var]}")
        else:
            raise SymbolicException("Could not substitute control variable %s",
                                    str(var))
        fx = sympy.sympify(fx_str)

        subs.append((u, fx))

    return subs


def reduce_model(linear_op, nonlinear_op, coordinates, size_tuple,
                 control_vars=None):
    """
    Simplifies the given system equation.

    Args:
        linear_op: Linear part of the constitutive relations.
        nonlinear_op: The corresponding nonlinear part; a symbolic vector with
        the same number of rows.
        coordinates: a list of all the relevant co-ordinates
        size_tuple:
        control_vars:

    Returns: a tuple describing the reduced system.

    The output of the reduced system is of the form :math:`(x, L, N, G)`
    such that the system dynamics satisfies

    .. math::

        Lx + N(x) = 0
        G(x) = 0
    """

    linear_op, nonlinear_op, constraints = smith_normal_form(
        matrix=linear_op,
        augment=nonlinear_op)

    rows_added = 0
    added_cvs = []
    cv_diff_dict = {}
    lin_dict = {}

    logger.debug("Handling algebraic constraints")

    ###
    # First; take care of control variables
    #

    #
    # Then substitute as much of the junction space as possible.
    #

    subs_list = _generate_substitutions(
        linear_op, nonlinear_op, constraints, coordinates, size_tuple
    )
    logger.debug("Applying substitutions")

    nonlinear_op = nonlinear_op.subs(subs_list)
    constraints = [c.subs(subs_list) for c in constraints]

    logger.debug("Reducing purely algebraic constraints")
    # second, reduce the order of all nonlinear constraints
    linear_op, nonlinear_op, constraints, coordinates, size_tuple =\
        _process_constraints(linear_op, nonlinear_op,
                             constraints, coordinates, size_tuple)
    logger.debug("Applying substitutions, round 2")
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
    # ## New Code
    state_size, network_size, inout_size, n = size_tuple
    offset = 2 * network_size + state_size
    for row in reversed(range(linear_op.rows, offset)):
        atoms = nonlinear_op[row].atoms()
        if not atoms & set(coordinates) and linear_op[row].nnz() > 1:
            logger.debug("Linear constraint in row %s", repr(row))
            for idx in range(state_size):
                v = linear_op[row, idx + offset]
                if v:
                    lin_dict.update({(rows_added, idx): v})
            for idx in range(inout_size):
                v = linear_op[row, idx + offset + state_size]
                if v:
                    cv_diff_dict.update({(rows_added, idx): v})

    for row in range(offset, linear_op.rows):
        logger.debug("Testing row %s: %s + %s", repr(row),
                     repr(linear_op[row, :] * sympy.Matrix(coordinates)),
                     repr(nonlinear_op[row]) if nonlinear_op else '')

        nonlinear_constraint = nonlinear_op[row]
        # F_args = set(coordinates[0:offset + state_size]) & \
        #    nonlinear_constraint.atoms()

        if linear_op[row, offset:-1].is_zero_matrix \
                and not nonlinear_constraint:
            continue

        state_constraint = linear_op[row, offset: offset + state_size]
        control_constraint = linear_op[row, offset + state_size:]

        row = state_constraint.row_join(
            sympy.SparseMatrix(1, offset + inout_size, {}))

        cv_dict = {}
        if not control_constraint.is_zero_matrix:
            logger.debug("Found higher order control constraint")
            for cv_col in range(control_constraint.cols):
                const = control_constraint[cv_col]
                if not const:
                    continue

                try:
                    idx = added_cvs.index(cv_col)
                except ValueError:
                    idx = len(added_cvs)
                    added_cvs.append(cv_col)
                    linear_op = linear_op.row_join(
                        sympy.SparseMatrix(linear_op.rows, 1, {}))
                    coord = coordinates[offset + state_size + cv_col]
                    d_coord = sympy.Symbol(f"d{str(coord)}")
                    coordinates.append(d_coord)
                    inout_size += 1
                    n += 1

                cv_dict[(0, idx)] = const

        row = row.row_join(sympy.SparseMatrix(1, len(added_cvs), cv_dict))
        jac_dx = [nonlinear_constraint.diff(c)
                  for c in coordinates[:state_size]]
        jac_junciton = [
            nonlinear_constraint.diff(c)
            for c in coordinates[state_size:offset]
        ]
        jac_x = [
            nonlinear_constraint.diff(c)
            for c in coordinates[offset:
                                 offset + state_size]
        ]
        jac_cv = [
            nonlinear_constraint.diff(c)
            for c in coordinates[offset + state_size:]
        ]

        nlin_row = sympy.S(0)

        if any(x != 0 for x in jac_dx):
            logger.warning("Second order constraint not implemented: %s",
                           jac_dx)

        elif any(x != 0 for x in jac_junciton):
            logger.warning(
                "First order junciton constraint not implemented: %s",
                str(jac_junciton)
            )

        elif any(x != 0 for x in jac_cv):
            logger.warning(
                "First order control constraint not implemented: %s",
                str(jac_cv)
            )

        elif any(x != 0 for x in jac_x):
            logger.debug("First order constriants: %s", jac_x)
            fx = sum(x * y for x, y in zip(jac_x, coordinates[:state_size]))
            logger.debug(repr(fx))
            p, q = sympy.fraction(sympy.simplify(fx))
            if row.is_zero_matrix:
                lin_dict, nlin = extract_coefficients(
                    p, {c: i for i, c in enumerate(coordinates)},
                    coordinates)

                for k, v in lin_dict.items():
                    row[0, k] += v

                nlin_row += nlin

            else:
                nlin_row += fx

        nonlinear_op = nonlinear_op.col_join(
            sympy.SparseMatrix(1, 1, [nlin_row]))

        linear_op = linear_op.col_join(row)
        rows_added += 1

    if rows_added:
        linear_op, nonlinear_op, constraints = \
            smith_normal_form(linear_op, nonlinear_op)

    return coordinates, linear_op, nonlinear_op, constraints


def flatten(sequence):
    """
    Gets a first visit iterator for the given tree.
    Args:
        sequence: The iterable that is to be flattened

    Returns: iterable
    """
    for item in sequence:
        if isinstance(item, (list, tuple)):
            for subitem in flatten(item):
                yield subitem
        else:
            yield item


def augmented_rref(matrix, augmented_rows=0):
    """ Computes the reduced row-echelon form (rref) of the given augmented
    matrix.

    That is for the augmented  [ A | B ], we fine the reduced row echelon form
    of A.

    Args:
        matrix (sympy.MutableSparseMatrix): The augmented matrix
        augmented_rows (int): The number of rows that have been augmented onto
         the matrix.

    Returns: a matrix M =  [A' | B'] such that A' is in rref.

    """

    m = matrix.cols - augmented_rows
    col = 0
    row = 0
    # forward pass
    while col < m and row < matrix.rows - 1:
        pivot = row
        while pivot < matrix.rows:
            if matrix[pivot, col] != 0:
                break
            pivot += 1

        if pivot == matrix.rows:
            col += 1
            continue    # no entry on this row
        elif pivot != row:
            matrix.row_swap(pivot, row)

        matrix[row, :] = sympy.expand(matrix[row, :] / matrix[row, col])
        for j in range(row + 1, matrix.rows):
            if matrix[j, col] != 0:
                matrix[j, :] = sympy.expand(
                    matrix[j, :] - matrix[j, col] * matrix[row, :])

        col += 1
        row += 1
    # reverse pass
    for row in range(1, matrix.rows):
        col = row
        while col < m:
            if matrix[row, col] != 0:
                break
            col += 1
        if col == m:
            break

        for i in range(row):
            a = matrix[i, col]
            if a != 0:
                matrix[i, :] = sympy.expand(matrix[i, :] - a * matrix[row, :])
    return matrix


def smith_normal_form(matrix, augment=None):
    """Computes the Smith normal form of the given matrix.


    Args:
        matrix:
        augment:

    Returns:
        n x n smith normal form of the matrix.
        Particularly for projection onto the nullspace of M and the orthogonal
        complement that is, for a matrix M,
        P = _smith_normal_form(M) is a projection operator onto the
        nullspace of M

    """

    if augment is not None:
        M = matrix.row_join(augment)
        k = augment.cols
    else:
        M = matrix
        k = 0
    m, n = M.shape
    M = augmented_rref(M, k)

    Mp = sympy.MutableSparseMatrix(n - k, n, {})

    constraints = []
    for row in range(m):
        leading_coeff = -1
        for col in range(row, n - k):
            if M[row, col] != 0:
                leading_coeff = col
                break
        if leading_coeff < 0:
            if not M[row, n - k:].is_zero_matrix:
                constraints.append(sum(M[row, :]))

        else:
            Mp[leading_coeff, :] = M[row, :]

    if augment:
        return Mp[:, :-k], Mp[:, -k:], constraints
    else:
        return Mp, sympy.SparseMatrix(m, k, {}), constraints


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


def get_relations_iterator(component, mappings, coordinates, io_map=None):
    local_tm, local_js, local_cv = component.basis_vectors
    inv_tm, inv_js, inv_cv = mappings

    num_ports = len(inv_js)
    num_state_vars = len(inv_tm)
    local_map = {}

    # todo: Fix this dirty hack; there has to be a better way to hand io ports
    for cv, value in local_cv.items():
        try:
            local_map[cv] = 2 * (num_ports + num_state_vars) + inv_cv[value]
        except KeyError:
            logger.debug("Could not find %s, trying the io_map", value)
            key = io_map[value]
            local_map[cv] = key
            logger.debug("Mapping %s to co-ord %s", cv, coordinates[key])

    for (x, dx), coord in local_tm.items():
        local_map[dx] = inv_tm[coord]
        local_map[x] = inv_tm[coord] + num_state_vars + 2 * num_ports

    for (e, f), port in local_js.items():
        local_map[e] = 2 * inv_js[port] + num_state_vars
        local_map[f] = 2 * inv_js[port] + num_state_vars + 1
    logger.debug("Getting relations iterator for %s", repr(component))
    for relation in component.constitutive_relations:
        if relation:
            yield extract_coefficients(relation, local_map, coordinates)
        else:
            yield {}, 0.0
