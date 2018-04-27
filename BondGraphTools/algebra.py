import sympy as sp


def extract_coefficients(equation, local_map, global_coords):

    coeff_dict = {}
    nonlinear_terms = sp.S(0)
    subs = [(k, global_coords[v]) for k,v in local_map.items()]

    for term in equation.expand().args:
        prod_iter = flatten(term.as_coeff_mul())
        coeff = next(prod_iter)
        base = list(prod_iter)

        if not base:
            coeff_dict[len(global_coords)] = coeff
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

