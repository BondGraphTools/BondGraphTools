import sympy as sp

def extract_coefficients(equation, constants=None):

    coeff_dict = {}

    for term in equation.expand().args:

        products = flatten(term.as_coeff_mul())

        coeff = sp.S(1)
        base = sp.S(1)

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