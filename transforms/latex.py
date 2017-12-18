import sympy as sp
from sympy.physics.mechanics import dynamicsymbols
import numpy as np



def bondgraph_to_sympy(graph):

    xdot, x, p, relations, labels = _construct_coordinates(graph)
    t = sp.symbols("t")
    Y = x + xdot
    ydot = xdot + [sp.diff(X, t) for X in xdot]

    N = len(Y)
    zero = sp.zeros(N,N)
    h_zero = sp.zeros(N/2, N/2)

    #L = -sp.eye(len(x)).col_insert(0,-sp.eye(len(x)))
    L = sp.Matrix(0, len(Y), [])
    NL = sp.Matrix(0, 1, [])

    for rel in relations:
        H = sp.hessian(rel, Y)
        if H == zero:
            row = sp.Matrix(1, N, [sp.diff(rel, y) for y in Y])
            L = L.row_insert(-1, row)
        else:
            if H[0:N/2, 0:N/2] == h_zero:
                row = sp.Matrix(1, 1, [rel])
            else:
                df = sp.Matrix([sp.diff(rel, y) for y in Y])
                row = df.dot(ydot)
            NL = NL.row_insert(-1, row)

    return L

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
        for port_id in node.ports:
            en = dynamicsymbols("e_{}".format(n))
            fn = dynamicsymbols("f_{}".format(n))
            pn = dynamicsymbols("p_{}".format(n))
            qn = dynamicsymbols("q_{}".format(n))
            coord_labels[(node_id, port_id)] = (node.name, [en,fn, pn,qn])
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
            for param in node.local_params:
                params_labels[an] = (node.name, param)
                p_an = sp.symbols("a_{}".format(an))
                A.append(p_an)
                subs.append((param, p_an))
                an += 1

        if node.global_params:
            for param in node.global_params:
                try:
                    p_bn = globals_map[param]
                except KeyError:
                    p_bn = sp.symbols("b_{}".format(bn))
                    globals_map[param] = p_bn
                    globals_labels[bn] = param
                    bn += 1

                subs.append((param, p_bn))

        node_cr = _build_relations(node)
        relations += [r.subs(subs) for r in node_cr]

    subs, (E, F, P, Q) = _construct_edge_subs(graph, coord_labels, (E, F, P, Q))

    relations = [r.subs(subs).simplify() for r in relations]
    print(E)
    B = list(globals_map.values())
    return E + F, P + Q, A + B, relations, (coord_labels, params_labels, globals_labels)



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
            (n1, p1), (n2, p2) = edge
            edge_list.remove(edge)

            if n1 == from_node:
                _,  (e1, f1, p1, q1) = coord_labels[(n1, p1)]
                _,  (e2, f2, p2, q2) = coord_labels[(n2, p2)]
                dirn = -1  # power leaving
            else:
                _, (e1, f1, p1, q1) = coord_labels[(n2, p2)]
                _, (e2, f2, p2, q2) = coord_labels[(n1, p1)]
                dirn = 1  # power coming in

            subs += [
                (e1, e2), (f1, -f2), (p1, p2), (q1, -q2)
            ]
            E.remove(e1)
            F.remove(f1)
            P.remove(p1)
            Q.remove(q1)

    return subs, (E, F, P, Q)