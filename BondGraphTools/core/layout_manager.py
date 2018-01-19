import math
import random
import logging
from functools import reduce
import numpy as np
from scipy.sparse.csgraph import floyd_warshall, laplacian
from collections import defaultdict

logger = logging.getLogger(__name__)
DIRS = [(-1, 0), (0, -1), (1, 0), (0, 1)]


def simulated_annealing(nodes, edges, N=10000):
    """
    e is a list of edge tuples [(i,j), (i,m), ... ]


    Returns:

    """
    dirs = [(a,b) for a in range(-1,2) for b in range(-1,2)]

    n = 0
    obj = 1000
    targets = set()
    w1 = 0.5
    w2 = 0.01

    x, y = zip(*nodes)
    lenx = len(x)
    step = math.ceil(math.sqrt(lenx))
    x = [int(step*math.cos(2*math.pi*k/lenx)) for k in range(lenx)]
    y = [int(step*math.sin(2*math.pi*k/lenx)) for k in range(lenx)]
    T = 1

    xt = [p for p in x]
    yt = [p for p in y]

    while n < N:
        # generate new candidate
        dist = 1 + step*int(math.exp(-n))
        for i in range(lenx):
            if i in targets:
                xt[i] = x[i] + random.randint(-step, step)
                yt[i] = y[i] + random.randint(-step, step)
            else:
                d = random.randint(0, 8)
                dx, dy = dirs[d]
                xt[i] = x[i] + dx*dist
                yt[i] = y[i] + dy*dist

        targets = set()
        crossings = 0
        zero_bonds = 0
        dist = 0
        node_dist = sum(
            sum(((x1-x2)**2 + (y1 - y2)**2 + w2)**(-2) for x1, y1 in zip(xt,yt))
                for x2,y2 in zip(xt,yt))

        bd_x = [xt[j] - xt[i] for i, j in edges]
        bd_y = [yt[j] - yt[i] for i, j in edges]

        for k, (i, j) in enumerate(edges):
            bx = bd_x[k]
            by = bd_y[k]
            d = bx**2 + by**2
            if d == 0:
                zero_bonds += 1
                targets &= {i, j}
            else:
                dist += d
                for l in range(k, len(edges)):
                    ip, jp = edges[l]

                    rs = bd_x[k]*bd_y[l] - bd_x[l]*bd_y[k]
                    if rs != 0:
                        t = ((x[ip] - x[i])*bd_y[l] - (y[ip] - y[i])*bd_x[l])/rs
                        u = ((x[ip] - x[i])*bx - (y[ip] - y[i])*by)/rs
                        if 0 <= u <= 1 and 0 <= t <= 1:
                            crossings += 1
                            targets &= {i, j, ip, jp}

        new_obj = dist**2*(1 + zero_bonds + crossings)**2 + w1 * node_dist

        delta = new_obj - obj
        if delta <= 0 or math.exp(-delta/T) > random.uniform(0, 1) or n == 0:
            x = [p for p in xt]
            y = [p for p in yt]
            obj = new_obj
        n += 1

    xm = int(sum(x)/lenx)
    ym = int(sum(y)/lenx)
    if abs(min(x) - max(x)) < abs(min(y) - max(y)):
        out = [(ym - yp, xp - xm) for xp, yp in zip(x, y)]
    else:
        out = [(xp - xm, yp - ym) for xp, yp in zip(x,y)]

    return out


def _make_planar_graph(bond_graph):
    nodes = list()
    adj_dict = defaultdict(lambda: dict())
    for node in bond_graph.nodes.values():
        try:
            x, y = node.pos
            nodes.append((x, y))
        except TypeError:
            nodes.append((0, 0))
    adj_matrix = np.zeros((len(nodes), len(nodes)))
    edges = []
    for k, (i, j, _, _) in enumerate(bond_graph.bonds):
        edges.append((i,j))
        adj_matrix[i,j] = 1
        adj_matrix[j,i] = 1

        adj_dict[i][j] = k
        adj_dict[j][i] = k

    return nodes, edges


def branch_and_bound(nodes, edges):

    n = len(nodes)
    A = np.zeros((n, n))

    for (i, j) in edges:
        A[i, j] = 1

    d = np.triu(
        floyd_warshall(A, directed=False, unweighted=True, overwrite=False)
    )

    M = np.triu(_distance_matrix(nodes))
    W = np.zeros((n, n))
    W[d != 0] = d[d != 0]**(-2)

    phi_0 = np.sum(W * (M - d)**2) + np.sum(d[M == 0])

    X0 = tuple((x, y) for x, y in nodes)
    tree = [(X0, phi_0)]
    new_list = []
    i = 0
    max_i = 10000

    seen_list = set(X0)
    while tree and i < max_i:
        X, phi = tree.pop()

        for j in range(n):
            for dx, dy in DIRS:
                Xjk = tuple((x + dx, y + dy) if l == j else (x, y)
                            for l, (x, y) in enumerate(X))

                if Xjk not in seen_list:
                    seen_list.add(Xjk)
                    M_jk = np.triu(_distance_matrix(Xjk))
                    phi_jk = np.sum(W * (M_jk - d) ** 2)
                    phi_jk += np.sum(
                        d[M_jk == 0]
                    )

                    new_list.append((Xjk, phi_jk))

        new_list.append((X, phi))
        new_list.sort(key=lambda x: x[1], reverse=True)
        tree += new_list
        new_list = []
        i += 1
        # if sub_tree[0][1] > phi:
        #     if len(leaves) > max_leaves:
        #         leaves.pop(-1)
        #     leaves.add((X, phi))
        # else:
        #     sub_tree.add((X, phi))
        #
        # tree |= sub_tree
        #
        # index = max_tree - len(tree)
        # if index < -1:
        #     print("clearing")
        #     del tree[index:-1]
    pos, _ = tree[-1]

    return pos


def _distance_matrix(nodes):

    n = len(nodes)
    M = np.zeros((n, n))
    for i, (x1, y1) in enumerate(nodes):
        for j in range(i, n):
            x2, y2 = nodes[j]
            #M[i, j] = max(abs(x1 - x2), abs(y1 - y2))
            M[i, j] = ((x1 - x2)**2 + (y1 - y2)**2)**0.5
    return M


def contract_graph(nodes, edges):

    deg = defaultdict(int)
    new_edges = []
    A = np.zeros((len(nodes), len(nodes)))
    for i, j in edges:
        deg[i] += 1
        deg[j] += 1
        A[i, j] = 1
        A[j, i] = 1
        new_edges.append(((i, j), 0))

    new_edges = _recurse_contract(deg, new_edges)
    inverse_map = {}
    forward_map = {}
    reduced_edges = []
    reduced_nodes = []

    for (i, j), weight in new_edges:
        for k in (i, j):
            if k not in forward_map:
                n = len(forward_map)
                forward_map[k] = n
                inverse_map[n] = k
                reduced_nodes.append(nodes[k])

        n_i = forward_map[i]
        n_j = forward_map[j]
        reduced_edges.append(((n_i, n_j), 2*weight + 1))

    n = len(forward_map)

    W_Adj = np.zeros((n, n))

    for (i, j), w in reduced_edges:
        if i < j:
            W_Adj[i, j] = w
        else:
            W_Adj[j, i] = w

    reduced_nodes = _layout_reduced_graph(reduced_nodes, W_Adj)
    return reduced_nodes
    # pre_interpolated_nodes = [None for _ in nodes]
    # for n_i, X in enumerate(reduced_nodes):
    #     i = inverse_map[n_i]
    #     pre_interpolated_nodes[i] = X
    #
    # return _interpolate(pre_interpolated_nodes
    #
    #                     )


def _interpolate(nodes, edges):

    return


def force_directed(nodes, edges):

    W_adj = np.zeros((len(nodes), len(nodes)))
    for i, j in edges:
        W_adj[i, j] = 1
        W_adj[j, i] = 1

    graph_dist = floyd_warshall(W_adj, directed=False)
    L = laplacian(W_adj)
    n = len(nodes)
    r0 = n**2
    Xnew = [
        (r0/L[i, i]*np.cos(2*np.pi*i/(n-1)),
         r0/L[i, i]*np.sin(2*np.pi*i/(n-1))) for i in range(n)
    ]
    euclid_dist = _distance_matrix(Xnew)
    weights = np.zeros(W_adj.shape)
    weights[graph_dist!=0] = graph_dist[graph_dist!=0]**(-2)
    wp = np.zeros(W_adj.shape)
    wp[W_adj > 1] = W_adj[W_adj > 1]
    wp += wp.transpose()
    wp2 = wp
    wp2[wp!=0] = 0.5
    EP = np.ones_like(W_adj)
    EYE = np.eye(n)
    sigma_new = np.sum(weights * (euclid_dist - graph_dist) ** 2
                           + wp2*(euclid_dist - wp) ** 2
                           + (graph_dist / (euclid_dist + EP) - EYE)
                           )
    sigma_old = sigma_new*10
    eps = 0.00001
    delta = 0.01
    its = 0
    max_its = 2500
    scale = 0
    while max_its > its and (scale == 0 or
                             abs(sigma_new - sigma_old) > eps*sigma_old):
        Xold = Xnew
        sigma_old = sigma_new
        Xnew = []
        coeff = weights * (euclid_dist - graph_dist) + wp2 *(euclid_dist - wp)
        for k in range(n):
            dx = 0
            dy = 0
            for i in range(n-1):
                for j in range(i, n):
                    xi, yi = Xold[i]
                    xj, yj = Xold[j]
                    rij = max(euclid_dist[i, j], eps)

                    if i == k:
                        dx += (coeff[i, j] - graph_dist[i,j]*rij**(-2)) * (xi - xj)/rij
                        dy += (coeff[i, j] - graph_dist[i,j]*rij**(-2)) * (yi - yj)/rij
                    elif j == k:
                        dx -= (coeff[i, j] - graph_dist[i,j]*rij**(-2)) * (xi - xj)/rij
                        dy -= (coeff[i, j] - graph_dist[i,j]*rij**(-2)) * (yi - yj)/rij
            x, y = Xold[k]
            x, y = x - delta*dx, y - delta*dy

            Xnew.append((x,y))

        euclid_dist = _distance_matrix(Xnew)
        scale = (euclid_dist + 2 * np.eye(n)).min()
        sigma_new = np.sum(weights * (euclid_dist - graph_dist) ** 2
                           + wp2 * (euclid_dist - wp) ** 2
                           + (graph_dist / (euclid_dist + EP) - EYE)
                           )
        its += 1

    xm, ym = reduce((lambda P, Q: (P[0]+Q[0], P[1]+Q[1])), Xnew)
    Xnew = [(x - xm/n, y - ym/n) for x,y in Xnew]
    euclid_dist = _distance_matrix(Xnew)
    scale = min(euclid_dist[euclid_dist!= 0])/2
    Xnew =[(np.ceil(x/scale) if x > 0 else np.floor(x),
            np.ceil(y/scale) if y > 0 else np.floor(y))
            for x, y in Xnew]

    return Xnew


def _recurse_contract(vertex_degree, edges):

    new_edges = []
    joins = 0
    while edges:
        (i, j), w = edges.pop()
        if vertex_degree[i] == 1 or vertex_degree[j] == 1:
            pass
        elif vertex_degree[i] > 2 and vertex_degree[j] > 2:
            new_edges.append(((i, j), w))
            # store and skip
        else:
            if vertex_degree[i] == 2:
                pivot = i
                base = j
            else:
                pivot = j
                base = i

            idx = 0
            while idx < len(edges):
                if pivot not in edges[idx][0]:
                    idx += 1
                else:
                    (k, l), w = edges.pop(idx)
                    if k == pivot:
                        new_edges.append(((base, l), w + 1))
                    else:
                        new_edges.append(((k, base), w + 1))
                    break
            joins += 1
    if joins == 0:
        return new_edges
    else:
        return _recurse_contract(vertex_degree, new_edges)


def _layout_algorithm(nodes, edges, constriants=None):

    N = len(nodes)
    K = len(edges)
    E = np.zeros((K, N))
    A = np.np.zeros((N, N))
    x, y = map(lambda q: np.array(q, ndmin=2).transpose(), zip(*nodes))
    # x y are always column vectors

    for k, (i, j) in enumerate(edges):
        A[i, j] = 1
        A[j, i] = 1
        E[k, j] = -1
        E[k, i] = 1

    degree = A.sum(1)
    D = floyd_warshall(A, directed=False)

    # compute Xx Xy; euclidian distances.


def _objective_function(x, y, adj, dist, deg):

    eps_m = 0.01
    n = len(x)
    xbar = x.sum()
    ybar = y.sum()
    Xx = x - x.transpose()
    Yy = y - y.transpose()

    M = (Xx**2 + Yy**2 + eps_m)**(-1)
    PHI = (M - np.eye(n)*eps_m).sum()/2
    Vx = Xx + n*np.diag(x) - np.eye(n)*xbar
    Vy = Yy + n*np.diag(y) - np.eye(n)*ybar
    dx = Vx.dot((M**(-2)).sum(1, keepdims=True))/2
    dy = Vy.dot((M**(-2)).sum(1, keepdims=True))/2

    PHI += ((adj*Xx)**2 + (adj*Yy)**2 - adj).sum()/2
    DB = adj*(1 - adj*(Xx**2 + Yy**2 + eps_m)**(-0.5))
    dx -= 2 * (DB*Xx).sum(1, keepdims=True)
    dy -= 2 * (DB*Yy).sum(1, keepdims=True)


def arrange(bond_graph,
            algorithm=branch_and_bound,
            **kwargs):
    nodes, edges = _make_planar_graph(bond_graph)
    args = []
    nodes = algorithm(nodes, edges, *args, **kwargs)

    for i in bond_graph.nodes:
        bond_graph.nodes[i].pos = nodes[i]





