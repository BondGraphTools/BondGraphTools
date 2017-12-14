import math
import random

def _optimise(nodes, edges, N=10000):
    """
    e is a list of edge tuples [(i,j), (i,m), ... ]


    Returns:

    """
    n = 0
    obj = 1000
    targets = set()
    x, y = zip(*nodes)
    T = 1
    xt = [p for p in x]
    yt = [p for p in y]
    lenx = len(x)
    step = math.ceil(math.sqrt(lenx))
    while n < N:
        # generate new candidate
        for i in range(lenx):
            if i in targets:
                xt[i] = x[i] + random.randint(-step, step)
                yt[i] = y[i] + random.randint(-step, step)
            else:
                xt[i] = x[i] + random.randint(-step, step)
                yt[i] = y[i] + random.randint(-step, step)

        # make sure we're not sitting on top of eachother
        if any(xt[i] == xt[j] and yt[i] == yt[j]
                for i in range(lenx) for j in range(lenx) if i != j):
            continue

        targets = set()
        crossings = 0
        dist = 0

        bd_x = [xt[j] - xt[i] for i, j in edges]
        bd_y = [yt[j] - yt[i] for i, j in edges]

        for k, (i, j) in enumerate(edges):
            bx = bd_x[k]
            by = bd_y[k]

            dist += bx**2 + by**2

            for l in range(k, len(edges)):
                ip, jp = edges[l]

                rs = bd_x[k]*bd_y[l] - bd_x[l]*bd_y[k]
                if rs != 0:
                    t = ((x[ip] - x[i])*bd_y[l] - (y[ip] - y[i])*bd_x[l])/rs
                    u = ((x[ip] - x[i])*bx - (y[ip] - y[i])*by)/rs
                    if 0 <= u <= 1 and 0 <= t <= 1:
                        crossings += 1
                        targets &= {i, j, ip, jp}

        new_obj = dist*(1 + crossings)**2
        delta = new_obj - obj
        if delta <= 0 or math.exp(-delta/T) > random.uniform(0, 1):
            x = [p for p in xt]
            y = [p for p in yt]
            obj = new_obj
        n += 1

    xm = int(sum(x)/lenx)
    ym = int(sum(y)/lenx)
    out = [(xp - xm, yp - ym) for xp, yp in zip(x,y)]

    return out

def arrange(bond_graph):
    nodes = list()
    for node in bond_graph.nodes.values():
        try:
            x, y = node.pos
            nodes.append((x,y))
        except:
            nodes.append((0, 0))
    edges = set()
    for i, j,_,_ in bond_graph.bonds:
        if i < j:
            edges.add((i, j))
        else:
            edges.add((j, i))
    pos = _optimise(nodes, list(edges))
    for n, node in enumerate(bond_graph.nodes.values()):
        node.pos = pos[n]




