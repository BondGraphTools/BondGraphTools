from collections import namedtuple
from itertools import permutations

import numpy as np
from matplotlib.offsetbox import AnchoredText
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle
import matplotlib.pyplot as plt

from scipy.sparse import dok_matrix
from scipy.sparse.csgraph import floyd_warshall

from .exceptions import InvalidComponentException

FONT = 14
FONT_SM = 10


def _build_graph(system):

    try:
        comp_map = {comp: i for i, comp in enumerate(system.components.values())}
        graph = dok_matrix((len(comp_map), len(comp_map)), dtype=np.int)
        for (c1, _), (c2, _) in system.bonds:
            graph[(comp_map[c1], comp_map[c2])] = 1
            graph[(comp_map[c2], comp_map[c1])] = 1

    except AttributeError as ex:
        raise InvalidComponentException(
            "Invalid System: has no components"
        ) from ex

    return graph.tocsr(copy=False)


def _metro_layout(graph):

    n, _ = graph.shape

    D = floyd_warshall(graph, directed=False)

    ecen = D.max(1)
    eta_0 = int(ecen.min())
    level_struct = {int(eta): np.where(ecen == eta)[0].tolist() for eta in
                    np.unique(ecen)}
    z = np.zeros((n, 1), dtype=np.complex_)
    z, R, N = initialise(z, eta_0, level_struct[eta_0], D)

    for k, eta in enumerate(range(eta_0, max(level_struct.keys()))):
        R += 1
        z, N = optimise_layer(z, N, R, eta, level_struct, graph, D)

    Z = z - z.T
    zmin = np.abs(Z[Z != 0]).min()

    return [(round(zp.real / zmin), round(zp.imag / zmin)) for zp in z.flatten()]


def permute(z, N, r, indicies):
    i_len = len(indicies)
    for idx in permutations(range(N), i_len):
        zp = z.copy()
        zp[indicies, 0] = [r * np.exp(w * 2 * np.pi * np.complex(0, 1) / N) for
                           w in idx]
        yield zp


def initialise(z, eta_0, level_set, d):
    if len(level_set) == 1:
        return z, 0, 0
    elif len(level_set) < 4:
        R_0 = 0.5
    else:
        R_0 = 1

    N_0 = int(2 * np.ceil(len(level_set) / 2))
    f = lambda x: init_obj(x, level_set, d, eta_0)

    for zp in permute(z, N_0, R_0, level_set):
        phi = f(zp)
        try:
            if phi_min > phi:
                phi_min = phi
                z_min = zp
        except NameError:
            z_min = zp
            phi_min = phi

    return z_min, R_0, N_0


def init_obj(z, level_set, d, eta_0):
    phi = sum(
        [(np.abs(z[i] - z[j]) / d[i, j] - 1) ** 2 for i in level_set for j in
         level_set if i != j])
    return phi


def optimise_layer(z, N_eta, R, eta, level_struct, adj, dist):
    N = int(max(N_eta, 2 * np.ceil(len(level_struct[eta + 1]) / 2)))

    # f = lambda x: objective(x, eta, level_struct, adj)
    f = lambda x: objective_f(x, eta, level_struct, dist)
    z_min = z
    for zp in permute(z, N, R, level_struct[eta + 1]):
        phi = f(zp)
        try:
            if phi_min > phi:
                z_min = zp
                phi_min = phi
        except NameError:
            z_min = zp
            phi_min = phi

    return z_min, N


def objective_f(z, eta, level_struct, d):
    target_set = [i for i in level_struct[eta] + level_struct[eta + 1]]

    phi = sum(
        [(np.abs(z[i] - z[j]) / d[i, j] - 1) ** 2 for i in level_struct[eta + 1]
         for j in target_set if i != j])
    return phi


class PortGlyph(Text):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = 0.075
        self.height = 0.075
        self.x = 0
        self.y = 0
    @property
    def pos(self):
        return (self.x, self.y)


class Glyph:
    def __init__(self, node):
        self._node = node
        self._axes = None
        self.x = 0
        self.y = 0
        self.string = None
        self.width = 0.1
        self.height = 0.1
        self._text = None
        self.ports = {
            'top':[],
            'right': [],
            'bottom': [],
            'left':[]
        }

    @property
    def pos(self):
        return self.x, self.y

    @pos.setter
    def pos(self, value):
        self.x, self.y = value

    @property
    def axes(self):
        return self._axes

    @axes.setter
    def axes(self, ax):
        self._axes = ax
        self._text = Text(
            self.x,
            self.y,
            f"${self.string}$",
            horizontalalignment='center',
            verticalalignment='center',
            size=FONT,
            usetex=True)

        ax.add_artist(self._text)

    def add_port(self, string, dir):

        dx,dy = dir
        text_dict = {
            's':string,
            'xy':self.pos,
            'size': FONT_SM
        }

        if dy > abs(dx):
            text_dict.update({
                'xytext':(self.x, self.y+self.height/2),
                'horizontalalignment':'center',
                'verticalalignment': 'bottom'
                }
            )
        elif -dy > abs(dx):
            text_dict.update({
                'xytext':(self.x, self.y-self.height/2),
                'horizontalalignment':'center',
                'verticalalignment': 'top'
                }
            )

        elif dx >= abs(dy):
            text_dict.update({
                'xytext':(self.x+self.width/2, self.y),
                'horizontalalignment': 'left',
                'verticalalignment': 'center'
                }
            )

        else:
            text_dict.update({
                'xytext': (self.x - self.width / 2, self.y),
                'horizontalalignment': 'right',
                'verticalalignment': 'center'
            }
            )
        text = Text(text_dict)
        self.axes.add_artist(text)
        port = PortGlyph()
        port.text = text

        return self




class Bond(Line2D):

    def __init__(self, port_1, port_2, *args, **kwargs):
        self.port_1 = port_1
        self.port_2 = port_2
        super().__init__([],[],*args, **kwargs)

    def calc_lines(self):
        x1, y1 = self.port_1.pos
        x2, y2 = self.port_2.pos

        r1 = max(self.port_1.height, self.port_1.width)
        r2 = max(self.port_2.height, self.port_2.width)

        dx = x2 - x1
        dy = y2 - y1
        x1 += r1*dx
        y1 += r1*dy
        x2 -= r2*dx
        y2 -= r2 * dy
        self.set_xdata([x1, x2])
        self.set_ydata([y1, y2])


class GraphLayout(Glyph):
    def draw(self, layout=_metro_layout):

        fig = plt.figure(
            figsize=(12, 9), dpi=80
        )

        ax = fig.gca()

        graph = _build_graph(self._node)

        points = layout(graph)
        bonds = []
        x_min = 0
        x_max = 0
        y_min = 0
        y_max = 0
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)

        for (local_name, component), (x,y) in zip(self._node.components.items(),
                                                  points):
            x_min = min(x, x_min)
            x_max = max(x, x_max)
            y_min = min(y, y_min)
            y_max = max(y, y_max)

            component.view.pos = (x,y)

            component.view.string = local_name
            component.view.axes = ax

        for (c1, port_1), (c2, port_2) in self._node.bonds:

            label_1 = f"[{port_1}]"
            label_2 = f"[{port_2}]"

            dx = c2.view.x - c1.view.x
            dy = c2.view.y - c1.view.y

            p1 = c1.view.add_port(label_1, (dx, dy))
            p2 = c2.view.add_port(label_2, (-dx, -dy))
            bond = Bond(p1, p2)
            ax.add_artist(bond)
            bonds.append(bond)

        for bond in bonds:
            bond.calc_lines()

        ax.axis([x_min - 1, x_max + 1, y_min - 1, y_max + 1])
        return fig


def find_renderer(fig):

    if hasattr(fig.canvas, "get_renderer"):
        #Some backends, such as TkAgg, have the get_renderer method, which
        #makes this easy.
        renderer = fig.canvas.get_renderer()
    else:
        #Other backends do not have the get_renderer method, so we have a work
        #around to find the renderer.  Print the figure to a temporary file
        #object, and then grab the renderer that was used.
        #(I stole this trick from the matplotlib backend_bases.py
        #print_figure() method.)
        import io
        fig.canvas.print_pdf(io.BytesIO())
        renderer = fig._cachedRenderer
    return(renderer)