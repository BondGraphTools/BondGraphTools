"""Tools for visualising bond graph.

This module contains temporary tools for producing visualizations of
bond graph network topologies.

"""

import logging

import numpy as np

from scipy.sparse import dok_matrix
from matplotlib.lines import Line2D
import networkx as nx

from .exceptions import InvalidComponentException

logger = logging.getLogger(__name__)
FONT = 14
FONT_SM = 10

__all__ = ["draw"]


def draw(system):
    """
    Produces a network layout of the system.

    Args:
        system: The system to visualise

    Returns:
        :obj:`matplotlib.Plot`
    """
    import matplotlib.pyplot as plt

    fig = plt.figure(
        figsize=(12, 9), dpi=80
    )
    plt.ioff()
    ax = fig.gca()
    ax.set_aspect("equal")
    ax.set_title(f"{system.name}")
    return system.view.draw(ax)


def _build_graph(system):

    try:
        comp_map = {comp: i for i, comp in enumerate(system.components)}
        graph = dok_matrix((len(comp_map), len(comp_map)), dtype=int)
        for (c1, _), (c2, _) in system.bonds:
            graph[(comp_map[c1], comp_map[c2])] = 1
            graph[(comp_map[c2], comp_map[c1])] = 1

    except AttributeError as ex:
        raise InvalidComponentException(
            "Invalid System: has no components"
        ) from ex

    return graph.tocsr(copy=False)


def _networkx_layout(graph):
    nx_graph = nx.Graph(graph)
    layout = nx.kamada_kawai_layout(nx_graph, scale=20)
    pos = [(pair[0], pair[1]) for pair in list(layout.values())]

    return pos


class PortGlyph:
    def __init__(self, ax, string, pos, dir, text_dict):

        from matplotlib.text import Annotation
        self.width = 0.1
        self.height = 0.1

        self.text = Annotation(
            string, pos, **text_dict
        )
        ax.add_artist(self.text)

        self.x, self.y = pos
        if dir == 'top':
            self.y += self.height / 2
        elif dir == 'bottom':
            self.y -= self.height / 2
        elif dir == 'right':
            self.x += self.width / 2
        else:
            self.x -= self.width / 2

    @property
    def pos(self):
        return self.x, self.y


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
            'top': [],
            'right': [],
            'bottom': [],
            'left': []
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
        from matplotlib.text import Text
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

        dx, dy = dir
        text_dict = {
            'size': FONT_SM
        }

        if dy > abs(dx):
            text_dict.update({
                'xytext': (self.x, self.y + self.height / 2),
                'horizontalalignment': 'center',
                'verticalalignment': 'bottom'
            }
            )
            dir = 'top'
        elif -dy > abs(dx):

            text_dict.update({
                'xytext': (self.x, self.y - self.height / 2),
                'horizontalalignment': 'center',
                'verticalalignment': 'top'
            }
            )
            dir = 'bottom'

        elif dx >= abs(dy):

            text_dict.update({
                'xytext': (self.x + self.width / 2, self.y),
                'horizontalalignment': 'left',
                'verticalalignment': 'center'
            }
            )
            dir = 'right'
        else:
            text_dict.update({
                'xytext': (self.x - self.width / 2, self.y),
                'horizontalalignment': 'right',
                'verticalalignment': 'center'
            }
            )
            dir = 'left'

        port = PortGlyph(self.axes, string, self.pos, dir, text_dict)
        self.ports[dir] = port

        return port


class BondView(Line2D):
    th = 3 * np.pi / 4
    R = np.array(((np.cos(th), -np.sin(th)), (np.sin(th), np.cos(th))))
    shortest_bond = None

    def __init__(self, port_1, port_2, *args, **kwargs):
        self.port_1 = port_1
        self.port_2 = port_2
        super().__init__([], [], *args, **kwargs)

    def calc_lines(self):
        x1, y1 = self.port_1.pos
        x2, y2 = self.port_2.pos

        r1 = max(self.port_1.height, self.port_1.width)
        r2 = max(self.port_2.height, self.port_2.width)

        dx = x2 - x1
        dy = y2 - y1
        x1 += r1 * dx
        y1 += r1 * dy
        x2 -= r2 * dx
        y2 -= r2 * dy

        lx, ly = x2 - x1, y2 - y1
        L = np.sqrt((lx)**2 + (ly)**2)

        if not self.shortest_bond or self.shortest_bond > L:
            self.shortest_bond = L

        lx /= L
        ly /= L

        headlength = self.shortest_bond / 5

        vect = np.array((lx, ly))

        assert abs(np.linalg.norm(vect) - 1) < 0.01
        x3, y3 = headlength * self.R.dot(vect) + (x2, y2)

        self.set_xdata([x1, x2, x3])
        self.set_ydata([y1, y2, y3])


class GraphLayout(Glyph):
    def draw(self, ax, layout=_networkx_layout):

        graph = _build_graph(self._node)

        points = layout(graph)
        bonds = []
        x_min = 0
        x_max = 0
        y_min = 0
        y_max = 0
        ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)

        for component, (x, y) in zip(self._node.components,
                                     points):
            x_min = min(x, x_min)
            x_max = max(x, x_max)
            y_min = min(y, y_min)
            y_max = max(y, y_max)

            component.view.pos = (x, y)
            if component.metamodel not in {'0', '1'}:
                try:
                    component.view.string = r"\mathbf{{{t}}}: {n}".format(
                        t=component.metamodel, n=component.name)
                except BaseException:
                    component.view.string = "{t}: {n}".format(
                        t=component.metamodel, n=component.name)
            else:
                try:
                    component.view.string = r"\mathbf{{{t}}}".format(
                        t=component.metamodel)
                except BaseException:
                    component.view.string = "{t}".format(
                        t=component.metamodel)
            component.view.axes = ax

        for tail, head in self._node.bonds:

            c1 = tail.component
            c2 = head.component

            try:
                label_1 = f"[{tail.name}]"
            except AttributeError:
                label_1 = ""

            try:
                label_2 = f"[{head.name}]"
            except AttributeError:
                label_2 = ""

            dx = c2.view.x - c1.view.x
            dy = c2.view.y - c1.view.y

            p1 = c1.view.add_port(label_1, (dx, dy))
            p2 = c2.view.add_port(label_2, (-dx, -dy))
            bond = BondView(p1, p2)
            ax.add_artist(bond)
            bonds.append(bond)

        for bond in bonds:
            bond.calc_lines()

        width = abs(x_max - x_min)
        height = abs(y_min - y_max)
        tweak = 0.1
        ax.axis([x_min - tweak * width,
                 x_max + tweak * width,
                 y_min - tweak * height,
                 y_max + tweak * height])


def find_renderer(fig):

    if hasattr(fig.canvas, "get_renderer"):
        # Some backends, such as TkAgg, have the get_renderer method, which
        # makes this easy.
        renderer = fig.canvas.get_renderer()
    else:
        # Other backends do not have the get_renderer method, so we have a work
        # around to find the renderer.  Print the figure to a temporary file
        # object, and then grab the renderer that was used.
        # (I stole this trick from the matplotlib backend_bases.py
        # print_figure() method.)
        import io
        fig.canvas.print_pdf(io.BytesIO())
        renderer = fig._cachedRenderer
    return(renderer)
