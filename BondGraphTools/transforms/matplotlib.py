"""
Matplotlib Interface
====================

make sure to set up your matplotlib interface prior to importing this module

You may also change the plot style by changing 
matplotlib.pyplot.stye  
in the usual manner
"""

import os
import logging
import numpy as np

import matplotlib.pyplot as pyplot

logger = logging.getLogger(__name__)

_mplstyle = os.path.join(os.path.dirname(__file__), "mtt2.mplstyle")
pyplot.style.use(_mplstyle)

phi = (5**0.5 - 1)/2


def draw(graph):
    """
    Uses Matplotlib to draw the bond graph

    Args:
        graph: The bond graph

    Returns:
        :obj:matplotlib.pyplot.figure
    """

    try:
        X = [node.pos[0] for node in graph.nodes.values()]
        Y = [node.pos[1] for node in graph.nodes.values()]
        min_x = min(X)
        max_x = max(X)
        min_y = min(Y)
        max_y = max(Y)
    except ValueError:
        logger.warning("Cannot draw empty graph")
        return None

    fig = pyplot.figure()

    width = abs(max_x - min_x)
    height = abs(max_y - min_y)
    height = max(height, width/phi)

    tweak_x = 0.25*width
    tweak_y = 0.25*height

    pyplot.axis([min_x - tweak_x,
                 max_x + tweak_x,
                 min_y - tweak_y,
                 max_y + tweak_y])
    renderer = find_renderer(fig)
    ax = pyplot.gca()
    eps = []

    for node_id, node in graph.nodes.items():
        x, y = node.pos
        text = pyplot.text(x, y, str(node),
                           horizontalalignment='center',
                           verticalalignment='center',
                           bbox=dict(facecolor='none', edgecolor='black')
        )

        bbox = text.get_window_extent(renderer).transformed(
            ax.transData.inverted())
        ep = np.sqrt(bbox.width**2 + bbox.height**2)
        eps.append(ep)

    for i, j, p_i, p_j in graph.bonds:
        x1, y1 = graph.nodes[i].pos
        x2, y2 = graph.nodes[j].pos

        vx1 = x2 - x1
        vy1 = y2 - y1
        vx1 = vx1 / np.sqrt(vx1**2 + vy1**2)
        vy1 = vy1 / np.sqrt(vx1**2 + vy1**2)
        vx2 = (vy1 - vx1)/np.sqrt(2)
        vy2 = -(vy1 + vx1)/np.sqrt(2)

        # eps1 = eps[i]*1.25
        # eps2 = eps[j]*1.25
        # arrow_l = 0.1
        # line_x = [x1 + eps1 * vx1, x2 - eps2 * vx1,
        #           x2 - eps2 * vx1 + arrow_l * vx2]
        # line_y = [y1 + eps1 * vy1, y2 - eps2 * vy1,
        #           y2 - eps2 * vy1 + arrow_l * vy2]

        eps1 = eps[i]
        eps2 = eps[j]
        line_x = [x1, x2]
        line_y = [y1, y2]
        pyplot.plot(line_x, line_y, color='k')

    fig.tight_layout(pad=0.2)
    fig.suptitle(graph.name)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')
    return fig


def find_renderer(fig):
    """
    source: https://stackoverflow.com/questions/22667224/matplotlib-get-text-bounding-box-independent-of-backend#22689498
    Args:
        fig:

    Returns:

    """
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

    return renderer
