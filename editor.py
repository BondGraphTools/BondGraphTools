import math

from PyQt5.QtWidgets import (QMainWindow, QWidget, QDockWidget,
                             QGridLayout, QPushButton, QGraphicsItem,
                             QGraphicsScene, QGraphicsTextItem, QGraphicsView,
                             QGraphicsPathItem)
from PyQt5.QtCore import Qt, QPointF
from PyQt5.QtGui import QPalette, QTransform
import component_manager as cm
import controller as co
import model as mo
B_SZ = 25
PAD = 15
TEXT_SIZE = 16

class EditorCanvas(QGraphicsScene):

    def __init__(self, model, controller, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.model = model
        self.controller = controller
        self._components = dict()


    def mousePressEvent(self, event):
        command = self.controller.context

        l_button = event.button() == Qt.LeftButton
        pos = event.scenePos()
        x = pos.x()
        y = pos.y()
        if l_button and command:
            self.controller.perform(
                command=command,
                pos=(x, y)
            )
        elif (l_button and event.modifiers() == Qt.ShiftModifier and
            self.selectedItems()):

            n0 = self.selectedItems()[0].node_id
            n1 = self.itemAt(x, y, QTransform())
            bond_id = self.controller.perform(
                (self.controller.add_bond, n0, n1)
            )
        else:
            super().mousePressEvent(event)

    def add_component(self, node_id, node_type, pos):
        c = ComponentWidget(
                node_type=node_type, node_name=node_id,
                node_id=node_id)
        self.addItem(c)
        x, y = pos
        c.setPos(x, y)
        self.update()


class BondWidget(QGraphicsPathItem):

    def __init__(self, start_node, end_node):
        super().__init__()
        self.start_node = start_node
        self.end_node = end_node
        self.update()

    def update(self, *__args):

        r2 = self.end_node.sceneBoundingRect()
        r1 = self.start_node.sceneBoundingRect()

        v = r2.center() - r1.center()
        modv2 = v.dotProduct(v)
        super().update(*__args)


class ComponentWidget(QGraphicsTextItem):
    def __init__(self, node_id, node_type, node_name):
        super().__init__(
            "{}: {}".format(node_type, node_name)
        )
        self.setHtml(
            "<b>{}:</b> {}".format(node_type, node_name)
        )
        self.node_type = node_type
        self.name = node_name
        self.node_id = node_id
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)


        # sp = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        # self.setSizePolicy(sp)
    #
    # def paintEvent(self, event):
    #
    #     qp = QPainter()
    #     qp.begin(self)
    #     self.paint(qp)
    #     qp.end()
    #
    # def paint(self, painter):
    #
    #     painter.setPen(Qt.black)
    #     painter.setRenderHint(QPainter.Antialiasing)
    #     print("painting")
    #     #bld_font = QFont("Computer Modern Roman", TEXT_SIZE, QFont.Bold)
    #     #base_font = QFont("Computer Modern Roman", TEXT_SIZE)
    #     bld_font = QFont("Arial", TEXT_SIZE, QFont.Bold)
    #     base_font = QFont("Arial", TEXT_SIZE)
    #     x, y = self.pos
    #     painter.drawPoint(x, y)
    #     height = TEXT_SIZE + 4
    #
    #     width = len(self.node_type)
    #     rect = QRect(x - width, y - int(height/2), width, height)
    #     painter.setFont(bld_font)
    #     painter.drawText(
    #         rect,
    #         Qt.AlignRight,
    #         "{}: ".format(self.node_type)
    #     )
    #
    #     rect = QRect(x, y - int(height / 2), width, height)
    #     painter.setFont(base_font)
    #
    #     painter.drawText(
    #         rect,
    #         Qt.AlignRight,
    #         "{}: ".format(self.name)
    #     )


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        model = mo.BondGraph()  # blank model
        controller = co.Controller(model)

        self.setWindowTitle("Bond Graph Editor")

        self.setBaseSize(640, 480)
        self.resize(640, 480)

        widget = BaseComponentPalette(controller)
        base_docker = QDockWidget("Base Library", self)
        base_docker.setWidget(widget)

        #area = QScrollArea()
        scene = EditorCanvas(model, controller)
        controller.view = scene
        view = QGraphicsView(scene)
        view.setMinimumSize(640, 480)


        p = QPalette()
        p.setColor(view.backgroundRole(), Qt.white)
        view.setPalette(p)
        view.setAutoFillBackground(True)

        self.setCentralWidget(view)
        base_docker.setAllowedAreas(Qt.AllDockWidgetAreas)
        self.addDockWidget(Qt.LeftDockWidgetArea, base_docker, Qt.Vertical)


class BaseComponentPalette(QWidget):
    def __init__(self, controller):

        super().__init__()
        self.controller = controller
        grid = QGridLayout()
        self.setLayout(grid)

        components = cm.get_components_list(cm.base_id)
        for n, (comp_id, name) in enumerate(components):
            new_button = QPushButton(comp_id, self)
            new_button.setFixedSize(B_SZ, B_SZ)
            new_button.setToolTip(name)
            new_button.command = (controller.add_component,
                                  (cm.base_id, comp_id))
            new_button.clicked.connect(self.button_clicked)

            grid.addWidget(new_button, int(n/3), n % 3)

        rows = math.ceil(len(components)/3)
        height = rows*B_SZ + PAD*rows
        width = 3*B_SZ + 2 * PAD

        self.setFixedSize(width, height)


    def button_clicked(self):

        sender = self.sender()
        self.controller.context =  sender.command