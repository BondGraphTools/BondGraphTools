import math
from sys import platform as _platform

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtPrintSupport import *
import component_manager as cm
import controller as co
import model as mo

_OSX = "darwin"


class config(object):
    font_small = 10
    font_size = 16
    grid_size = 10
    grid_line_weight = 0.25
    base_view_size = (640, 480)
    class toolbar(object):
        pad = 14
        button_size = 25


def snap(x):
    return round(x/config.grid_size)*config.grid_size


class EditorCanvas(QGraphicsScene):
    shift_down = pyqtSignal(bool)

    def __init__(self, model, controller, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.model = model
        self.controller = controller
        self._components = dict()
        self.show_grid = True

    def keyPressEvent(self, event):

        super().keyPressEvent(event)

        if event.key() == Qt.Key_Shift:
            self.shift_down.emit(True)

    def keyReleaseEvent(self, event):
        super().keyReleaseEvent(event)
        if event.key() == Qt.Key_Shift:
            self.shift_down.emit(False)

    def mousePressEvent(self, event):
        command = self.controller.context

        l_button = event.button() == Qt.LeftButton
        pos = event.scenePos()
        x = pos.x()
        y = pos.y()
        if l_button and command:
            self.controller.perform(
                command=command,
                pos=(snap(x), snap(y))
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

    def mouseReleaseEvent(self, event):

        super().mouseReleaseEvent(event)

    def add_component(self, node_id, node_type, pos):
        if node_type == "0":
            c = ZeroWidget(node_id=node_id)
        else:
            c = ComponentWidget(node_type=node_type, node_id=node_id)
        self.addItem(c)
        self.shift_down.connect(c.shift)
        x, y = pos
        c.setPos(x, y)

        self.update()

    def drawForeground(self, painter, rect):
        super().drawForeground(painter, rect)

    def drawBackground(self, painter, rect):
        super().drawBackground(painter, rect)
        if self.show_grid:
            t = rect.top()
            b = rect.bottom()
            l = rect.left()
            r = rect.right()

            y = int(t / config.grid_size) * config.grid_size
            x = int(l / config.grid_size) * config.grid_size
            lines = []

            while y <= b:
                lines.append(QLineF(l, y, r, y))
                y += config.grid_size

            while x <= r:
                lines.append(QLineF(x,b,x,t))
                x += config.grid_size

            pen = QPen(Qt.DashLine)
            pen.setColor(Qt.gray)
            pen.setWidthF(config.grid_line_weight)

            painter.setPen(pen)
            painter.drawLines(lines)


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


class ComponentWidget(QGraphicsWidget):
    def __init__(self, node_id, node_type, ports=None, **kwargs):

        super().__init__(**kwargs)
        self.node_type = node_type
        self.node_id = node_id
        self.setFlag(QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)
        self.ports = ports
        self.center = ComponentTextWidget(
            "{}: {}".format(node_type, node_id))
        layout = PortLayoutManager(self.center)
        self.setLayout(layout)
        self.setFocusProxy(self.center)

    def itemChange(self, change_type, value):

        if change_type == self.ItemPositionChange:
            value.setX(snap(value.x()))
            value.setY(snap(value.y()))

        return super().itemChange(change_type, value)

    @pyqtSlot(bool)
    def shift(self, is_down):
        if is_down:
            for child in self.childItems():
                if child is not self.center:
                    child.show()
            self.setFlag(QGraphicsItem.ItemIsMovable, False)
            self.setFlag(QGraphicsItem.ItemIsSelectable, False)
        else:
            self.setFlag(QGraphicsItem.ItemIsMovable, True)
            self.setFlag(QGraphicsItem.ItemIsSelectable, True)
            for child in self.childItems():
                if child is not self.center:
                    child.hide()


class PortLayoutManager(QGraphicsAnchorLayout):
    layout_sequence = {
        0: [(Qt.AnchorTop, Qt.AnchorBottom)],
        1: [(Qt.AnchorRight, Qt.AnchorLeft)],
        2: [(Qt.AnchorHorizontalCenter, Qt.AnchorRight)],
        3: [(Qt.AnchorHorizontalCenter, Qt.AnchorLeft)]}

    def __init__(self, center, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.center = center
        self.addAnchor(self, Qt.AnchorVerticalCenter,
                       center, Qt.AnchorVerticalCenter)
        self.addAnchor(self, Qt.AnchorHorizontalCenter,
                       center, Qt.AnchorHorizontalCenter)
        self.items = 0

    def add_port(self, widget):
        if self.items == 0:
            self.addAnchor(self.center, Qt.AnchorTop,
                           widget, Qt.AnchorVerticalCenter)
            self.addAnchor(self.center, Qt.AnchorLeft, widget, Qt.AnchorRight)
        elif self.items == 1:
            self.addAnchor(self.center, Qt.AnchorTop,
                           widget, Qt.AnchorVerticalCenter)
            self.addAnchor(self.center, Qt.AnchorRight, widget, Qt.AnchorLeft)
        else:
            raise NotImplementedError
        self.items += 1

    def remove_port(self, n=None):

        pass


class ComponentTextWidget(QGraphicsWidget):
    def __init__(self, text, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.text = text
        self.font = QFont("Computer Modern Roman", config.font_size)
        self.metrics = QFontMetricsF(self.font)

    def sizeHint(self, *args, **kwargs):
        return self.metrics.size(Qt.AlignCenter | Qt.TextSingleLine, self.text)

    def paint(self, painter, style, widget=None):
        painter.setFont(self.font)
        painter.setPen(Qt.black)
        painter.drawText(self.boundingRect().adjusted(0, -2, 0, -1),
                         Qt.AlignCenter, self.text)
        painter.drawRect(self.boundingRect().adjusted(0, -2, 0, -1))


class ZeroWidget(ComponentWidget):
    def __init__(self, node_id):
        super().__init__(node_id, node_type="0")
        self.center.text = "0"
        self.n_ports = 0
        self.add_port()

    def add_port(self, show=False):
        w = PortWidget(port_id=self.n_ports, label=" ")
        w.setParent(self)
        w.setVisible(show)
        self.layout().add_port(w)
        self.n_ports += 1


class PortWidget(QGraphicsWidget):
    def __init__(self, port_id, label, domain=None, restrictions=None, **kwargs):

        super().__init__(**kwargs)
        self.port_id = port_id
        self.text = "[{}]".format(label)
        self.font = QFont("Computer Modern Roman", config.font_small)
        self.metrics = QFontMetricsF(self.font)

    def sizeHint(self, *args, **kwargs):
        return self.metrics.size(Qt.TextSingleLine, self.text)

    def paint(self, painter, style, widget=None):
        painter.setFont(self.font)
        painter.setPen(Qt.red)
        painter.drawText(self.boundingRect().adjusted(0, -1, 0, 0),
                         Qt.AlignCenter, self.text)
        painter.drawRect(self.boundingRect())

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            print("node {} port {} clicked".format(self.parent().node_id,
                                                   self.port_id))


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        model = mo.BondGraph()  # blank model
        self.controller = co.Controller(model)

        w, h = config.base_view_size
        self.setWindowTitle("Bond Graph Editor")

        self.setBaseSize(w, h)
        self.resize(w, h)

        self.build_menu()

        widget = BaseComponentPalette(self.controller)
        base_docker = QDockWidget("Base Library", self)
        base_docker.setWidget(widget)

        # area = QScrollArea()
        scene = EditorCanvas(model, self.controller)
        self.controller.view = scene
        scene.setSceneRect(-w/2, -h/2, w/2, h/2)
        view = QGraphicsView(scene)
        view.setMinimumSize(w, h)
        view.setResizeAnchor(QGraphicsView.AnchorViewCenter)
        view.setTransformationAnchor(QGraphicsView.AnchorViewCenter)
        p = QPalette()
        p.setColor(view.backgroundRole(), Qt.white)
        view.setPalette(p)
        view.setAutoFillBackground(True)

        self.setCentralWidget(view)
        base_docker.setAllowedAreas(Qt.AllDockWidgetAreas)
        base_docker.setFeatures(QDockWidget.DockWidgetFloatable |
                                QDockWidget.DockWidgetMovable |
                                QDockWidget.DockWidgetVerticalTitleBar)
        self.addDockWidget(Qt.LeftDockWidgetArea, base_docker, Qt.Vertical)

    def build_menu(self):
        if _platform == _OSX:
            bar = self.menuBar()
            bar.setNativeMenuBar(False)

        self._build_file_menu()
        self._build_edit_menu()
        # edit_menu = bar.addMenu("Edit")
        # tool_menu = bar.addMenu("Tools")
        #
        # help_menu = bar.addMenu("Help")

    def _build_file_menu(self):

        bar = self.menuBar()

        file = bar.addMenu("&File")
        new = QAction("&New", self)
        open_file = QAction("&Open", self)
        save = QAction("&Save", self)
        saveas = QAction("Save &As", self)
        print_action = QAction("&Print", self)
        print_action.triggered.connect(self.menu_print)
        file_quit = QAction("&Quit", self)
        file_quit.triggered.connect(qApp.quit)

        file.addAction(new)
        file.addAction(open_file)
        file.addAction(save)
        file.addAction(saveas)
        file.addAction(print_action)
        file.addAction(file_quit)

    def menu_print(self):
        printer = QPrinter()
        scene = self.centralWidget().scene()

        if QPrintDialog(printer).exec() == QDialog.Accepted:
            has_grid = scene.show_grid
            scene.show_grid = False
            painter = QPainter(printer)
            painter.setRenderHint(QPainter.Antialiasing)
            scene.render(painter)
            painter.end()
            scene.show_grid = has_grid

    def _build_edit_menu(self):
        bar = self.menuBar()

        undo = QAction("&Undo", self)
        redo = QAction("&Redo", self)
        edit = bar.addMenu("&Edit")
        edit.addAction(undo)
        edit.addAction(redo)


    def resizeEvent(self, event):
        """
        todo: resize canvas
        """
        super().resizeEvent(event)
        view = self.centralWidget()
        rect = QRectF(view.rect())
        view_rect = rect.adjusted(-rect.width()/2,
                                  -rect.height()/2,
                                  -rect.width()/2,
                                  -rect.height()/2)

        bounding_rect = view.scene().itemsBoundingRect()
        x1 = min(view_rect.left(), bounding_rect.left())
        y1 = min(view_rect.top(), bounding_rect.top())
        x2 = max(view_rect.right(), bounding_rect.right())
        y2 = max(view_rect.bottom(), bounding_rect.bottom())
        view.setSceneRect(x1, y1, x2-x1, y2-y1)


class BaseComponentPalette(QWidget):
    def __init__(self, controller):

        super().__init__()
        self.controller = controller
        grid = QGridLayout()
        self.setLayout(grid)

        components = cm.get_components_list(cm.base_id)
        for n, (comp_id, name) in enumerate(components):
            new_button = ComponentLibraryWidget(comp_id, self)
            size = config.toolbar.button_size
            new_button.setFixedSize(size, size)
            new_button.setToolTip(name)
            new_button.command = (controller.add_component,
                                  (cm.base_id, comp_id))

            grid.addWidget(new_button, int(n/3), n % 3)

        rows = math.ceil(len(components)/3)
        height = rows*size + config.toolbar.pad*rows
        width = 3*size + 2 * config.toolbar.pad

        self.setFixedSize(width, height)


class ComponentLibraryWidget(QPushButton):
    """
    Todo:
    Implement Drag and Drop functionality
    """
    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        self.parent().controller.context = self.command
