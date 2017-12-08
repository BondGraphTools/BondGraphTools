import math

from PyQt5.QtWidgets import (QMainWindow, QWidget, QDockWidget,
                             QGridLayout, QPushButton)
from PyQt5.QtCore import Qt
import component_manager as cm

add_component = "add_component"

B_SZ = 25
PAD = 15

class EditorCanvas(QWidget):
    def __init__(self, model, controller):
        super().__init__()
        self.model = model
        self.controller = controller

    def mousePressEvent(self, event):
        command = self.controller.context
        try:
            if event.button() == Qt.LeftButton and command:
                pos = event.pos()
                print((command, pos))
        except AttributeError:
            pass
class Controller(object):
    add_component = 0
    add_bond = 1

    def __init__(self, model):
        self.history = list()
        self.undo_history = list()
        self.context = None
    def undo(self):
        pass

    def set_context(self, context):
        self.context = context


class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        model = None # blank model
        controller = Controller(model)

        view = EditorCanvas(model, controller)
        self.setBaseSize(640, 480)
        self.resize(640, 480)
        view.show()
        widget = BaseComponentPalette(controller)
        base_docker = QDockWidget("Base", self)
        base_docker.setWidget(widget)

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
            new_button.command = (controller.add_component, cm.base_id, comp_id)
            new_button.clicked.connect(self.button_clicked)

            grid.addWidget(new_button, int(n/3), n % 3)

        rows = math.ceil(len(components)/3)
        height = rows*B_SZ + PAD*rows
        width = 3*B_SZ + 2 * PAD

        self.setFixedSize(width, height)


    def button_clicked(self):

        sender = self.sender()
        self.controller.set_context(sender.command)