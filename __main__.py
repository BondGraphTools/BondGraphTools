#! python3

import wx
import math

SIZE = 40

class EditorFrame(wx.Frame):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.status_bar = self.CreateStatusBar(1)
        self.panel = Editor(
            parent=self, size=self.Size, pos=(-1,-1)
        )

        self.make_menu_bar()

    def make_menu_bar(self):

        file_menu = wx.Menu()
        new_item = file_menu.Append(wx.ID_NEW, "&New")
        save_item = file_menu.Append(wx.ID_SAVE, "&Save")
        exit_item = file_menu.Append(wx.ID_EXIT, "&Quit")

        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, "&File")

        self.SetMenuBar(menu_bar)

        self.Bind(wx.EVT_MENU, self.on_exit, exit_item)
        self.Bind(wx.EVT_MENU, self.noop, new_item)
        self.Bind(wx.EVT_MENU, self.noop, save_item)

    def noop(self, event):
        pass

    def on_exit(self, event):
        self.Close(True)

def _pabs(point):

    return sum(
        x**2 for x in point.Get()
    )**0.5


class Editor(wx.Panel):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.SetBackgroundColour("white")

        self.nodes = dict()
        self.edges = []

        self._dragged_obj = None
        self._next_obj_cls = None
        self._adding_bonds = False
        self._bond_start = None
        self.Bind(wx.EVT_LEFT_UP, self.on_release)
        self.Bind(wx.EVT_LEFT_DOWN, self.on_click)
        self.Bind(wx.EVT_MOTION, self.on_drag)
        self.Bind(wx.EVT_PAINT, self.on_paint)

        self.Bind(wx.EVT_KEY_UP, self.on_key_up)

    def on_key_up(self, event):
        key = chr(event.GetUnicodeKey())
        if key in component_keymap:
            self._next_obj_cls = component_keymap[key]
            self.Parent.status_bar.SetStatusText(
                "Adding: {0}".format(self._next_obj_cls.desc)
            )
        elif key == " ":
            self.Parent.status_bar.SetStatusText("Adding Bonds")
            self._adding_bonds = True
            self._bond_start = None
        else:
            self._adding_bonds = False
            self.Parent.status_bar.SetStatusText("waiting")
            self._next_obj_cls = None

    def on_release(self, event):
        if self._dragged_obj:
            self._dragged_obj = None
            return

    def _get_node(self, point, distance=SIZE):

        for node_id, node in self.nodes.items():
            dist = _pabs(node.GetPosition() - point)
            if dist < distance:
                return node_id
        return None

    def on_drag(self, event):
        if self._dragged_obj:
            x, y = event.GetPosition().Get()
            self._dragged_obj.SetPosition((x - SIZE/2, y - SIZE/2))
            self.Refresh()

    def on_click(self, event):

        pos = event.GetPosition()
        node_id = self._get_node(pos)

        if node_id:
            self._dragged_obj = self.nodes[node_id]
        elif self._next_obj_cls:
            dialog = wx.TextEntryDialog(
                self, "Node Name", caption="Node Name",
                value="", style=wx.TextEntryDialogStyle, pos=pos)

            if dialog.ShowModal():
                node_id = dialog.GetValue()

                new_node = self._next_obj_cls(
                    (pos.x - SIZE/2, pos.y - SIZE/2), (SIZE, SIZE))
                self._next_obj_cls = None
                self.nodes[node_id] = new_node
                self.Parent.Refresh()
        elif self._adding_bonds:
            if not node_id:
                return
            elif not self._bond_start:
                self._bond_start = node_id




    def on_paint(self, event):
        event.Skip()

        dc = wx.PaintDC(self)

        for node_id, node in self.nodes.items():
            dc.DrawLabel(
                node.label.format(node_id=node_id),
                rect=node,
                alignment=wx.ALIGN_CENTER
            )

        for start_node, end_node, _ in self.edges:
            line_start, line_end = get_line(start_node, end_node)

            v = (y - x for x, y in zip(line_start, line_end))
            mod_v = sum(x**2 for x in v)**0.5
            v = rotate(v, 135)
            arrow_end = (ARROW*(x + y/mod_v) for x, y in zip(line_end, v))

            dc.DrawLine(line_start, line_end)
            dc.DrawLine(line_end, arrow_end)

ARROW= 20


def rotate(tup, angle_d):
    sin_d = math.sin(angle_d)
    cos_d = math.cos(angle_d)

    x = tup[0]*cos_d - tup[1]*sin_d
    y = tup[0]*sin_d + tup[1]*cos_d
    return (x, y)

def get_line(start_rect, end_rect):
    diff_y = end_rect.GetPosition().y - start_rect.GetPosition().y
    diff_x = end_rect.GetPosition().x - start_rect.GetPosition().x

    if abs(diff_x) > abs(diff_y): # ie latch to horizonta;s.
        if diff_x > 0:
            start = start_rect.GetRight()
            end = end_rect.GetLeft()
        else:
            start = start_rect.GetLeft()
            end = end_rect.GetRight()
    else:
        if diff_y > 0:
            start = start_rect.GetBottom()
            end = end_rect.GetTop()
        else:
            start = start_rect.GetTop()
            end = end_rect.GetBottom()

    return start, end






class Node(wx.Rect):
    label = "?"
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CCompenent(Node):
    ports = (1, 0)
    desc = "Accumulator Component"
    label = "C: {node_id}"

class RComponent(Node):
    ports = (1, 0)
    desc = "Resistive Component"
    label = "R: {node_id}"

class IComponent(Node):
    desc = "Intertial Storage Component"
    ports = (1, 0)
    label = "I: {node_id}"

class ZeroComponent(Node):
    desc = "Common Effort Connector"
    ports = None #None is shorthand for any number
    label = "0"

class OneComponent(Node):
    desc = "Common Flow Connector"
    ports = None
    label = "1"

class TFComponent(Node):
    desc = "Transformer Component"
    ports = (1, 1)
    label = "TF: {node_id}"

class GYComponent(Node):
    desc = "Gyrator Component"
    ports = (1,1)
    label = "GY: {node_id}"

class SeComponent(Node):
    desc = "Effort Source"
    ports = (1, 0)
    label = "Se: {node_id}"

class SfComponent(Node):
    desc = "Flow Source"
    ports = (1, 0)
    label = "Sf: {node_id}"

class OeComponent(Node):
    desc = "Effort Sink"
    ports = (0, 1)
    label = "Se: {node_id}"

class OfComponent(Node):
    desc = "Flow Sink"
    ports = (0, 1)
    label = "Sf: {node_id}"


component_keymap = {
    "C": CCompenent,
    "I": IComponent,
    "T": TFComponent,
    "0": ZeroComponent,
    "1": OneComponent,
    "R": RComponent,
    "E": SeComponent,
    "F": SfComponent,
    "G": GYComponent
}

if __name__ == '__main__':
    app = wx.App()

    frame = EditorFrame(None, title='Bond Graph Editor')

    frame.Show()

    app.MainLoop()