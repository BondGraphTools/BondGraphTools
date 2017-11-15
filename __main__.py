
import wx


SIZE = 20

class EditorFrame(wx.Frame):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

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
        self.Bind(wx.EVT_LEFT_DOWN, self.on_click)
        self.Bind(wx.EVT_PAINT, self.on_paint)
        self.nodes = []

    def on_click(self, event):

        pos = event.GetPosition()
        nearest_node = None
        for node in self.nodes:
            distance = _pabs(node.GetPosition() - pos)
            if distance < SIZE:
                nearest_node = node
                return


        new_node = Node(pos, (SIZE, SIZE))
        self.nodes.append(new_node)
        self.Parent.Refresh()




    def on_paint(self, event):
        dc = wx.PaintDC(self)

        for node in self.nodes:
            dc.SetPen(wx.Pen("black"))
            dc.SetBrush(wx.Brush(wx.TRANSPARENT_BRUSH))
            dc.DrawRectangle(node.GetPosition(), node.GetSize())


class Node(wx.Rect):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)





if __name__ == '__main__':
    app = wx.App()

    frame = EditorFrame(None, title='Bond Graph Editor')

    frame.Show()

    app.MainLoop()