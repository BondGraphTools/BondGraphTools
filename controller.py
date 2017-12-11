import component_manager as cm


class Controller(object):
    add_component = 0
    add_bond = 1

    def __init__(self, model):
        self.history = list()
        self.undo_history = list()
        self.context = None
        self.bondgraph = model
        self.view = None

    def undo(self):
        pass

    def perform(self, command, pos=None):
        """

        Args:
            command: a tuple
            pos:

        Returns:

        """
        action, params = command
        if action == self.add_component:
            library, component = params
            result = self._add_node(library, component, pos)
        elif action == self.add_bond:
            from_node, to_node = params

        self.history.append(
            (action, params, pos, result)
        )

        self.context = None

    def _add_node(self, library, component, pos, node_id=None):

        new_node = cm.build_node(library, component)
        new_node.pos = pos
        result_id = self.bondgraph.add_node(new_node, node_id)
        self.view.add_component(
            node_id=result_id,
            node_type=new_node.node_type,
            pos=pos
        )
        return result_id

    def _delete_node(self, node_id):
        self.bondgraph.remove_node(node_id)
        self.view.remove_component(node_id)

    def undo_action(self):

        action, params, pos, result = self.history.pop()

        if action == self.add_component:
            self._delete_node(result)
            self.undo_history.append(
                (action, params, pos, result)
            )

    def redo_action(self):

        action, params, pos, result = self.undo_history.pop()

        if action == self.add_component:
            library, component = params
            self._add_node(library, component, pos, result)

        self.history.append(
            (action, params, pos, result)
        )






