import copy

from .component_manager import get_component, base_id
from .exceptions import *
from .base import BondGraphBase, Bond, Port, FixedPort, PortExpander


def disconnect(target, other):
    """
    Args:
        target: BondGraphComponent or (Port or tuple)
        other: BondGraphComponent or (Port or tuple)
    """
    if isinstance(target, BondGraphBase):
        model = target.parent
    elif isinstance(target, Port):
        model =  target.component.parent
    else:
        model = target[0].parent

    if isinstance(other, BondGraphBase):
        model_prime = other.parent
    elif isinstance(other, Port):
        model_prime = other.component.parent
    else:
        model_prime = other[0].parent

    if not model or not model_prime or (model is not model_prime):
        raise InvalidComponentException(f"Could not find components")

    def _filter(item):
        if isinstance(item, BondGraphBase):
            filtered = {bond for bond in model._bonds if
                        bond.tail.component is item
                        or bond.head.component is item}
        elif isinstance(target, (Bond, tuple)):
            filtered = {bond for bond in model._bonds if bond.tail is item
                        or bond.head is item}
        else:
            raise InvalidComponentException(
                f"{item} is not a valid target for disconnection"
            )
        return filtered

    targets = _filter(target) & _filter(other)

    for bond in targets:
        model._bonds.remove(bond)


def connect(source, destination):
    """
    Connects two components or ports with a bond.
    We assume that either the source and/or destination is or has a free
    port.

    raises:
        InvalidPortException, InvalidComponentException
    """

    if isinstance(source, BondGraphBase):
        model = source.parent
    elif isinstance(source, Port):
        model =  source.component.parent
    else:
        model = source[0].parent

    if isinstance(destination, BondGraphBase):
        model_prime = destination.parent
    elif isinstance(destination, Port):
        model_prime = destination.component.parent
    else:
        model_prime = destination[0].parent

    if not model or not model_prime or (model is not model_prime):
        raise InvalidComponentException(f"Could not find components")

    def find_port(component):
        if isinstance(component, BondGraphBase):
            port = component.get_port()
        elif isinstance(component, Port):
            port = component.component.get_port(component)
        elif isinstance(component, (tuple, list)):
            item, port_idx = component
            port = item.get_port(port_idx)
        else:
            raise InvalidComponentException()
        return port

    tail = find_port(source)
    head = find_port(destination)
    bond = Bond(tail, head)
    model._bonds.add(bond)


def swap(old_component, new_component):
    """
    Replaces the old component with a new component.
    Components must be of compatible classes; 1 one port cannot replace an
    n-port, for example.
    The old component will be completely removed from the system model.

    Args:
        old_component: The component to be replaced. Must already be in the
         model.
        new_component: The substitute component which must not be in the
         model
    """

    # TODO: More validation required
    #
    def is_swap_valid(old_comp, new_comp):
        if not isinstance(new_comp, BondGraphBase):
            return False

        # Dirty Hack because 'in' for BondSet double counts
        num_bonds = len({b for b in old_comp.parent.bonds if old_comp in b})/2

        if isinstance(new_comp, FixedPort) and len(new_comp.ports) < num_bonds:
            return False

        return True

    model = old_component.parent
    if not model:
        return

    if new_component in model.components:
        raise InvalidComponentException(
            "Component is already in the model"
        )

    elif not is_swap_valid(old_component, new_component):
        raise InvalidComponentException("Cannot swap components")

    model.add(new_component)

    swaps = []
    for bond in model.bonds:
        try:
            if bond.tail.component is old_component:
                tail = new_component.get_port()
                head = bond.head
            elif bond.head.component is old_component:
                tail = bond.tail
                head = new_component.get_port()
            else:
                continue
        except InvalidPortException:
            raise InvalidComponentException(
                "Cannot swap components: Incompatible ports"
            )

        swaps.append((bond, Bond(tail, head)))

    for old_bond, new_bond in swaps:
        disconnect(old_bond.tail, old_bond.head)
        connect(new_bond.tail, new_bond.head)

    model.remove(old_component)


def new(component=None, name=None, library=base_id, value=None, **kwargs):
    """
    Creates a new Bond Graph from a library component.

    Args:
        component(str or obj): The type of component to create.
         If a string is specified, the the component will be created from the
         appropriate libaray. If an existing bond graph is given, the bond
         graph will be cloned.
        name (str): The name for the new component
        library (str): The library from which to find this component (if
        component is specified by string).
        value:

    Returns: instance of :obj:`BondGraph`

    """
    if not component:
        cls = _find_subclass("BondGraph", BondGraphBase)
        return cls(name=name)
    elif isinstance(component, str):
        build_args = get_component(component, library)

        if name:
            build_args.update({"name": name})
        if value or isinstance(value, (int, float, complex)):
            _update_build_params(build_args, value, **kwargs)
        cls =_find_subclass(
            build_args["class"], BondGraphBase
        )
        del build_args["class"]

        return cls(**build_args)

    elif isinstance(component, BondGraphBase):
        obj = copy.copy(component)
        if name:
            obj.name = name
        if value:
            _update_build_params(obj.__dict__, value)

        return obj

    else:
        raise NotImplementedError(
            "New not implemented for object {}", component
        )


def _update_build_params(build_args, value, **kwargs):

    if isinstance(value, (list, tuple)):
        assignments = zip(build_args["params"].keys(), value)
        for param, v in assignments:
            build_args["params"][param]["value"] = v
    elif isinstance(value, dict):
        for param, v in value.items():
            if isinstance(build_args["params"][param], dict):
                build_args["params"][param]["value"] = v
            else:
                build_args["params"][param] = v
    else:
        p = next(iter(build_args["params"]))
        build_args["params"][p] = value


def _find_subclass(name, base_class):

    for c in base_class.__subclasses__():
        if c.__name__ == name:
            return c
        else:
            sc = _find_subclass(name, c)
            if sc:
                return sc


