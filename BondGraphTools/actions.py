from .exceptions import *
from .base import BondGraphBase, Bond, Port


def disconnect(target, other):
    """
    Args:
        target: BondGraphComponent or (Port or tuple)
        other: BondGraphComponent or (Port or tuple)
    """
    if isinstance(target, BondGraphBase):
        model = target.parent
    elif isinstance(target,Port):
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

    def _find_port(component):
        if component not in model.components:
            raise InvalidComponentException(
                "Could not find %s: not contained in %s", component, model)
        used_ports = {p for bond in model.bonds for (c, p) in bond
                      if c is component}

        free_ports = set(component.ports) - used_ports

        if not free_ports:
            try:
                port = component.make_port()
            except AttributeError:
                raise InvalidPortException(
                    "Could not find a free port on %s", component)
        elif len(free_ports) > 1:
            raise InvalidPortException(
                "Could not find a unique free port on %s: "
                "specify a port ", component)
        else:
            port = free_ports.pop()
        return component, port


    def _validate_port(component, port):
        if component not in model.components:
            raise InvalidComponentException(f"Component {component} "
                                            f"not found ")
        elif (component, port) in {p for bond in model.bonds for p in bond}:
            raise InvalidPortException("Could not connect port: already in"
                                       "use")
        if port not in component.ports:
            component.make_port(port)

    if isinstance(source, BondGraphBase):
        src, src_port = _find_port(source)

    elif isinstance(source, (tuple, list)):
        src, src_port = source
        _validate_port(src, src_port)
    else:
        raise InvalidComponentException()

    if isinstance(destination, BondGraphBase):
        dest, dest_port = _find_port(destination)

    elif isinstance(destination, (tuple, list)):
        dest, dest_port = destination
        _validate_port(dest, dest_port)
    else:
        raise InvalidComponentException()

    bond = Bond(Port(src, src_port),Port(dest, dest_port))

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

    #TODO: More validation required
    def is_swap_valid(old_comp, new_comp):
        if not isinstance(new_comp, BondGraphBase):
            return False
        elif new_comp.max_ports and new_comp.max_ports < len(old_comp.ports):
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

    for bond in model.bonds:
        if bond.tail.component is old_component:
            tail = Port(component=new_component, port=bond.tail.port)
            head = bond.head
        elif bond.head.component is old_component:
            tail = bond.tail
            head = Port(component=new_component, port=bond.head.port)
        else:
            continue
        disconnect(*bond)
        connect(tail, head)

    model.remove(old_component)




